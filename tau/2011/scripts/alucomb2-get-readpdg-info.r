#!/usr/bin/env Rscript

##
## 
##

require(stringr, quietly=TRUE)

## ////////////////////////////////////////////////////////////////////////////
## definitions

##--- substitute after evaluating arg
esub = function(expr, sublist) do.call("substitute", list(expr, sublist))
esub.expr = function(expr, sublist) {
  sapply(expr, function(call) as.expression(esub(call, sublist)))
}

##--- evaluates expression recursively substituting defs
proc = function(e, env = parent.frame()) {
  for(nm in all.vars(e)) {
    if (exists(nm, env) && is.language(g <- get(nm, env))) {
      if (is.expression(g)) g <- g[[1]]
      g <- Recall(g, env)
      L <- list(g)
      names(L) <- nm
      e <- esub(e, L)
    }
  }
  e
}

##--- return just the matches, string can be a vector
str.match.just = function(string, pattern) {
  rc = str_match(string, pattern)
  return(rc[, -1, drop=FALSE])
}

##--- return just the matches, assume string is of length one
str.match.just.single = function(string, pattern) {
  rc = str_match(string, pattern)
  return(rc[1, -1, drop=TRUE])
}

##
## return numeric id for sorting labels like "Gamma5", "Gamma3by5", ...
##
alucomb2.gamma.num.id = function(gamma.name) {
  gamma.name = ifelse(gamma.name == "GammaAll", "Gamma999", gamma.name)
  num1 = as.numeric(sub("^(\\D+)(\\d+)(|by.*)$", "\\2", gamma.name, perl=TRUE))
  tmp = sub("^\\D+\\d+(by|)(\\d*)$", "\\2", gamma.name)
  num2 = ifelse(tmp=="", 0, as.numeric(tmp))
  return(10000*num1 + num2)
}

##
## print matrix with even formatting
##
alu.rbind.print = function(x, width=13, num.columns=NULL) {
  number.width = width
  rn.max = max(nchar(rownames(x)))
  cn.max = max(nchar(colnames(x)))
  cn.max = max(number.width, cn.max)
  fmt = paste("% -", number.width, ".", number.width-7, "g", sep="")
  width = getOption("width")
  digits = getOption("digits")
  if (is.null(num.columns)) {
    items.per.row = max(1, floor((width-rn.max) / (cn.max+1))) 
  } else {
    items.per.row = num.columns
  }
  for (i.first in seq(1, ncol(x), by=items.per.row)) {
    ## print(format(x), quote=FALSE)
    i.last = min(i.first + items.per.row - 1, ncol(x))
    cat(format("", width=rn.max+1), paste(format(colnames(x)[i.first:i.last], width=cn.max), collapse=" "), "\n", sep="")
    mapply(function(label, vec) {
      cat(format(label, width=rn.max), " ", paste(format(sprintf(fmt, unlist(vec)), width=cn.max), collapse=" "), "\n", sep="")
    }, rownames(x), apply(x[,i.first:i.last, drop=FALSE], 1, list))
  }
  return(invisible(NULL))
}

##
## test if pattern matches string irrespective of letters case
##
match.nocase = function(pattern, str) {
  return(regexpr(pattern, str, ignore.case=TRUE) != -1)
}

##/////////////////////////////////////////////////////////////////////////////
##
## - read output file produced by readpdg.cc
##   - quantity description
##   - PDG node
##   - PDG Gamma number
## - replace numbers with expressions with symbolic parameters
##
read.info = function(file) {

  quantities = list()
  constr.lin.val = numeric(0)
  constr.lin.comb = list()
  constr.nl.val = numeric(0)
  constr.nl.expr = list()
  
  if (!file.exists(file)) {
    stop("cannot find file ", file, "\n")
  }
  dir.base = dirname(file)
  lines = readLines(file)
  cat("read file", file, "\n", file=stderr())
  
  status = "s.idle"
  num.list = character(0)

  iline = 0
  for (line in lines) {
    iline = iline+1

    line = sub("\\s+$", "", line)
    empty.line = grepl("^\\s*$", line)

    if (status == "s.idle") {
      rc = str.match.just.single(line, "^[*]\\s+(\\S+)\\s+:\\s+(.*\\S)\\s+:\\s+(.*\\S)\\s*$")
      if (!is.na(rc[1])) {
        ##
        ## beginning of quantity definition, stode gamma, node, descr
        ##
        node = rc[1]
        descr = rc[2]
        gamma = rc[3]
        
        quantities[[gamma]]$node = node
        quantities[[gamma]]$descr = descr
        
        status = "s.begin.def"
        next
      }
    }
    
    if (status == "s.begin.def") {
      rc = str.match.just.single(line, "^(\\S+)\\s+=\\s+(.*\\S)\\s*$")
      if (!is.na(rc[1])) {
        ##
        ## get equation that defines gamma
        ##
        status == "s.begin.formula"
        gamma.eq = line
        status = "s.idle"

        ##--- add spaces around to facilitate matches
        rc = paste(" ", line, " ", sep="")
        ##--- replace whitespace with single space
        rc = gsub("\\s+", " ", rc, perl=TRUE)
        ##--- remove spaces after braket
        rc = gsub("\\(\\s+", "(", rc, perl=TRUE)
        ##--- remove trailing zeros after a dot and some digits
        rc = gsub("([.]\\d+[1-9])0+(\\D)", "\\1\\2", rc, perl=TRUE)
        ##--- remove trailing zeros immediately after a dot
        rc = gsub("([.])0+(\\D)", "\\1\\2", rc, perl=TRUE)
        ##--- remove trailing "."
        rc = gsub("[.](\\D)", "\\1", rc, perl=TRUE)
        ##--- remove "*1"
        rc = gsub("[*]1([^\\d.])", "\\1", rc, perl=TRUE)
        ##--- remove "1*"
        rc = gsub("(\\D)1[*]([^\\d.])", "\\1\\2", rc, perl=TRUE)
        ##--- remove unneeded brakets
        rc = gsub("\\(([^()+-]+)\\)", "\\1", rc, perl=TRUE)
        ##--- remove added spaces
        rc = gsub("^\\s", "", rc)
        rc = gsub("\\s+$", "", rc)

        ##
        ## replace numbers with expressions using symbolic parameters
        ##
        if (TRUE) {
          ##--- substitute numeric values

          if (gamma == "Gamma1") {
            rc = gsub("Gamma44*0.5000000000", "Gamma44*BRA_Kzbar_KL_KET", rc, fixed=TRUE)
          } else if (gamma == "Gamma7") {
            rc = gsub("Gamma35*0.5000000000", "Gamma35*BRA_Kzbar_KL_KET", rc, fixed=TRUE)
            rc = gsub("Gamma37*0.5000000000", "Gamma37*BRA_Kz_KL_KET", rc, fixed=TRUE)
          } else if (gamma == "Gamma33") {
            rc = gsub("Gamma35*0.5000000000", "Gamma35*BRA_Kzbar_KS_KET", rc, fixed=TRUE)
            rc = gsub("Gamma37*0.5000000000", "Gamma37*BRA_Kz_KS_KET", rc, fixed=TRUE)
            rc = gsub("Gamma40*0.5000000000", "Gamma40*BRA_Kzbar_KS_KET", rc, fixed=TRUE)
            rc = gsub("Gamma42*0.5000000000", "Gamma42*BRA_Kz_KS_KET", rc, fixed=TRUE)
            rc = gsub("Gamma44*0.5000000000", "Gamma44*BRA_Kzbar_KS_KET", rc, fixed=TRUE)
          }
          
          rc = gsub("0.0153", "BR_om_pimpip", rc, fixed=TRUE)
          rc = gsub("0.0828", "BR_om_pizgamma", rc, fixed=TRUE)
          rc = gsub("0.09418761", "(BR_KS_2piz*BR_KS_2piz)", rc, fixed=TRUE)
          rc = gsub("0.1263054152", "(BR_phi_KSKL*BR_KS_2piz)/(BR_phi_KmKp+BR_phi_KSKL)", rc, fixed=TRUE)
          rc = gsub("0.15345", "(BRA_Kz_KS_KET*BR_KS_2piz)", rc, fixed=TRUE)
          rc = gsub("0.2274", "BR_eta_pimpippiz", rc, fixed=TRUE)
          rc = gsub("0.26128673", "(BRA_Kzbar_KL_KET*BR_eta_pimpippiz + BRA_Kzbar_KS_KET*BR_KS_2piz*BR_eta_pimpippiz + BRA_Kzbar_KS_KET*BR_KS_pimpip*BR_eta_3piz)", rc, fixed=TRUE)
          rc = gsub("0.281", "BR_eta_charged", rc, fixed=TRUE)
          rc = gsub("0.2847942238", "(BR_phi_KSKL*BR_KS_pimpip)/(BR_phi_KmKp+BR_phi_KSKL)", rc, fixed=TRUE)
          rc = gsub("0.3069", "BR_KS_2piz", rc, fixed=TRUE)
          rc = gsub("0.3257", "BR_eta_3piz", rc, fixed=TRUE)
          rc = gsub("0.346", "(BRA_Kz_KS_KET*BR_KS_pimpip)", rc, fixed=TRUE)
          rc = gsub("0.3595", "(BRA_Kzbar_KS_KET*BR_eta_neutral)", rc, fixed=TRUE)
          rc = gsub("0.4115523466", "BR_phi_KSKL/(BR_phi_KmKp+BR_phi_KSKL)", rc, fixed=TRUE)
          rc = gsub("0.4247496", "(2*BR_KS_pimpip*BR_KS_2piz)", rc, fixed=TRUE)
          rc = gsub("0.5531", "(BR_eta_3piz+BR_eta_pimpippiz)", rc, fixed=TRUE)
          rc = gsub("0.5884476534", "BR_phi_KmKp/(BR_phi_KmKp+BR_phi_KSKL)", rc, fixed=TRUE)
          rc = gsub("0.65345", "(BRA_Kzbar_KS_KET*BR_KS_2piz+BRA_Kzbar_KL_KET)", rc, fixed=TRUE)
          rc = gsub("0.692", "BR_KS_pimpip", rc, fixed=TRUE)
          rc = gsub("0.719", "BR_eta_neutral", rc, fixed=TRUE)
          rc = gsub("0.8732418773", "(BR_phi_KmKp + BR_phi_KSKL*BR_KS_pimpip)/(BR_phi_KmKp+BR_phi_KSKL)", rc, fixed=TRUE)
          rc = gsub("0.892", "BR_om_pimpippiz", rc, fixed=TRUE)
          rc = gsub("0.9073", "(BR_om_pimpippiz+BR_om_pimpip)", rc, fixed=TRUE)
          rc = gsub("1.09418761", "((BR_KS_2piz*BR_KS_2piz) + BRA_KzKzbar_KLKL_KET_by_BRA_KzKzbar_KSKS_KET)", rc, fixed=TRUE)
        }
        
        gamma.eq = rc
        parts = str_match(rc, "(.*\\S)\\s*=\\s*(\\S.*)")
        
        ##--- by default no definition
        quantities[[gamma]]$def = ""
        if (parts[2] != parts[3]) {
          ##--- if the equation is not an identity, set definition
          quantities[[gamma]]$def = parts[3]
          ##+++ remove wrong def for Gamma115
          if (gamma == "Gamma115") quantities[[gamma]]$def = ""
        }
      }
      next
    }
  }

  return(invisible(quantities))
}

##
## process quantities read from readpdg.cc output file
##
process.quantities = function(quantities) {
  quant.names = names(quantities)
  quant.order = order(alucomb2.gamma.num.id(quant.names))
  quant.names = quant.names[quant.order]
  num.list = numeric(0)
  constr.lin.val = numeric(0)
  constr.lin.expr = list()
  constr.nl.val = numeric(0)
  constr.nl.expr = list()

  ##
  ## list containing parameters initializazion
  ## can be used to evaluate symbolic expressions
  ## is not used in the present code
  ##
  params = list()
  params = within(params, {
    BR_eta_2gam        = 39.31e-2 # ( 39.31 +- 0.20 ) %
    BR_eta_neutral     = 71.90e-2 # ( 71.90 +- 0.34 ) %
    BR_eta_3piz        = 32.57e-2 # ( 32.57 +- 0.23 ) %
    BR_eta_pimpippiz   = 22.74e-2 # ( 22.74 +- 0.28 ) %
    BR_eta_charged     = 28.10e-2 # ( 28.10 +- 0.34 ) %
    BR_KS_2piz         = 30.69e-2 # ( 30.69 +- 0.05 ) %
    BR_KS_pimpip       = 69.20e-2 # ( 69.20 +- 0.05 ) %
    BR_om_pimpippiz    = 89.2e-2  # ( 89.2  +- 0.7 ) %
    BR_om_pimpip       = 1.53e-2  # ( 1.53  +  0.11 - 0.13 ) %
    BR_om_pizgamma     = 8.28e-2  # ( 8.28  +- 0.28 ) %
    BR_phi_KmKp        = 48.9e-2  # ( 48.9  +- 0.5 ) %
    BR_phi_KSKL        = 34.2e-2  # ( 34.2  +- 0.4 ) %
    
    BRA_Kz_KS_KET = 1/2 # |<Kz | KS>|^2
    BRA_Kz_KL_KET = 1/2 # |<Kz | KS>|^2
    BRA_Kzbar_KS_KET = 1/2 # |<Kzbar | KS>|^2
    BRA_Kzbar_KL_KET = 1/2 # |<Kzbar | KS>|^2
    
    BRA_KzKzbar_KLKL_KET_by_BRA_KzKzbar_KSKS_KET = 1
  })
  
  for(quant.name in quant.names) {
    quant = quantities[[quant.name]]
    constr.card = ""
    quant.def.vars = character(0)
    
    quant.def = quant$def
    if (quant.def != "") {
      quant.def.expr = parse(text=quant.def)
      ## quant.def.expr.subs = esub.expr(quant.def.expr, params)    
      ## quant.def = paste(str_trim(deparse(quant.def.expr.subs)), collapse="")
    }
    
    if (quant.def != "") {
      nlconstr.str =  paste("-", quant.name, " + ", quant.def, sep="")
      nlconstr.expr = parse(text=nlconstr.str)
      
      ##--- equation is linear if there is no variable dependence after derivation
      linear.flag = TRUE
      coeffs = list()
      vars = all.vars(nlconstr.expr)
      for (var in vars) {
        nlconstr.deriv = D(nlconstr.expr, var)
        if (length(all.vars(nlconstr.deriv)) > 0) {
          linear.flag = FALSE
        } else {
          coeffs[[var]] = eval(nlconstr.deriv)
        }
      }
      ##--- save all quantities in quantity def
      quant.def.vars = setdiff(unique(vars), quant.name)
      quant.def.vars = quant.def.vars[quant.def.vars %in% quant.names]
      quant.def.vars = quant.def.vars[order(alucomb2.gamma.num.id(quant.def.vars))]
      
      if (FALSE && linear.flag) {
        ##--- if linear constraint, write COMBOFQUANT card
        constr.card = paste("COMBOFQUANT", vars[1], paste(coeffs[vars[-1]], vars[-1], collapse=" "))
        
        constr.lin.val = c(constr.lin.val, 0)
        lin.comb = coeffs[vars[-1]]
        names(lin.comb) = vars[-1]
        constr.lin.comb[[paste(vars[1], "coq", sep=".")]] = lin.comb
      } else {
        ##--- write NLCONSTRAINT non-linear constraint
        constr.card = paste(
          "NLCONSTRAINT",
          paste(quant.name, "c" ,sep="."),
          "0",
          paste("\"", nlconstr.str, "\"", sep="")
          )
      
        constr.nl.val = c(constr.nl.val, 0)
        constr.nl.expr[[paste(vars[1], "c", sep=".")]] = nlconstr.str
      }
      
      ##--- get all the numbers in the expression
      matches = gregexpr("[*]([0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?)", quant.def, perl=TRUE)
      all.numbers = sub("^[*]", "", unlist(regmatches(quant.def, matches)))
      if (length(all.numbers) > 0) {
        all.numbers = unique(all.numbers)
        all.numbers = all.numbers[!(all.numbers %in% num.list)]
        num.list = c(num.list, all.numbers)
      }
    }

    ##
    ## print collected information as alubomb2.r cards
    ##
    cat("#\n")
    cat("# ", quant.name, " = ", quant$descr, "\n", sep="")
    cat("#\n")
    cat("QUANTITY ", quant.name, " node ", quant$node, " descr \"", quant$descr, "\"\n", sep="")
    if (constr.card != "") {
      cat("#\n")
      cat("# ", quant.name, " = ", quant$def, "\n", sep="")
      cat("#\n")
      if (TRUE) {
        for (var in quant.def.vars) {
          descr = quantities[[var]]$descr
          if (!is.null(descr)) {
            cat("# ", var, " = ", descr, "\n", sep="")
          }
        }
        cat("#\n")
      }
      cat(constr.card, "\n", sep="")
    }
    cat("\n")
  }
  
  num.list.order = order(as.numeric(num.list))
  num.list = num.list[num.list.order]

  ##--- print list of numbers used in the constraint equations
  if (FALSE && length(num.list)>0) {
    cat("numeric constants in equations\n\n")
    cat(num.list, sep="\n")
  }

  ##--- print collected information
  if (FALSE) {
    cat("\n")
    print(quantities)
    
    print(constr.lin.val)
    print(constr.lin.comb)
    
    print(constr.nl.val)
    print(constr.nl.expr)
  }
  
  return(num.list)
}

##
## read file produced by readpdg.cc containing all quantities definitions
##
quantities = read.info("../../2009/TauFit_Mar2011/all_node_def.txt")

##
## process quantities and prinf alucomb2.r cards
##
num.list = process.quantities(quantities)
