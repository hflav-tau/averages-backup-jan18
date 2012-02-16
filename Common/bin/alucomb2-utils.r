#!/usr/bin/env Rscript

## ////////////////////////////////////////////////////////////////////////////
## definitions

##
## create diagonal matrix also for vectors of length one
##
diag.m <- function(vec) {
  if (length(vec) <= 1) {
    rc = as.matrix(vec)
  } else {
    rc = diag(vec)
  }
  rc
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

##--- measurement names matching a specific name
meas.match = function(gamma) {
  return( regexpr(paste("[.]", gamma, "[.]", sep=""), meas.names) != -1 )
}

##--- sum in quadrature
quadrature = function(x) {
  return(sqrt(sum(x^2)))
}

##--- normalize publication state
alu.norm.pubstate = function(str) {
  str = sub("^pub[^.]*", "pub", str, perl=TRUE)
  str = sub("^sub[^.]*", "sub", str, perl=TRUE)
  str = sub("^acc[^.]*", "acc", str, perl=TRUE)
  str = sub("^prelim[^.]*", "prelim", str, perl=TRUE)
  str = sub("^PUB[^.]*", "PUB", str, perl=TRUE)
  str = sub("^SUB[^.]*", "SUB", str, perl=TRUE)
  str = sub("^ACC[^.]*", "ACC", str, perl=TRUE)
  str = sub("^PRELIM[^.]*", "PRELIM", str, perl=TRUE)
  return(str)
}

##
## cards
##

##
## measurement block:
##
## old style
##   BEGIN <exp> <quantity> <pub|prelim> <reference>.[<other.tags>]
##   MEASUREMENT    m_<quantity>  statistical  systematic 
##   DATA           m_<quantity>  statistical  systematic 
##                  <val>         <stat>       <syst>
##   ...
##   STAT_CORR_WITH 0.235119 BaBar Gamma9by5 pub AUBERT.10P
##   ...
##   ##--- systematics terms
##   DATA      
##     term1 <signed val>
##     term2 <signed val>
##   END
## new style
##   BEGIN MEASUREMENT <exp> <quantity> <pub|prelim> <reference> [<other tags>]
##   VALUE <val> <stat> <syst>
##   ...
##   SYSTEMATICS
##     term1 <signed val>
##     term2 <signed val>
##   ...
##   END
## 
## a "measurement" is stored in a list with fields
##
## - tags = character array: measurement name fields
## - quant = character: measured quantity
## - value = numeric,
## - stat = numeric,
## - syst = numeric,
## - params = list of numeric arrays: parameter value, pos/neg uncertainty
## - syst.terms = numeric array: syst. contributions related to external parameters
## - corr.terms = numeric array: stat. corr. with other measurements
## - corr.terms.tot = numeric array: error corr. with other measurements
## 

##
## combination block:
##
## old style
##   BEGIN <tag> [<tag 2> ...]
##   COMBINE * * *
##   PARAMETERS !      assumed_value  pos_excurs   neg_excurs
##     par1            0.919          +0.003       -0.003
##   ...
##   MEASUREMENT m_<quantity> statistical systematic [seed <value>] [scale <value>]
##   ...
##   END
# new style
##   BEGIN COMBINATION <tag> [<tag 2> ...]
##   PARAMETERS !      assumed_value  pos_excurs   neg_excurs
##     par1            0.919          +0.003       -0.003
##   ...
##   QUANTITY <quantity> [seed <value>] [scale <value>]
##   ...
##   END
##
## a "combination" is stored as a list with fields:
##
## - tags = character array: combination name fields
## - quantities = character array: measured quantities to be fitted
## - params = list of numeric arrays: parameter value, pos/neg uncertainty
## - constr.lin.comb = list of numeric arrays: coeffs of quantities in linear constraints
## - constr.lin.val = list or numeric: values of linear constraint equations
## - constr.nl.expr = list of numeric arrays: coeffs of quantities in linear constraints
## - constr.nl.val = list or character: expression of non-linear constraint
##

## ////////////////////////////////////////////////////////////////////////////
## code

##
## alucomb.read
##

alucomb.read = function(file = "") {

if (match.nocase("^\\s*$", file)) {
  stop("alucomb: please provide as argument the card input file\n")
}

flag.replace.corr = FALSE
flag.build.delta = FALSE

if (!file.exists(file)) {
  stop("cannot find file ", file, "\n")
}
dir.base = dirname(file)
lines = readLines(file)
cat("read file", file, "\n")

##
## deal with INCLUDE lines
##
iline = 1
lines.len = length(lines)
repeat {
  if (iline > lines.len) break
  fields = unlist(strsplit(lines[iline], "\\s+", perl=TRUE))
  if (length(fields) > 0 &&
      match.nocase("^INCLUDE$", fields[1])) {
    file.inc = paste(dir.base, fields[2], sep="/")
    if (!file.exists(file.inc)) {
      stop("cannot find included file ", file.inc, "\n")
    }
    lines.inc = readLines(file.inc)
    cat("read file", file.inc, "\n")
    lines = c(
      if (iline-1>=1) {lines[1:(iline-1)]} else {NULL},
      lines.inc,
      if (iline+1<=lines.len) {lines[(iline+1):(lines.len)]} else {NULL}
      )
    lines.len = length(lines)
    next
  }
  iline = iline + 1
}

##
## init global storage for cards
##
measurements = list()
combination = list()

s.inblock = FALSE
s.inclause = FALSE

t.inblock = FALSE
t.inclause = FALSE
t.continue = FALSE
t.endclause = FALSE
t.endblock = FALSE
t.endfile = FALSE

iline = 1
lines.len = length(lines)
repeat {
  ##
  ## read new line
  ##
  if (iline > lines.len) {
    t.endfile = TRUE
  } else {
    line = lines[iline]
    iline = iline+1

    ##--- skip comment and empty lines
    if (regexpr("^\\s*$", line, perl=TRUE) != -1 ||
        regexpr("^[\\#!*;]", line, perl=TRUE) != -1) {
      next
    }
    ##--- remove end-of-line comments
    line = sub("\\s*[\\#!].*", "", line, perl=TRUE)
    fields = unlist(strsplit(line, "\\s+", perl=TRUE))
    fields[1] = toupper(fields[1])
    
    ##--- block begin
    t.inblock = (fields[1] == "BEGIN")
    ##--- block end
    t.endblock = (fields[1] == "END")
    ##--- a line beginning with space is the continuation of the last clause
    t.continue = (fields[1] == "")
    ##--- begin a clause if: 1) does not begin/end a block and 2) does not begin with blank
    t.inclause = !(t.inblock || t.endblock || t.continue)
  }

  ##--- when continuing a clause collect all its fields
  if (t.continue) {
    if (s.inclause) {
      ##--- continuation of previous clause
      clause.fields = c(clause.fields, fields[-1])
    } else {
      t.continue = FALSE
      cat("warning, continuation line outside of clause, ignored, line ...\n")
      cat("  '", line, "'\n", sep="")
    }
    next
  }
  
  ##
  ## end of a clause -> interpret its fields
  ##
  if (s.inclause && !t.continue) {
    ##--- if no continuation then the current clause ends
    s.inclause = FALSE
    if (t.inblock) {
      cat("warning, clause ended by BEGIN block, line ...\n")
      cat("  '", line, "'\n", sep="")
    }
    if (t.endfile) {
      cat("warning, clause ended by end of file, clause...\n")
      cat("  ", paste(clause.fields, collapse=" "), "\n", sep="")
    }
    ##
    ## handle a complete clause here
    ##
    ## cat("END clause", unlist(clause.fields), "\n\n")
    clause.keyw = toupper(clause.fields[1])
    clause.fields = clause.fields[-1]

    ##--- values are numbers possible preceded by "+", "-", "+-"
    clause.fields.nosign = sub("([+]|-|[+]-)", "", clause.fields)
    clause.fields.signtag = sub("^([+]|-|[+]-|).*$", "\\1", clause.fields)
    clause.fields.val = suppressWarnings(as.numeric(clause.fields))
    ##--- labels are all non-numeric words
    clause.labels = clause.fields[is.na(clause.fields.val)]
    
    ##--- indices to numeric values preceded by plus sign
    plus.ind = grep("+", clause.fields.signtag, fixed=TRUE)
    plus.ind = plus.ind[!is.na(clause.fields.val[plus.ind])]
    ##--- indices to numeric values preceded by minus sign
    minus.ind = grep("-", clause.fields.signtag, fixed=TRUE)
    minus.ind = minus.ind[!is.na(clause.fields.val[minus.ind])]
    ##--- apply minus signs
    clause.fields.val[minus.ind] = -clause.fields.val[minus.ind]

    ##--- indices to sequences +<numeric value 1> -<numeric value 2>
    plus.minus.ind = intersect(plus.ind, minus.ind-1)
    ##--- add negative value as attribute to positive value
    plus.minus.val = mapply(
      function(first, second) {
        attr(first, "negval", second)
        return(first)
      }, clause.fields.val[plus.minus.ind], clause.fields.val[plus.minus.ind+1])
    ##--- replace +<numeric value 1> -<numeric value 2> with a single number
    clause.fields.val = as.list(clause.fields.val)
    clause.fields.val[plus.minus.ind] = plus.minus.val
    clause.fields.val[plus.minus.ind+1] = NA
    pm.ind = grep("+-", clause.fields.signtag, fixed=TRUE)
    clause.fields.val[pm.ind] = lapply(clause.fields.val[pm.ind], function(val) {attr(val, "pm") = TRUE; val})
    clause.values = clause.fields.val[!is.na(clause.fields.val)]

    if (clause.keyw == "COMBINE") {
      ##
      ## COMBINE
      ##

      ##+++ combos
      block.type = "COMBINATION"

    } else if (clause.keyw == "PARAMETERS") {
      ##
      ## PARAMETERS (same treatment in MEASUREMENT and COMBINATION block)
      ##

      ##--- after collapsing +devp -devm into devp with attr negval expect 2x values as labels
      if (length(clause.values) != 2*length(clause.labels)) {
        stop("wrong parameter data in line...\n  ", paste(c(clause.keyw, clause.fields), collapse=" "))
      }
      param.val = clause.values[seq(1, 2*length(clause.labels), by=2)]
      names(param.val) = clause.labels
      param.delta.p = clause.values[seq(2, 2*length(clause.labels), by=2)]
      param.delta.n = lapply(param.delta.p, function(delta) {
        ifelse(is.null(attr(delta, "negval")), -delta, attr(delta, "negval"))
      })
      block$params = mapply(
        function(val, delta.p, delta.n) {
          c("value"=val, "delta_pos"=delta.p, "delta_neg"=delta.n)
        }, param.val, param.delta.p, param.delta.n, SIMPLIFY=FALSE)
      print(block$params)
      print(param.val)
      
    } else if (clause.keyw == "QUANTITY" || (block.type == "COMBINATION" && clause.keyw == "MEASUREMENT")) {
      ##
      ## QUANTITY or (+++combos) MEASUREMENT in COMBINE block
      ##
      
      ##+++ combos
      meas.name = sub("^m_", "", clause.fields[1], ignore.case=TRUE)

      if (!is.null(block$quantities[[meas.name]])) {
        stop("duplicate quantity ", meas.name, " in same block")
      }
      block$quantities[[meas.name]] = list()
      clause.labels = clause.labels[-1]

      ##+++ combos, remove labels named "stat*", "syst*"
      clause.labels = clause.labels[!match.nocase("^(syst|stat)", clause.labels)]

      if (length(clause.labels) != length(clause.values)) {
        stop("mismatch between labels and numeric values in line...\n  ", paste(c(clause.keyw, clause.fields), collapse=" "))
      }
      if (length(clause.values) > 0) {
        names(clause.values) = clause.labels
        block$quantities[[meas.name]] = clause.values
      }

    } else if (clause.keyw == "MEASUREMENT") {
      ##
      ##+++ MEASUREMENT clause, combos compatibility
      ##
      block.type = "MEASUREMENT"
      meas.name = sub("^m_", "", clause.fields[1], ignore.case=TRUE)
      if (meas.name != block.fields[2]) {
        stop("MEASUREMENT quantity does not match BLOCK quantity...\n  ", paste(c(clause.keyw, clause.fields), collapse=" "))
      }
      measurement.in.data = TRUE
    } else if (clause.keyw == "VALUE" || (clause.keyw == "DATA" && measurement.in.data)) {
      ##
      ## VALUE --> enter measurement value
      ##+++ DATA after MEASUREMENT clause, combos compatibility
      ##
      
      if (!is.null(block.meas$val)) {
        stop("duplicated measurement in line...", paste(c(clause.keyw, clause.fields), collapse=" "))
      }

      if (length(clause.values) == 3) {
        if (!is.null(attr(clause.values[[1]], "negval"))) {
          stop("+-values for measurements value in line...\n  ", paste(c(clause.keyw, clause.fields), collapse=" "))
        }
        ##--- use values if existing
        block.meas$value = clause.values[[1]]

        block.meas$stat.p = clause.values[[2]]
        negval = attr(clause.values[[2]], "negval")
        if (!is.null(negval)) {
          block.meas$stat.n = negval
          block.meas$stat = sqrt((block.meas$stat.p^2 + block.meas$stat.n^2)/2)
        } else {
          block.meas$stat.n = block.meas$stat.p
          block.meas$stat = block.meas$stat.p
        }

        block.meas$syst.p = clause.values[[3]]
        negval = attr(clause.values[[3]], "negval")
        if (!is.null(negval)) {
          block.meas$syst.n = negval
          block.meas$syst = sqrt((block.meas$syst.p^2 + block.meas$syst.n^2)/2)
        } else {
          block.meas$syst.n = block.meas$syst.p
          block.meas$syst = block.meas$syst.p
        }
        ##+++ set end of MEASEREMENT DATA (combos compatibility)
        measurement.in.data = FALSE
      } else {
        stop("wrong number of numeric data in line...\n  ", paste(c(clause.keyw, clause.fields), collapse=" "))
      }
    } else if (clause.keyw == "DATA" ||
               clause.keyw == "SYSTEMATICS" || clause.keyw == "SYSTLOCAL" || clause.keyw == "SYSTPAPER" ) {
      ##
      ## SYSTEMATICS or (+++combos) DATA
      ## list systematic terms with a label that possibly refers
      ## to a global parameter or to the same effect in another
      ## measurement
      ##
      ## SYSTLOCAL like SYSTEMATICS but add measurement name to make it local
      ## SYSTPAPER like SYSTEMATICS but add reference name to make it paper-wide
      ##

      ##
      ## - "%" on either name or value means percent of measurement value
      ## - "&" on either name or value means relative of measurement value
      ##
      clause.df = data.frame(
        field = clause.fields,
        rel = grepl("^(.*)[&]$", clause.fields, perl=TRUE),
        perc = grepl("^(.*)[%]$", clause.fields, perl=TRUE),
        num = suppressWarnings(as.numeric(clause.fields)),
        stringsAsFactors=FALSE)
      clause.labels = sub("^(.*)[&%]$", "\\1", with(clause.df, field[is.na(num)]))
      clause.values = with(clause.df, num[!is.na(num)])
      if (length(clause.values) != length(clause.labels)) {
        stop("mismatch between labels and numeric values in line...\n  ", paste(c(clause.keyw, clause.fields), collapse=" "))
      }
      clause.perc = with(clause.df, perc[is.na(num)]) | with(clause.df, perc[!is.na(num)])
      clause.rel = with(clause.df, rel[is.na(num)]) | with(clause.df, rel[!is.na(num)])
      clause.values =
        ifelse(clause.perc, clause.values*block.meas$value/100,
               ifelse(clause.rel, clause.values*block.meas$value,
                      clause.values))

      if (!is.null(block.meas$syst.terms)) {
        stop("duplicated systematic terms in line...\n  ", paste(c(clause.keyw, clause.fields), collapse=" "))
      }
      if (clause.keyw == "SYSTLOCAL") {
        clause.labels = paste(paste(block.fields, collapse="."), clause.labels, sep=".")
      } else if (clause.keyw == "SYSTPAPER") {
        clause.labels = paste(paste(block.fields[-2], collapse="."), clause.labels, sep=".")
      }
      names(clause.values) = clause.labels
      block.meas$syst.terms = unlist(clause.values)
      
    } else if (clause.keyw == "STAT_CORR_WITH" || clause.keyw == "ERROR_CORR_WITH") {
      ##
      ## STAT_CORR_WITH, ERROR_CORR_WITH
      ##

      if (length(clause.values) != 1) {
        stop("exactly one numeric value required in line...\n  ", paste(c(clause.keyw, clause.fields), collapse=" "))
      }
      corr = clause.values[[1]]

      clause.labels[3] = alu.norm.pubstate(clause.labels[3])
      names(corr) =  paste(clause.labels, collapse=".")
      if (clause.keyw == "STAT_CORR_WITH") {
        block.meas$corr.terms = c(block.meas$corr.terms, corr)
      } else {
        block.meas$corr.terms.tot = c(block.meas$corr.terms.tot, corr)
      }
    } else if (clause.keyw == "SUMOFQUANT" || clause.keyw == "COMBOFQUANT" || clause.keyw == "CONSTRAINT") {
      ##
      ## SUMOFQUANT, COMBOFQUANT, CONSTRAINT
      ##
      
      if (clause.keyw != "CONSTRAINT") {
        constr.name = paste(clause.fields[1], "coq", sep=".")
        constr.val = 0
      }

      if (clause.keyw == "SUMOFQUANT") {
        constr.comb = c(-1, rep(1, length(clause.fields)-1))
        names(constr.comb) = clause.fields
      } else if (clause.keyw == "COMBOFQUANT") {
        if (length(clause.labels)-1 != length(clause.values)) {
          stop("error: mismatch between labels and data values in line...\n  ", paste(c(clause.keyw, clause.fields), collapse=" "))
        }
        constr.comb = c(-1, unlist(clause.values))
        names(constr.comb) = clause.labels
      } else if (clause.keyw == "CONSTRAINT") {
        if (length(clause.labels) != length(clause.values)) {
          stop("error: mismatch between labels and data values in line...\n  ", paste(c(clause.keyw, clause.fields), collapse=" "))
        }
        constr.name = clause.labels[1]
        constr.val = clause.values[[1]]
        constr.comb = unlist(clause.values[-1])
        names(constr.comb) = clause.labels[-1]
      }

      block$constr.lin.val[[constr.name]] = constr.val
      block$constr.lin.comb[[constr.name]] = constr.comb

    } else if (clause.keyw == "NLCONSTRAINT") {
      ##
      ## NLCONSTRAINT
      ##
      
      if (is.na(suppressWarnings(as.numeric(clause.fields[2])))) {
        stop("error: missing numeric value in non-linear constraint in...\n  ", paste(c(clause.keyw, clause.fields), collapse=" "))
      }
      constr.name = clause.labels[1]
      constr.val = clause.values[[1]]
      constr.expr = paste(clause.fields[-(1:2)], collapse=" ")
      constr.expr = sub("^\\s*\"?([^\"]*)\"?\\s*$", "\\1", constr.expr, perl=TRUE)
      block$constr.nl.val[[constr.name]] = constr.val
      block$constr.nl.expr[[constr.name]] = constr.expr
      
      ## } else if (clause.keyw == "") {
    }
  }
  
  if (t.inclause) {
    ##--- new clause (previous one has been handled already)
    t.inclause = FALSE
    if (s.inblock) {
      s.inclause = TRUE
      clause.fields = fields
    } else {
      cat("warning, begin clause outside of block, ignored, line ...\n")
      cat("  '", line, "'\n", sep="")
    }
  }
  
  if (t.inblock) {
    t.inblock = FALSE
    if (s.inblock) {
      cat("error, nested BEGIN block, ignored, line ...\n")
      cat("  '", line, "'\n", sep="")
      next
    }
    s.inblock = TRUE

    ##--- remove "BEGIN"
    block.fields = fields[-1]
    
    block.type = ""
    if (length(block.fields) > 0) {
      block.type = toupper(block.fields[1])
      if (block.type == "COMBINATION" || block.type == "MEASUREMENT") {
        block.fields = block.fields[-1]
      } else {
        ##+++ temporarily just warning for COMBOS-like cards compatibility
        if (TRUE) {
          cat("warning, unknown BEGIN block in line...\n")
          cat(" ", fields, "\n")
        }
      }
    }
    block.meas = list()
    ##+++ block.meas$syst.terms = numeric(0)
    block.meas$corr.terms = numeric(0)
    block.meas$corr.terms.tot = numeric(0)

    block = list()
    block$params = list()
    block$quantities = list()
    block$constr.lin.val = list()
    block$constr.lin.comb = list()
    block$constr.nl.val = list()
    block$constr.nl.expr = list()
  }
  
  if (s.inblock && (t.endblock || t.endfile)) {
    ##--- end block
    t.endblock = FALSE
    s.inblock = FALSE
    ## cat("END block", unlist(block.fields), "\n\n")
    if (block.type == "MEASUREMENT") {
      ##
      ## MEASUREMENT BLOCK
      ##
      
      ##+++combos
      if (FALSE && length(block.fields) == 3) {
        block.fields = c(block.fields[-3], strsplit(block.fields[3], "[.]", perl=TRUE)[[1]])
      }

      if (length(block.fields) < 3) {
        stop("too few fields in line...\n  ", paste(c("BEGIN", block.type, block.fields), collapse=" "))
      }
      block.fields[3] = alu.norm.pubstate(block.fields[3])

      meas.tag = paste(block.fields, collapse=".")
      block.meas$tags = block.fields

      if (!is.null(measurements[[meas.tag]])) {
        stop("duplicate measurement in line...\n  ", paste(c("BEGIN", block.type, block.fields), collapse=" "))
      }
      block.meas$params = block$params
      ##
      ## measured quantity = 2nd tag in the BEGIN MEASUREMENT statement
      ## BEGIN [MEASUREMENT] <esperiment> <quantity> <pub|prelim> <reference> [...]
      ##
      block.meas$quant = block.fields[2]
      measurements[[meas.tag]] = block.meas
      
    } else if (block.type == "COMBINATION") {
      ##
      ## COMBINE BLOCK
      ##
      block$tags = unlist(block.fields)
      combination = block
      
    } else {
      stop("unknown block type in line...\n  ", paste(c("BEGIN", block.type, block.fields[-1]), collapse=" "))
    }
  }
  
  if (t.endblock && !s.inblock) {
    t.endblock = FALSE
    cat("error, END block without BEGIN, ignored, line ...\n")
    cat("  '", line, "'\n", sep="")
  }

  if (t.endfile) {
    t.endfile = FALSE
    break
  }
}

list(combination=combination, measurements=measurements)
} ##--- end of alucomb.read

##
## transform the name of a measurement in the HFAG-tau format into a format usable by Root
##
hfag.to.root = function(str) {
  if (str == "TauMass") {
    return("m_{#tau} [MeV]")
  }
  
  str.orig = str

  str = gsub("Numbar", "#bar{#nu}_{#mu}",str)  

  str = gsub("(^|[^A-Z])([A-Z][a-y]*)([z+-]*)b", "\\1#bar{\\2}\\3", str, perl=TRUE)
  str = gsub("F1", "f_{1}",str)
  str = gsub("Pi", "#pi",str)

  str = gsub("Nue", "#nu_{e}",str)
  str = gsub("Num", "#nu_{#mu}",str)
  str = gsub("Nut", "#nu_{#tau}",str)
  str = gsub("Nu", "#nu",str)

  str = gsub("M", "#mu",str)
  str = gsub("H", "h",str)
  str = gsub("m($|[#}A-Zh])", "^{-}\\1", str, perl=TRUE)
  str = gsub("p($|[#}A-Zh])", "^{+}\\1", str, perl=TRUE)
  str = gsub("z($|[#}A-Zh])", "^{0}\\1", str, perl=TRUE)

  if (str.orig %in% c("HmHmHpNu", "PimKmPipNu", "PimPimPipNu")) {
    str = paste(str, "(ex.K^{0})")
  }

  str = paste("B(#tau^{-} #rightarrow ", str, ")", sep="")
  return(str)
}

##
## compute root of symmetric real positive semi-definite matrix
##
alu.matr.sqrt.symm.semipos = function(X, tol = sqrt(.Machine$double.eps)) {
  if (!is.matrix(X)) 
    X <- as.matrix(X)
  if (length(dim(X)) > 2L || !(is.numeric(X))) 
    stop("'X' must be a numeric or complex matrix")

  X = (X + t(X))/2
  rc.e = eigen(X)
  comp.positive = rc.e$values > tol*rc.e$values[1]
  rc.vals = ifelse(comp.positive, rc.e$values, 0)
  rc = rc.e$vectors %*% diag.m(sqrt(rc.vals)) %*% t(rc.e$vectors)
  rc = (rc + t(rc))/2
  rownames(rc) = rownames(X)
  colnames(rc) = colnames(X)

  return(rc)
}

##
## compute pseudo-inverse square root of symmetric real positive semi-definite matrix
## zero eigenvalues produce zero eigenvalues in the inverse-square matrix
##
alu.matr.inv.sqrt.symm.semipos = function(X, tol = sqrt(.Machine$double.eps)) {
  if (!is.matrix(X)) 
    X <- as.matrix(X)
  if (length(dim(X)) > 2L || !(is.numeric(X))) 
    stop("'X' must be a numeric or complex matrix")

  X = (X + t(X))/2
  rc.e = eigen(X)
  ##--- force to zero computationally zero or negative eigenvalues
  comp.positive = rc.e$values > tol*rc.e$values[1]
  rc.vals = ifelse(comp.positive, rc.e$values, 0)
  rc.vecs.inv = ifelse(rc.vals > 0, 1/rc.vals, 0)
  
  rc.inv = rc.e$vectors %*% diag.m(rc.vecs.inv) %*% t(rc.e$vectors)
  rc.inv = (rc.inv + t(rc.inv))/2
  rownames(rc.inv) = rownames(X)
  colnames(rc.inv) = colnames(X)
  rc = rc.e$vectors %*% diag.m(sqrt(rc.vecs.inv)) %*% t(rc.e$vectors)
  rc = (rc + t(rc))/2
  rownames(rc) = rownames(X)
  colnames(rc) = colnames(X)
  
  attr(rc, "pos.eigen.num") = sum(comp.positive)
  attr(rc, "inv") = rc.inv

  attr(rc, "dof.eff") = apply(rc.e$vectors, 1, function(x) sum(x^2 * (rc.vals != 0)))
  ## attr(rc, "dof.factors") = ifelse(attr(rc, "dof.factors") != 0, 1/sqrt(attr(rc, "dof.factors")), 0)
  names(attr(rc, "dof.eff")) = rownames(X)
  
  return(rc)
}

##
## compute pseudo-inverse square root of symmetric real positive semi-definite matrix
## zero eigenvalues produce zero eigenvalues in the inverse-square matrix
## apply provided weight vector to matrix before inversion in order to reduce computational errors
## the resulting matrix Y satisfies: t(Y) %*% Y = X^-1
##
alu.matr.inv.sqrt.symm.semipos.norm = function(X, X.norm, tol = sqrt(.Machine$double.eps)) {
  if (!is.matrix(X)) 
    X <- as.matrix(X)
  if (length(dim(X)) > 2L || !is.numeric(X))
    stop("'X' must be a numeric matrix")

  if (dim(X)[1] != dim(X)[2])
    stop("'X' must be a square matrix")

  if (!is.numeric(X.norm))
    stop("'X.norm' must be a numeric vector/matrix")

  if (length(X.norm) != dim(X)[1])
    stop("X.norm vector dimension must match X matrix rank")

  XX = diag.m(1/X.norm) %*% X %*% diag.m(1/X.norm)
  rownames(XX) = rownames(X)
  colnames(XX) = colnames(X)
  XX.inv.sqrt.attr = alu.matr.inv.sqrt.symm.semipos(XX, tol)
  XX.inv.sqrt = XX.inv.sqrt.attr %*% diag.m(1/X.norm)
  rownames(XX.inv.sqrt) = rownames(X)
  colnames(XX.inv.sqrt) = colnames(X)
  
  attr(XX.inv.sqrt, "inv") = diag.m(1/X.norm) %*% attr(XX.inv.sqrt.attr, "inv") %*% diag.m(1/X.norm)
  attr(XX.inv.sqrt, "inv") = (attr(XX.inv.sqrt, "inv")+ t(attr(XX.inv.sqrt, "inv")))/2
  rownames(attr(XX.inv.sqrt, "inv")) = rownames(X)
  colnames(attr(XX.inv.sqrt, "inv")) = colnames(X)
  
  attr(XX.inv.sqrt, "pos.eigen.num") = attr(XX.inv.sqrt.attr, "pos.eigen.num")
  attr(XX.inv.sqrt, "dof.eff") = attr(XX.inv.sqrt.attr, "dof.eff")

  return(XX.inv.sqrt)
}
