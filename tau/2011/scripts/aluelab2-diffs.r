#!/usr/bin/env Rscript

##
## aluelab-diffs.r [flags] [<.rdata file>]
## elaborate hfag tau 2011 results
##
## tau 2012 proceedings:
## tau/2001/TauFit> ../scripts/aluelab2-diff.r average2-aleph-hcorr-ref.rdata
##
## - compute fit residual significance
## - print a latex table of the BaBar and Belle significances
##
## the covariance of the fit residuals is computed, diagonalized,
## then the inverse square root is taked nf non-zero eigen-values
## multiplying the fit residuals by this matrix one obtains signed
## fit residual significances. The sum of the squared significances
## is equal to the global fit chi square.
## 

require(stringr, quietly=TRUE)
source("../../../Common/bin/aluelab2.r")
source("../../../Common/bin/alucomb2-hfag-tau.r")

## ////////////////////////////////////////////////////////////////////////////
## functions

##
## get list of quantity names in the definition of a quantity
## using the relevant constraint equation
## used since Oct 2010
## since March 2012, use aluelab.get.str.expr, needs alucomb2 format
##
aluelab.get.quant.names = function(quant.name, combination) {
  var.comb = combination$constr.all.comb[[paste(quant.name, ".c", sep="")]]
  if (is.null(var.comb) || is.na(var.comb)) {
    var.comb = combination$constr.all.comb[[paste(quant.name, ".coq", sep="")]]
  }
  if (is.null(var.comb)) {
    var.comb = combination$constr.lin.comb[[paste(quant.name, ".coq", sep="")]]
  }
  if (is.null(var.comb)) {
    var.comb = combination$constr.all.comb[[quant.name]]
  }
  if (is.null(var.comb)) {
    var.comb = combination$constr.lin.comb[[quant.name]]
  }
  var.comb = var.comb[names(var.comb) != quant.name]
  var.names = names(var.comb[var.comb != 0])
  return(var.names)
}

##
## get string expression for a quantity, from its constraint equations
## needs alucomb2 format
##
aluelab.get.str.expr = function(quant.name, combination) {
  str.expr = combination$constr.all.str.expr[[paste(quant.name, ".c", sep="")]]
  str.expr = str_match(str.expr, "-([[:alnum:]]+)\\s+[+]\\s+(.*\\S)\\s*$")[,3]
  ##--- remove outer braces
  str.expr = gsub("(\\(((?:[^()]++|(?1))+)*+\\))", "\\2", str.expr, perl=TRUE)
  return(str.expr)
}

##
## return tex label such as \Gamma_1 or \frac{\Gamma_1}{\Gamma_2}
##
alucomb2.gamma.texlabel.nv = function(gamma.name) {
  if (gamma.name == "GammaAll") return("1")
  str = str_match(gamma.name, "Gamma(\\d+)(by(\\d+))?")[1,]
  if (is.na(str[1])) return(gamma.name)
  if (str[4] == "") {
    return(paste("\\Gamma_{", str[2], "}", sep=""))
  }
  return(paste("\\frac{\\Gamma_{", str[2], "}}{\\Gamma_{", str[4], "}}", sep=""))
}
alucomb2.gamma.texlabel = Vectorize(alucomb2.gamma.texlabel.nv)

##
##
##
get.tex.quant.descr = function(quant) {
  repeat {
    texdescr = quant$texdescr
    if (!is.null(texdescr) && texdescr != "") {
      ##+++ check for ratio of BR
      texdescr = paste("\\BRF{\\tau^-}{", texdescr, "}", sep="")
      texdescr = gsub("mathrm", "text", texdescr, fixed=TRUE)
      break
    }
    descr = quant$descr
    if (is.null(descr) || descr == "") {
      texdescr = ""
      break
    }
    descr = sub("\\s+/\\s+G\\(total\\)\\s*", "", descr)
    descr = gsub("ex[.]\\s+K\\(S\\)0 --> pi- pi[+]", "ex. K0", descr, perl=TRUE)
    descr = gsub("\\s*\\(``\\d-prong''\\)", "", descr, perl=TRUE)
    descr = gsub("-->", "\\to", descr, fixed=TRUE)
    descr = gsub("([+-])", "^\\1", descr)
    descr = gsub("\\^-prong", "\\\\text{-prong}", descr)
    descr = gsub("G(\\(((?:[^()]++|(?1))+)*+\\))", "\\\\BRF{tau^-}{\\2}", descr, perl=TRUE)
    descr = gsub(">=", "\\\\ge{}", descr)
    descr = gsub("nubar\\((mu|tau|e)\\)", "\\\\bar{nu}_\\1", descr)
    descr = gsub("nu\\((mu|tau|e)\\)", "nu_\\1", descr)
    descr = gsub("pi0", "pi^0", descr)
    descr = gsub("Kstar", "K^*", descr)
    descr = gsub("Kbar0", "\\\\bar{K}^0", descr)
    descr = gsub("K0", "K^0", descr)
    descr = gsub("K\\(S\\)0", "K_S^0", descr)
    descr = gsub("K\\(L\\)0", "K_L^0", descr)
    descr = gsub("(nu|tau|mu|pi|omega|eta|gamma)", "\\\\\\1", descr)
    descr = gsub("\\s+\\(ex[.]", "\\\\;(\\\\text{ex.}", descr)
    descr = gsub("(neutrals)", "\\\\text{\\1}", descr)
    descr = gsub("(particles|strange|total)", "\\\\text{\\1}", descr)
    descr = gsub("(.*\\S)\\s*/\\s*(\\S.*)$", "\\\\frac{\\1}{\\2}", descr, perl=TRUE)
    texdescr = descr
    break
  }

  return(texdescr)
}

##
## return numeric value formatted val +- stat in a string
## use automatic precision and order of magnitude
##
aluelab.fmt = function(val) {
  precision = 4.2
  val.order = ifelse(val == 0, 0, floor(log(abs(val)*1.01)/log(10)))
  if (val.order > 1 && val.order < 4) precision = precision - 0.1*(val.order-1)
  ifelse((val.order <= -3 || val.order >= 4) && FALSE,
         sprintf(paste("\\ensuremath{%", precision, "f\\cdot 10^{%.0f}}", sep=""), val / 10^val.order, val.order),
         sprintf(paste("%", precision, "f", sep=""), val))
}

## ////////////////////////////////////////////////////////////////////////////
## code

aluelab.results = function(args) {
  if(any(args == "-h")) {
    ##--- help
    cat("aluelab-diffs.r [flags] [<.rdata file>]\n")
    cat("  compute deviations of B-factories results\n")
    args = args[args != "-s"]
    stop()
  }

  if (length(args) > 0) {
    file.name = args[1]
  } else {
    file.name = "average.rdata"
  }

  ##--- get alucomb results and data
  load(file.name)

  quant = StatComb$new(quant.val, quant.cov)
  quant.names = names(quant.val)
  comb.params = lapply(combination$params, function(x) unname(x["value"]))
  quant$param.add.list(comb.params)
  quant$param.add.list(c(pi=pi))

  tol = .Machine$double.eps
  
  ##--- pull = measurements minus fitted quantities
  pull.val = drop(meas.val - delta %*% quant.val)
  ##--- covariance of pulls
  pull.cov = meas.cov - delta %*% quant.cov %*% t(delta)
  
  ##--- compute pseudo square root of inverse of a singular matrix
  rc.e = eigen(pull.cov)
  ##--- the non-zero eigen-values are just the ones corresponding to the degrees of freedome of the fit
  eigen.nonzero.num = meas.num - quant.num + constr.num
  ##--- explicitly zero the eigen-values beyond the ones corresponding to real dof
  pull.cov.eigen.val = c(rc.e$values[1:eigen.nonzero.num], rep(0, meas.num - eigen.nonzero.num))
 
  ##--- invert just the non-zero eigen values, then take square root
  pull.cov.eigen.val.inv.sqrt = ifelse( pull.cov.eigen.val != 0, 1/sqrt(pull.cov.eigen.val), 0)
  ##--- get inverse square root of covariance for measurements pulls
  pull.cov.inv.sqrt = rc.e$vectors %*% diag(pull.cov.eigen.val.inv.sqrt) %*% t(rc.e$vectors)
  ##--- compute significance of measuremnts pulls
  pull.chisq = drop(pull.val %*% pull.cov.inv.sqrt)
  names(pull.chisq) = names(meas.val)
  
  ##--- compute equivalent number of dof for measurements pulls
  pull.cov.eigen.nonzero = ifelse(pull.cov.eigen.val != 0, 1, 0)
  pull.dof = drop(rc.e$vectors^2 %*% pull.cov.eigen.nonzero)
  names(pull.dof) = names(meas.val)
  ##--- zero dof numerically compatible with zero
  pull.dof = ifelse(pull.dof > length(pull.cov.eigen.nonzero)^2*tol, pull.dof, 0)

  ##--- explicitly set to zero the chisq contribution of pulls for single measurements with zero ndof
  pull.chisq = ifelse(pull.dof != 0, pull.chisq, 0)
  
  ##--- normalize pull significance with square root of effective dof (probably wrong)
  pull.signif.a = ifelse(pull.dof != 0, pull.chisq/sqrt(pull.dof), 0)

  ##--- get pulls in the base that diagonalizes the covariance matrix
  pull.diag.val = drop(t(rc.e$vectors) %*% pull.val)
  ##--- significance in the diag. base (set to 1 if the eigenvalue is zero)
  pull.diag.signif2 = ifelse(pull.cov.eigen.val != 0, pull.diag.val^2 / pull.cov.eigen.val, 1)
  ##--- pull significances square are sum of diag base significances square times coeff square
  pull.signif2.b = drop(rc.e$vectors^2 %*% pull.diag.signif2)
  ##--- explicitly set to zero square signif of pulls for single measurements with zero ndof
  pull.signif2.b = ifelse(pull.dof != 0, pull.signif2.b, 0)
  ##--- use sign of chisq residuals to sign the pull significances
  pull.signif.b = sign(pull.chisq)*sqrt(pull.signif2.b)

  ##--- pulls weighted with inverse sqrt of measurements covariance
  rc.e = eigen(meas.cov)
  meas.cov.inv.sqrt = rc.e$vectors %*% diag(1/sqrt(rc.e$values)) %*% t(rc.e$vectors)
  pull.signif.c = drop(pull.val %*% meas.cov.inv.sqrt)
  pull.signif.c = ifelse(pull.dof != 0, pull.signif.c, 0)
  names(pull.signif.c) = names(meas.val)

  ##--- choose signif definition
  pull.signif = pull.signif.c
  
  ##--- reorder by abs value of pulls significance
  pull.signif.order = order(abs(pull.signif))
  pull.signif.sorted = pull.signif[pull.signif.order]
  pull.val.sorted = pull.val[pull.signif.order]
  meas.err.sorted = meas.err[pull.signif.order]
  pull.dof.sorted = pull.dof[pull.signif.order]

  ## print(cbind(pull.val.sorted, meas.err.sorted, pull.signif.sorted, pull.dof.sorted))

  ##--- get fit-quantity name corresponding to a measurement
  meas.to.quant.name = function(name) {
    ## sub("^[^.]+[.]([^.]+)[.].*", "\\1", name, perl=TRUE)
    quant.names[apply(delta[name, , drop=FALSE], 1, function(x) which(x != 0))]
  }

  meas.names = names(meas.val)

  ##--- select BaBar measurements
  meas.babar = meas.names[grep("^BaBar[.]", names(meas.val), perl=TRUE)]
  ##--- remove non overconstrained measurements
  meas.babar = setdiff(meas.babar, meas.names[pull.dof == 0])
  names(meas.babar) = meas.to.quant.name(meas.babar)

  ##--- select Belle measurements
  meas.belle = meas.names[grep("^Belle[.]", names(meas.val), perl=TRUE)]
  ##--- remove non overconstrained measurements
  meas.belle = setdiff(meas.belle, meas.names[pull.dof == 0])
  names(meas.belle) = meas.to.quant.name(meas.belle)

  ##--- prepare ordered listing by quantity
  quant.bfactories = unique(c(names(meas.babar), names(meas.belle)))
  quant.bfactories = quant.bfactories[order(alucomb2.gamma.num.id(quant.bfactories))]

  cat("\\begin{tabular}{lrr}\n")
  cat("\\toprule\n")
  cat("tau branching fraction\ & \\multicolumn{1}{c}{\\babar} & \\multicolumn{1}{c}{Belle} \\\\\n")
  cat("\\midrule\n")
  pull.check.meas = character()
  for (quant.name in quant.bfactories) {
    ## cat(alucomb2.gamma.texlabel(quant.name),"\n")
    meas.babar.quant = meas.babar[quant.name == names(meas.babar)]
    meas.belle.quant = meas.belle[quant.name == names(meas.belle)]
    if (length(meas.babar.quant)>1 || length(meas.belle.quant)>1) {
      stop("more than one measurement for", quant.name)
    }
    if (length(meas.babar.quant)==1) pull.check.meas = c(pull.check.meas, meas.babar.quant)
    if (length(meas.belle.quant)==1) pull.check.meas = c(pull.check.meas, meas.belle.quant)
    cat("\\ensuremath{",
        get.tex.quant.descr(combination$quantities[[quant.name]]),
        "} & ",
        ifelse(length(meas.babar.quant)==1, aluelab.fmt(pull.signif[meas.babar.quant]), ""),
        " & ",
        ifelse(length(meas.belle.quant)==1, aluelab.fmt(pull.signif[meas.belle.quant]), ""),
        " \\\\\n", sep="")
  }
  cat("\\bottomrule\n")
  cat("\\end{tabular}\n")
  
  rc = print(cbind(
    pull=pull.val[pull.check.meas],
    meas.err=meas.err[pull.check.meas],
    pull.chisq=pull.chisq[pull.check.meas],
    signif.a=pull.signif.a[pull.check.meas],
    signif.b=pull.signif.b[pull.check.meas],
    signif.c=pull.signif.c[pull.check.meas],
    eff.dof=pull.dof[pull.check.meas]))
}

args = commandArgs(TRUE)
aluelab.results(args)
