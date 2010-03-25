#!/usr/bin/env Rscript

##
## test 2nd average of tau -> hhh nu by Swagato
## - including statistical correlations between 4 results in Combos
## - including tau -> hhh nu measurements from PDG
## - takes all measurement values and errors from ../TauToHmHmHpNu/ input cards
## - uses total correlation coefficients from ../TauToHmHmHpNu/ .pl files
##

library(methods)
library(maxLik)

## ////////////////////////////////////////////////////////////////////////////
## definitions

##
## return list of matching strings or NULL is no match
##
get.matches <- function(expr, text, ...) {
  match.pos.list = gregexpr(expr, text, ...)
  matches.list = list()
  temp = mapply(function(line, match.pos) {
    if (match.pos[1] != -1) {
      matches = mapply(function(pos, count) {
        substr(line, pos, pos+count-1)
      }, match.pos, attr(match.pos, "match.length"))
    } else {
      matches = NULL
    }
    return(matches)
  }, text, match.pos.list, USE.NAMES=FALSE)
  return(temp)
}

##
## return list of lines from file
##
get.file.lines <- function(fname) {
  fh <- file(fname)
  lines  <- readLines(fh)
  close(fh)
  return(lines)
}

##
## class for a measurement external parameter
##
rc = setClass("parameter",
  representation(value = "numeric",
                 error = "numeric"
                 ),
  prototype(value = numeric(0),
            error = numeric(0)
            )
  )

##
## class for a measurement
##
rc = setClass("measurement",
  representation(value = "numeric",
                 stat = "numeric",
                 syst = "numeric",
                 bibitem = "character",
                 tag = "character",
                 params = "list",
                 syst.terms = "list",
                 correlations = "list"
                 ),
  prototype(value = numeric(0),
            stat = numeric(0),
            syst = numeric(0)
            )
  )

## ////////////////////////////////////////////////////////////////////////////
## code

kLocalCards = TRUE

lines = NULL
if (kLocalCards) {
  lines = c(lines, get.file.lines("belle.input"))
  lines = c(lines, get.file.lines("babar.input"))
  lines = c(lines, get.file.lines("pdg.input"))
  lines = c(lines, get.file.lines("pdg_hhh.input"))
} else {
  lines = c(lines, get.file.lines("../TauToHmHmHpNu/belle.input"))
  lines = c(lines, get.file.lines("../TauToHmHmHpNu/babar.input"))
  lines = c(lines, get.file.lines("../TauToHmHmHpNu/pdg.input"))
  lines = c(lines, get.file.lines("../TauToHmHmHpNu/pdg_hhh.input"))
}
  
measurements = list()

kIdle = 1
kBegin = 2
kMeasurement = 3
kMeasData = 4

status = kIdle
for (line in lines) {
  if (regexpr("^\\s*$", line) != -1 ||
      regexpr("^[*cC#;]", line) != -1) {
    next
  }
  fields = get.matches("(\\S+)", line)
  if (regexpr("^BEGIN$", fields[1], ignore.case=TRUE) != -1) {
    meas = new("measurement")
    meas@bibitem = as.character(fields[-1])
    meas@tag = paste(as.character(fields[2:4]), collapse=".")
    status = kBegin
    ## cat(line,"\n")
  } else if (status == kBegin && regexpr("^MEASUREMENT$", fields[1], ignore.case=TRUE) != -1) {
    status = kMeasurement
    ## cat("measurement\n")
  } else if (status == kMeasurement && regexpr("^DATA$", fields[1], ignore.case=TRUE) != -1) {
    status = kMeasData
    ## cat("measurement data\n")
  } else if (status == kMeasData) {
    warning = FALSE
    withCallingHandlers({fields.num <- as.numeric(fields)}, warning = function(w) {
      warning <<- TRUE; invokeRestart("muffleWarning")})
    if (warning) {
      cat("Cannot interpret the following data as a measurement value, stat, syst\n")
      cat("  '",line,"'\n")
    } else {
      meas@value = fields.num[1]
      meas@stat = fields.num[2]
      meas@syst = fields.num[3]
    }
    status = kBegin
  } else if (regexpr("^STAT_CORR_WITH$", fields[1], ignore.case=TRUE) != -1) {
    corr = as.numeric(fields[5])
    names(corr) =  paste(as.character(fields[2:4]), collapse=".")
    meas@correlations = c(meas@correlations,corr)
  } else if (regexpr("^END$", fields[1], ignore.case=TRUE) != -1) {
    status = kIdle
    measurements = c(measurements, meas)
  }
}

##-- assign names to measurements array elements equal to their tags
names(measurements) = unlist(lapply(measurements, function(x) x@tag))
meas.num = length(measurements)
meas.names = names(measurements)
## cat(names(measurements),"\n")

##-- init correlation matrix
## meas.corr = matrix(0, meas.num, meas.num)
meas.corr = diag(rep(1,meas.num))
rownames(meas.corr) = names(measurements)
colnames(meas.corr) = names(measurements)

##--- build off-diagonal correlation matrix
for (meas in measurements) {
  mapply(function(other.tag, other.corr) {
    meas.corr[meas@tag,other.tag] <<- other.corr
    ## cat("other = ", other.corr, "\n")
  }, names(meas@correlations), meas@correlations)
}

##
## temporary code to print the non zero off-diagonal elements of the correlation matrix
##
if (FALSE) {
index = 0
for (meas.i in meas.names) {
  index = index+1
  if (index+1 > length(meas.names)) next
  ## cat(meas.i,"\n")
  ## cat(meas.names[(index+1):length(meas.names)],"\n")
  for (meas.j in meas.names[(index+1):length(meas.names)]) {
    if (meas.corr[meas.i,meas.j] != 0) {
      cat("meas.corr[\"",meas.i,"\",\"",meas.j,"\"] = ",meas.corr[meas.i,meas.j],"\n",sep="")
    }
  }
}
}

if (!kLocalCards) {
##
## replace correlation coefficients computed by Swagato with the original ones
## which are total correlation coefficients
##
meas.corr = 0*meas.corr
meas.corr["Belle.PimPimPipNu.published","Belle.PimKmPipNu.published"] = 0.1749885
meas.corr["Belle.PimPimPipNu.published","Belle.PimKmKpNu.published"] = 0.04948276
meas.corr["Belle.PimPimPipNu.published","Belle.KmKmKpNu.published"] = -0.05346557
meas.corr["Belle.PimKmPipNu.published","Belle.PimKmKpNu.published"] = 0.08026913
meas.corr["Belle.PimKmPipNu.published","Belle.KmKmKpNu.published"] = 0.03505142
meas.corr["Belle.PimKmKpNu.published","Belle.KmKmKpNu.published"] = -0.008312885
meas.corr["BaBar.PimPimPipNu.published","BaBar.PimKmPipNu.published"] = 0.543535
meas.corr["BaBar.PimPimPipNu.published","BaBar.PimKmKpNu.published"] = 0.390346
meas.corr["BaBar.PimPimPipNu.published","BaBar.KmKmKpNu.published"] = 0.031469
meas.corr["BaBar.PimKmPipNu.published","BaBar.PimKmKpNu.published"] = 0.177495
meas.corr["BaBar.PimKmPipNu.published","BaBar.KmKmKpNu.published"] = 0.0931907
meas.corr["BaBar.PimKmKpNu.published","BaBar.KmKmKpNu.published"] = 0.0870484
##--- symmetrize (we only entered the upper part)
meas.corr = meas.corr + t(meas.corr)
##--- add diagonal unity matrix
meas.corr = meas.corr + diag(rep(1,meas.num))
}

##
## get list of stat, syst errors, compute total error
##
meas.stat = unlist(lapply(measurements, function(x) x@stat))
meas.syst = unlist(lapply(measurements, function(x) x@syst))
meas.error = sqrt(meas.stat^2 + meas.syst^2)

##--- build covariance matrix using errors and correlation coefficients
meas.cov = meas.corr * (meas.error %o% meas.error)

##--- quantities we want to determinee from measurements
quant.names = c("PimPimPipNu", "PimKmPipNu", "PimKmKpNu", "KmKmKpNu")
quant.num = length(quant.names)

delta = matrix(0, meas.num, quant.num)
colnames(delta) = quant.names
rownames(delta) = meas.names

##
## temporary code to print the non-zero delta matrix coefficients
##
if (FALSE) {
for (meas in meas.names) {
  for (quant in quant.names) {
    if (regexpr(quant, meas, fixed=TRUE) != -1) {
      cat("delta[\"",meas,"\",\"",quant,"\"] = 1\n",sep="")
    }
    if (regexpr("HmHmHpNu", meas, fixed=TRUE) != -1) {
      cat("delta[\"",meas,"\",\"",quant,"\"] = 1\n",sep="")
    }
  }
}
}

##
## build delta matrix
## - measurements are experimental results or external PDG averages
## - quantities are the results of the HFAG averaging procedure
## measurements are linear combinations of quantities, meas_i = delta_ij * quant_j
## if a measurement i corresponds to a quantity j then delta_ij = 1
## some measurements are actually a sum of quantities: meas_i = quant_j1 + quant_j2 + ...
## in this case delta_i,j1 = 1, delta_i,j2 = 1, ...
## all remaining delta matrix terms are zero
## one can generalize the above concepts to measurements that are linear combinations
## of quantities by using proper coefficients different from 1
##
delta["Belle.PimPimPipNu.published","PimPimPipNu"] = 1
delta["Belle.PimKmPipNu.published","PimKmPipNu"] = 1
delta["Belle.PimKmKpNu.published","PimKmKpNu"] = 1
delta["Belle.KmKmKpNu.published","KmKmKpNu"] = 1
delta["BaBar.PimPimPipNu.published","PimPimPipNu"] = 1
delta["BaBar.PimKmPipNu.published","PimKmPipNu"] = 1
delta["BaBar.PimKmKpNu.published","PimKmKpNu"] = 1
delta["BaBar.KmKmKpNu.published","KmKmKpNu"] = 1
delta["CLEO3.PimPimPipNu.published","PimPimPipNu"] = 1
delta["OPAL.PimKmPipNu.published","PimKmPipNu"] = 1
delta["CLEO3.PimKmPipNu.published","PimKmPipNu"] = 1
delta["CLEO.PimKmPipNu.published","PimKmPipNu"] = 1
delta["ALEPH.PimKmPipNu.published","PimKmPipNu"] = 1
delta["CLEO3.PimKmKpNu.published","PimKmKpNu"] = 1
delta["OPAL.PimKmKpNu.published","PimKmKpNu"] = 1
delta["CLEO.PimKmKpNu.published","PimKmKpNu"] = 1
delta["AlEPH.PimKmKpNu.published","PimKmKpNu"] = 1
delta["DELPHI.HmHmHpNu.published","PimPimPipNu"] = 1
delta["DELPHI.HmHmHpNu.published","PimKmPipNu"] = 1
delta["DELPHI.HmHmHpNu.published","PimKmKpNu"] = 1
delta["DELPHI.HmHmHpNu.published","KmKmKpNu"] = 1
delta["OPAL.HmHmHpNu.published","PimPimPipNu"] = 1
delta["OPAL.HmHmHpNu.published","PimKmPipNu"] = 1
delta["OPAL.HmHmHpNu.published","PimKmKpNu"] = 1
delta["OPAL.HmHmHpNu.published","KmKmKpNu"] = 1
delta["CLEO.HmHmHpNu.published","PimPimPipNu"] = 1
delta["CLEO.HmHmHpNu.published","PimKmPipNu"] = 1
delta["CLEO.HmHmHpNu.published","PimKmKpNu"] = 1
delta["CLEO.HmHmHpNu.published","KmKmKpNu"] = 1

##--- measurement values
meas = unlist(lapply(measurements, function(x) {x@value}))

##
## solve for quantities
##

invcov = solve(meas.cov)

quant.cov = solve(t(delta) %*% invcov %*% delta)
rownames(quant.cov) = quant.names
colnames(quant.cov) = quant.names
quant.err = sqrt(diag(quant.cov))

quant = drop(quant.cov %*% t(delta) %*% (invcov %*% meas))
names(quant) = quant.names

quant.corr = quant.cov / (quant.err %o% quant.err)

chisq = t(meas - delta %*% quant) %*% invcov %*% (meas - delta %*% quant)

cat("##\n")
cat("## exact solution, chi-square = ",chisq,"\n")
cat("##\n")
cat("averages\n")
show(quant)
cat("errors\n")
show(quant.err)
## cat("covariance\n")
## show(quant.cov)
cat("correlation\n")
show(quant.corr)

##
## solve for quantities with iterative chi-square minimization
##

logLik.average = function(par) {
  chisq = t(meas - delta %*% par) %*% invcov %*% (meas - delta %*% par)
  return(-1/2*chisq)
}

fit = maxLik(logLik.average, start=quant*1.5)

quant = coef(fit)
quant.cov = vcov(fit)
quant.err = sqrt(diag(quant.cov))
quant.corr = quant.cov / (quant.err %o% quant.err)

chisq = t(meas - delta %*% quant) %*% invcov %*% (meas - delta %*% quant)
chisq.fit = -2*logLik(fit)

cat("##\n")
cat("## minimum chi square fit, chi-square= ", chisq, ", from fit= ", chisq.fit, "\n",sep="")
cat("##\n")
cat("averages\n")
show(quant)
cat("errors\n")
show(quant.err)
## cat("covariance\n")
## show(quant.cov)
cat("correlation\n")
show(quant.corr)
cat("\n")
cat("## begin fit summary\n")
show(fit)
cat("## end fit summary\n")

if (FALSE) {
##
## output Mathematica code to solve the same problem
##

matrix.to.math = function( matr ) {
  cov.out = character(0)
  for(row in 1:dim(matr)[1]) {
    cov.out = c(cov.out, paste("{",paste(matr[row,],collapse=","),"}",sep=""))
  }
  return(paste("{",paste(cov.out,collapse=","),"}", sep=""))
}

vector.to.math = function( vect ) {
  return(paste("{",paste(vect,collapse=","),"}",sep=""))
}

cat("cov = ", matrix.to.math(covariance), ";\n", sep="")
cat("delta = ", matrix.to.math(delta), ";\n", sep="")
cat("meas = ", vector.to.math(meas), ";\n", sep="")

cat("invcov = Inverse[cov];\n")
cat("covQuant = Inverse[Transpose[delta] invcov delta];\n")
cat("errQuant = Sqrt[Diag[covQuant]];\n")
cat("quant = -covQuant Transpose[delta] (invcov meas);\n")
}
