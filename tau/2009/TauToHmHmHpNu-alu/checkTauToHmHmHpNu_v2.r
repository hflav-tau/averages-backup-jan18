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
## corr = matrix(0, meas.num, meas.num)
corr = diag(rep(1,meas.num))
rownames(corr) = names(measurements)
colnames(corr) = names(measurements)

##--- build off-diagonal correlation matrix
for (meas in measurements) {
  mapply(function(other.tag, other.corr) {
    corr[meas@tag,other.tag] <<- other.corr
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
    if (corr[meas.i,meas.j] != 0) {
      cat("corr[\"",meas.i,"\",\"",meas.j,"\"] = ",corr[meas.i,meas.j],"\n",sep="")
    }
  }
}
}

if (!kLocalCards) {
##
## replace correlation coefficients computed by Swagato with the original ones
## which are total correlation coefficients
##
corr = 0*corr
corr["Belle.PimPimPipNu.published","Belle.PimKmPipNu.published"] = 0.1749885
corr["Belle.PimPimPipNu.published","Belle.PimKmKpNu.published"] = 0.04948276
corr["Belle.PimPimPipNu.published","Belle.KmKmKpNu.published"] = -0.05346557
corr["Belle.PimKmPipNu.published","Belle.PimKmKpNu.published"] = 0.08026913
corr["Belle.PimKmPipNu.published","Belle.KmKmKpNu.published"] = 0.03505142
corr["Belle.PimKmKpNu.published","Belle.KmKmKpNu.published"] = -0.008312885
corr["BaBar.PimPimPipNu.published","BaBar.PimKmPipNu.published"] = 0.543535
corr["BaBar.PimPimPipNu.published","BaBar.PimKmKpNu.published"] = 0.390346
corr["BaBar.PimPimPipNu.published","BaBar.KmKmKpNu.published"] = 0.031469
corr["BaBar.PimKmPipNu.published","BaBar.PimKmKpNu.published"] = 0.177495
corr["BaBar.PimKmPipNu.published","BaBar.KmKmKpNu.published"] = 0.0931907
corr["BaBar.PimKmKpNu.published","BaBar.KmKmKpNu.published"] = 0.0870484
##--- symmetrize (we only entered the upper part)
corr = corr + t(corr)
##--- add diagonal unity matrix
corr = corr + diag(rep(1,meas.num))
}

##
## get list of stat, syst errors, compute total error
##
meas.stat = unlist(lapply(measurements, function(x) x@stat))
meas.syst = unlist(lapply(measurements, function(x) x@syst))
meas.error = sqrt(meas.stat^2 + meas.syst^2)

##--- build covariance matrix using errors and correlation coefficients
covariance = corr * (meas.error %o% meas.error)

##--- quantities we want to determinee from measurements
quant.names = c("PimPimPipNu", "PimKmPipNu", "PimKmKpNu", "KmKmKpNu")
quant.num = length(quant.names)

delta = matrix(0, quant.num, meas.num)
rownames(delta) = quant.names
colnames(delta) = meas.names

##
## temporary code to print the non-zero delta matrix coefficients
##
if (FALSE) {
for (quant in quant.names) {
  for (meas in meas.names) {
    if (regexpr(quant, meas, fixed=TRUE) != -1) {
      cat("delta[\"",quant,"\",\"",meas,"\"] = ","1\n",sep="")
    }
    if (regexpr("HmHmHpNu", meas, fixed=TRUE) != -1) {
      cat("delta[\"",quant,"\",\"",meas,"\"] = ","1\n",sep="")
    }
  }
}
}

##
## build delta matrix
## measurements are linear combinations of quantities we want to average
## for measurements of a quantity, a -1 must be set in the slot that
## relates the quantity with the measurement
## for measurements that are the sum of some (or all) the quantities,
## a -1 must be set in all the slots that relates the quantity with the measurements
##
delta["PimPimPipNu","Belle.PimPimPipNu.published"] = -1
delta["PimPimPipNu","BaBar.PimPimPipNu.published"] = -1
delta["PimPimPipNu","CLEO3.PimPimPipNu.published"] = -1
delta["PimPimPipNu","DELPHI.HmHmHpNu.published"] = -1
delta["PimPimPipNu","OPAL.HmHmHpNu.published"] = -1
delta["PimPimPipNu","CLEO.HmHmHpNu.published"] = -1
delta["PimKmPipNu","Belle.PimKmPipNu.published"] = -1
delta["PimKmPipNu","BaBar.PimKmPipNu.published"] = -1
delta["PimKmPipNu","OPAL.PimKmPipNu.published"] = -1
delta["PimKmPipNu","CLEO3.PimKmPipNu.published"] = -1
delta["PimKmPipNu","CLEO.PimKmPipNu.published"] = -1
delta["PimKmPipNu","ALEPH.PimKmPipNu.published"] = -1
delta["PimKmPipNu","DELPHI.HmHmHpNu.published"] = -1
delta["PimKmPipNu","OPAL.HmHmHpNu.published"] = -1
delta["PimKmPipNu","CLEO.HmHmHpNu.published"] = -1
delta["PimKmKpNu","Belle.PimKmKpNu.published"] = -1
delta["PimKmKpNu","BaBar.PimKmKpNu.published"] = -1
delta["PimKmKpNu","CLEO3.PimKmKpNu.published"] = -1
delta["PimKmKpNu","OPAL.PimKmKpNu.published"] = -1
delta["PimKmKpNu","CLEO.PimKmKpNu.published"] = -1
delta["PimKmKpNu","AlEPH.PimKmKpNu.published"] = -1
delta["PimKmKpNu","DELPHI.HmHmHpNu.published"] = -1
delta["PimKmKpNu","OPAL.HmHmHpNu.published"] = -1
delta["PimKmKpNu","CLEO.HmHmHpNu.published"] = -1
delta["KmKmKpNu","Belle.KmKmKpNu.published"] = -1
delta["KmKmKpNu","BaBar.KmKmKpNu.published"] = -1
delta["KmKmKpNu","DELPHI.HmHmHpNu.published"] = -1
delta["KmKmKpNu","OPAL.HmHmHpNu.published"] = -1
delta["KmKmKpNu","CLEO.HmHmHpNu.published"] = -1

##--- measurement values
meas = unlist(lapply(measurements, function(x) {x@value}))

##
## solve for quantities
##

invcov = solve(covariance)

covariance.quant = solve(delta %*% invcov %*% t(delta))
rownames(covariance.quant) = quant.names
colnames(covariance.quant) = quant.names
quant.err = sqrt(diag(covariance.quant))

quant = drop(-covariance.quant %*% delta %*% (invcov %*% meas))
names(quant) = quant.names

corr.quant = covariance.quant / (quant.err %o% quant.err)

chisq = t(meas + t(delta) %*% quant) %*% invcov %*% (meas + t(delta) %*% quant)

cat("##\n")
cat("## exact solution, chi-square = ",chisq,"\n")
cat("##\n")
cat("averages\n")
show(quant)
cat("errors\n")
show(quant.err)
## cat("covariance\n")
## show(covariance.quant)
cat("correlation\n")
show(corr.quant)

##
## solve for quantities with iterative chi-square minimization
##

cat("##\n")
cat("## minimum shi-square fit\n")
cat("##\n")

logLik.average = function(par) {
  chisq = t((meas + t(delta)) %*% par) %*% invcov %*% ((meas + t(delta)) %*% par)
  return(-1/2*chisq)
}

fit = maxLik(logLik.average, start=quant*1.5)

show(fit)
covariance.quant = vcov(fit)
quant.err = sqrt(diag(covariance.quant))
corr.quant = covariance.quant / (quant.err %o% quant.err)

chisq = t(meas + t(delta) %*% quant) %*% invcov %*% (meas + t(delta) %*% quant)

cat("##\n")
cat("## minimum chi square fit solution = ",chisq,"\n")
cat("##\n")
cat("averages\n")
show(quant)
cat("errors\n")
show(quant.err)
## cat("covariance\n")
## show(covariance.quant)
cat("correlation\n")
show(corr.quant)

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
cat("covQuant = Inverse[delta invcov Transpose[delta]];\n")
cat("errQuant = Sqrt[Diag[covQuant]];\n")
cat("quant = -covQuant delta (invcov meas);\n")
}
