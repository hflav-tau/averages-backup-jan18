#!/usr/bin/env Rscript

##
## alucomb
##
## Copyright Alberto Lusiani 2010. All rights reserved.
## an open source license will be set once tested and working
##
## - averages measurements in a way that is mostly compatible with Combos
## - can read Combos input files
## - can average multiple quantities related to multiple statistically correlated measurements
##

library(methods)

flag.no.maxLik = FALSE
rc = try(library(maxLik))
if (inherits(rc, "try-error")) flag.no.maxLik = TRUE

## ////////////////////////////////////////////////////////////////////////////
## definitions

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
## test if pattern matches string irrespective of letters case
##
match.nocase = function(pattern, str) {
  return(regexpr(pattern, str, ignore.case=TRUE) != -1)
}
  
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
                 syst.terms = "numeric",
                 corr.terms = "list",
                 corr.terms.tot = "list"
                 ),
  prototype(value = numeric(0),
            stat = numeric(0),
            syst = numeric(0),
            params = list(),
            syst.terms = numeric(0)
            )
  )

##
## class for a combination
##
rc = setClass("combination",
  representation(value = "numeric",
                 error = "numeric",
                 bibitem = "character",
                 tag = "character",
                 quantities = "character",
                 params = "list",
                 meas.lin.combs = "list"
                 ),
  prototype(value = numeric(0),
            error = numeric(0),
            quantities = character(0),
            params = list(),
            meas.lin.combs = list()
            )
  )

## ////////////////////////////////////////////////////////////////////////////
## code

file = "average.input"

##++alucomb = function(file = "") {

if (match.nocase("^\\s*$", file)) {
  stop("alucomb: please provide as argument the card input file\n")
}

flag.replace.corr = FALSE
flag.build.delta = FALSE

if (!file.exists(file)) {
  stop("cannot find file ", file, "\n")
}
dir.base = dirname(file)
lines = get.file.lines(file)
cat("read file", file, "\n")

iline = 1
repeat {
  lines.len = length(lines)
  if (iline > lines.len) break
  fields = unlist(strsplit(lines[iline], "\\s+", perl=TRUE))
  if (length(fields) > 0 &&
      match.nocase("^INCLUDE$", fields[1])) {
    file.inc = paste(dir.base, fields[2], sep="/")
    if (!file.exists(file.inc)) {
      stop("cannot find included file ", file.inc, "\n")
    }
    lines.inc = get.file.lines(file.inc)
    cat("read file", file.inc, "\n")
    lines = c(
      if (iline-1>=1) {lines[1:(iline-1)]} else {NULL},
      lines.inc,
      if (iline+1<=lines.len) {lines[(iline+1):(lines.len)]} else {NULL}
      )
    next
  }
  iline = iline + 1
}

##
## define storage for measurements and one combination
##
measurements = list()
combination = new("combination")

flag.in.meas = FALSE
flag.in.data = FALSE
flag.in.params = FALSE
flag.in.combine = FALSE
flag.in.sumofmeas = FALSE

for (line in lines) {
  ## cat(line,"\n")
  if (regexpr("^\\s*$", line, perl=TRUE) != -1 ||
      regexpr("^[*#;]", line, perl=TRUE) != -1) {
    next
  }
  line = gsub("\\s*!.*", "", line, perl=TRUE)
  fields = unlist(strsplit(line, "\\s+", perl=TRUE))
  if (match.nocase("^BEGIN$", fields[1])) {
    if (flag.in.meas) {
      cat("error, BEGIN keyword inside definition of", meas@tag, "\n")
    } else {
      meas = new("measurement")
      meas@bibitem = as.character(fields[-1])
      meas@tag = paste(as.character(fields[2:4]), collapse=".")
      meas.labels = character(0)
      data.labels = character(0)
      data.values = numeric(0)
      sumofmeas.values = numeric(0)
      measlincombs.list = list()
      flag.in.meas = TRUE
    }
    next
  }
  if (!flag.in.meas) {
    cat("error, ", fields[1], "outside a measurement definition (BEGIN..END)\n")
    next
  }
  if (!match.nocase("^\\s*$", fields[1])) {
    if (flag.in.sumofmeas) {
      ##-- vector with one for each quantity whose sum corresponds to the measurement
      val = rep(1,length(sumofmeas.values[-(1:3)]))
      names(val) = sumofmeas.values[-(1:3)]
      ##-- list of measurements with the coefficients corresponding to the related quantities
      measlincombs.list = c(measlincombs.list, list(val))
      ##-- the first field corresponds to the measurement tag
      names(measlincombs.list)[length(measlincombs.list)] = paste(sumofmeas.values[1:3],collapse=".")
      ##-- reset vector of SUMOFMEAS parameter
      sumofmeas.values = numeric(0)
    }
    flag.in.data = FALSE
    flag.in.params = FALSE
    flag.in.sumofmeas = FALSE
  }
  if (match.nocase("^END$", fields[1])) {
    if (!flag.in.combine) {
      ##
      ## measurement cards
      ##
      names(data.values) = data.labels
      if (length(meas.labels) != 3) {
        stop("wrong number of MEASUREMENT labels: ",length(meas.labels),"instead of 3\n")
      }
      for (i in 1:3) {
        if (is.na(data.values[meas.labels[1]])) {
          stop("missing MEASUREMENT data for ", meas.labels[1], "\n")
        }
      }
      meas.labels[1] = sub("^m_", "", meas.labels[1], ignore.case=TRUE)
      data.labels[1] = sub("^m_", "", data.labels[1], ignore.case=TRUE)
      if (meas@bibitem[2] != data.labels[1]) {
        ##-- when combining multiple quantities, methos should be set as the quantity
        cat("warning: measurement method '", meas@bibitem[2], "' does not match value '", data.labels[1], "'\n", sep="")
      }
      names(data.values)[1] = data.labels[1]
      meas@value = data.values[meas.labels[1]]
      patt.perc = "([[:alnum:]]+[^[:alnum:]]*)(%)([^[:alnum:]]*)$"
      for (i in 2:length(data.values)) {
        if (regexpr(patt.perc, data.labels[i]) != -1) {
          data.labels[i] = gsub(patt.perc, "\\1\\3", data.labels[i])
          data.values[i] = meas@value * data.values[i] /100
        }
      }
      names(data.values)[-1] = data.labels[-1]
      meas@stat = data.values[meas.labels[2]]
      names(meas@stat) = sub("statistical", "stat", names(meas@stat), ignore.case=TRUE)
      meas@syst = data.values[meas.labels[3]]
      names(meas@syst) = sub("systematic", "syst", names(meas@syst), ignore.case=TRUE)
      meas@syst.terms = data.values[!data.labels %in% meas.labels]
      measurements = c(measurements, meas)
    } else {
      ##
      ## combination cards
      ##
      combination@tag = meas@tag
      combination@bibitem = meas@bibitem
      combination@quantities = meas.labels
      combination@params = meas@params
      combination@meas.lin.combs = measlincombs.list
    }
    flag.in.meas = FALSE
    flag.in.combine = FALSE
    next
  }
  if (match.nocase("^MEASUREMENT$", fields[1])) {
    if (!flag.in.combine) {
      ##-- get labels for value, stat.error, syst.error
      meas.labels = fields[2:4]
    } else {
      ##-- get list of measurements to combine
      meas.labels = c(meas.labels, sub("^m_","",fields[2]))
    }
    next
  }
  if (match.nocase("^COMBINE$", fields[1])) {
    flag.in.combine = TRUE
    next
  }
  if (match.nocase("^DATA$", fields[1]) ||
      (flag.in.data && match.nocase("^\\s*$", fields[1]))) {
    flag.in.data = TRUE
    data.labels = c(data.labels,
      unlist(lapply(fields[-1], function(elem) {if (is.na(suppressWarnings(as.numeric(elem)))) {elem}})))
    data.values = c(data.values,
      unlist(lapply(fields[-1], function(elem) {if (!is.na(suppressWarnings(as.numeric(elem)))) {as.numeric(elem)}})))
    next
  }
  if (match.nocase("^STAT_CORR_WITH$", fields[1])) {
    corr = as.numeric(fields[5])
    names(corr) =  paste(as.character(fields[2:4]), collapse=".")
    meas@corr.terms = c(meas@corr.terms,corr)
    next
  }
  if (match.nocase("^ERROR_CORR_WITH$", fields[1])) {
    corr = as.numeric(fields[5])
    names(corr) =  paste(as.character(fields[2:4]), collapse=".")
    meas@corr.terms.tot = c(meas@corr.terms.tot, corr)
    next
  }
  if (match.nocase("^PARAMETERS$", fields[1]) ||
      (flag.in.params && match.nocase("^\\s*$", fields[1]))) {
    flag.in.params = TRUE
    params.data = lapply(fields[-1], function(x) type.convert(x, as.is=TRUE))
    if (length(params.data) == 0) next
    if (!is.character(params.data[[1]]) ||
        !is.numeric(params.data[[2]]) ||
        !is.numeric(params.data[[3]]) ||
        !is.numeric(params.data[[4]])) {
      stop("wrong parameter data ", paste(params.data,sep=","), "\n")
    }
    names(params.data)[2:4] = c("value", "delta_pos", "delta_neg")
    meas@params = c(meas@params, list(unlist(params.data[2:4])))
    names(meas@params)[length(meas@params)] = params.data[[1]]
    next
  }
  if (match.nocase("^SUMOFMEAS$", fields[1]) ||
      (flag.in.sumofmeas && match.nocase("^\\s*$", fields[1]))) {
    flag.in.sumofmeas = TRUE
    sumofmeas.values = c(sumofmeas.values, as.character(fields[-1]))
  }
}

##-- assign names to measurements list elements equal to their tags
names(measurements) = unlist(lapply(measurements, function(x) x@tag))

##-- get quantities measured by each experiment
meas.quantities = unlist(lapply(measurements, function(x) names(x@value)))

##-- build list of all measurements mentioned in the COMBINE section
meas.list = rep(FALSE, length(measurements))
names(meas.list) = names(measurements)
for (quant in combination@quantities) {
  meas.list = meas.list | (quant == meas.quantities)
}

##-- include measurements that correspond to combination of quantities
##++ should probably check also that _all_ quantities are in the combination
for (meas in names(combination@meas.lin.combs)) {
  if (sum(combination@quantities %in% names(combination@meas.lin.combs[[meas]])) != 0) {
    cat("meas", meas, "included, as linear combination\n")
    meas.list[meas] = TRUE
  }
}

##-- retain only measurements that are related to the spec. quantities
if (sum(!meas.list) > 0) {
  cat("warning: the following measurements are discarded\n")
  cat(paste("  ",names(measurements)[!meas.list], collapse="\n",sep=""), "\n")
  cat("end warning\n")
}
##-- update
measurements = measurements[meas.list]
meas.quantities = meas.quantities[meas.list]
meas.num = length(measurements[meas.list])
meas.names = names(measurements)

##-- check that the sum of syst. terms does not exceed the syst. error
for (meas in names(measurements)) {
  syst.contribs = sqrt(sum(measurements[[meas]]@syst.terms^2))
  if (syst.contribs > (1+1e-3)*measurements[[meas]]@syst) {
    stop("error: sum of syst. contributions larger than syst. error ",
         syst.contribs, ", ", measurements[[meas]]@syst)
  } else if (syst.contribs > (1+1e-5)*measurements[[meas]]@syst) {
    cat("warning: sum of syst. terms slightly larger than syst. error\n  ",
        syst.contribs, " vs. ", measurements[[meas]]@syst, "\n", sep="")
  }
}

##
## shift measurements according to updated external parameter dependencies
## update systematic terms according to updated external parameter errors
##
for (param.upd in names(combination@params)) {
  for (meas in names(measurements)) {
    for (param.orig in names(measurements[[meas]]@params)) {
      if (param.orig == param.upd) {
        measurements[[meas]]@value =
          (measurements[[meas]]@value
           + (combination@params[[param.upd]]["value"] - measurements[[meas]]@params[[param.orig]]["value"])
           * measurements[[meas]]@syst.terms[param.orig] / measurements[[meas]]@params[[param.orig]]["delta_pos"])
        measurements[[meas]]@syst.terms[param.orig] =
          (measurements[[meas]]@syst.terms[param.orig]
           * combination@params[[param.upd]]["delta_pos"] / measurements[[meas]]@params[[param.orig]]["delta_pos"])
      }
    }
  }
}

##-- collect what measurements are affected by each syst. term
syst.terms.list = list()
for (meas in names(measurements)) {
  for (syst.term in names(measurements[[meas]]@syst.terms)) {
    ##-- add measurement to the list of the currect syst. term
    syst.terms.list[[syst.term]] = c(syst.terms.list[[syst.term]], meas)
  }
}
##-- retain just the syst. contributions that affect at least two measurements
syst.terms.corr = lapply(syst.terms.list, length) >= 2
if (length(syst.terms.corr) > 0) {
  syst.terms.corr = names(syst.terms.corr)[syst.terms.corr]
} else {
  syst.terms.corr = character(0)
}

##
## add correlated syst. terms to the measurements vector
## as fake measurements equal to zero with error one
##
meas.quantities.true = meas.quantities
meas.num.true = meas.num
meas.names.true = meas.names
meas.num.fake = length(syst.terms.corr)
meas.num = meas.num.true + meas.num.fake
meas.names.fake = character(0)
if (meas.num.fake > 0) {
  meas.names.fake = paste(syst.terms.corr, "m", sep=".")
}
meas.names = c(meas.names.true, meas.names.fake)

##
## in the following, the covariance matrix for measurements is assembled
##

##-- get list of stat and syst errors
meas.stat.true = unlist(lapply(measurements, function(x) x@stat))
names(meas.stat.true) = meas.names.true
meas.syst.true = unlist(lapply(measurements, function(x) x@syst))
names(meas.syst.true) = meas.names.true

##-- fake measurements have stat. error = 1, syst. error = 0
meas.stat.fake = rep(1, meas.num.fake)
names(meas.stat.fake) = meas.names.fake
meas.syst.fake = rep(0, meas.num.fake)
names(meas.syst.fake) = meas.names.fake

##-- combine true and fake info, compute total error
meas.stat = c(meas.stat.true, meas.stat.fake)
meas.syst = c(meas.syst.true, meas.syst.fake)
meas.error = sqrt(meas.stat^2 + meas.syst^2)

##
## build correlation matrix
##
meas.corr = diag(rep(0,meas.num))
rownames(meas.corr) = meas.names
colnames(meas.corr) = meas.names
meas.corr.stat = meas.corr

##-- set off-diagonal statistical correlation matrix coefficients from cards
for (meas in measurements) {
  mapply(function(other.tag, other.corr) {
    meas.corr.stat[meas@tag,other.tag] <<- other.corr
  }, names(meas@corr.terms), meas@corr.terms)
}

##
## set off-diagonal total correlation matrix coefficients from cards
## total correlation terms are to be multiplied by the total errors
##
for (meas in measurements) {
  mapply(function(other.tag, other.corr) {
    meas.corr[meas@tag,other.tag] <<- other.corr
  }, names(meas@corr.terms.tot), meas@corr.terms.tot)
}

##-- not handled and forbidden to enter both total and stat. only correlations
flag.ok = TRUE
for (i in 1:meas.num) {
  for (j in i:meas.num) {
    if (meas.corr[i,j] != 0 && meas.corr.stat[i,j] != 0) {
      flag.ok = FALSE
      cat(paste("error: both total and statistical correlation specified for measurements:\n  ",
                meas.names[i], ", ", meas.names[j], "\n", collapse=""))
    }
  }
}
if (!flag.ok) {
  stop("aborted because of above errors\n")
}

##
## build covariance matrix using errors and correlation coefficients
## diagonal terms are the errors squared
## stat. correlation is multiplied by stat. errors
## total correlation is multiplied by total errors
##
meas.cov = meas.corr.stat * (meas.stat %o% meas.stat)
##
## keep note of total correlation coefficients that might have to be subtracted
## to remove the fraction of total correlation actually due to correlated systematics
##
meas.cov.tot = meas.corr * (meas.error %o% meas.error)
meas.cov = meas.cov + meas.cov.tot
meas.cov = meas.cov + diag(meas.error^2)

##
## from variance and covariance terms between true measurements
## we subtract the correlated systematic terms contributions
##
## for off-diagonal terms the subtraction is only done if
## total correlation terms were specified
##
meas.cov.orig = meas.cov
for (meas.i in meas.names.true) {
  syst.i = measurements[[meas.i]]@syst.terms
  for (meas.j in meas.names.true) {
    if (meas.i != meas.j && meas.cov.tot[meas.i,meas.j] == 0) next
    syst.j = measurements[[meas.j]]@syst.terms
    correl.i.j = intersect(names(syst.i), names(syst.j))
    correl.i.j = intersect(correl.i.j, syst.terms.corr)
    if (length(correl.i.j) == 0) next
    cov.contrib = sum(syst.i[correl.i.j] * syst.j[correl.i.j])
    meas.cov[meas.i,meas.j] = meas.cov[meas.i,meas.j] - cov.contrib
  }
}

##
## quantities we want to determine from measurements
##
quant.names.true = combination@quantities
quant.num.true = length(quant.names.true)
quant.num.fake = length(syst.terms.corr)
quant.num = quant.num.true + quant.num.fake
quant.names.fake = character(0)
if (quant.num.fake > 0) {
  quant.names.fake = paste(syst.terms.corr, "q", sep=".")
}
quant.names = c(quant.names.true, quant.names.fake)
##-- set quantity - measurement correspondence for fake ones
meas.quantities.true = meas.quantities
meas.quantities = c(meas.quantities, quant.names.fake)

delta = matrix(0, meas.num, quant.num)
colnames(delta) = quant.names
rownames(delta) = meas.names

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

##
## build column by column
## a 1 is set for each row where a measurement measures a quantity
##
for (quant in quant.names) {
  delta[,quant] = as.numeric(quant == meas.quantities)
}

##
## for measurements that are linear combination of quantities
## set the delta matrix coefficients as specified
##
for (meas in names(combination@meas.lin.combs)) {
  quants = names(combination@meas.lin.combs[[meas]])
  delta[meas,quants] = combination@meas.lin.combs[[meas]]
}

##
## for each measurement, systematic term = fake quantity,
## the delta coefficient is the respective syst. term
##
for (syst.term.name in names(syst.terms.list[syst.terms.corr])) {
  quant.name.fake = paste(syst.term.name,"q",sep=".")
  for (meas in syst.terms.list[syst.terms.corr][[syst.term.name]]) {
    delta[meas, quant.name.fake] = measurements[[meas]]@syst.terms[syst.term.name]
  }
}

##-- get measurement values
meas.true = unlist(lapply(measurements, function(x) {x@value}))
meas.fake = rep(0, meas.num.fake)
names(meas.fake) = meas.names.fake
meas = c(meas.true, meas.fake)

##-- print corrected measurements
if (FALSE) {
  show(meas.true)
  show(meas.stat.true)
  show(meas.syst.true)
  show(meas.error.true)
}

##
## analytical minimum chi-square solution for quantities
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
cat("## exact solution, chisq/d.o.f. = ",chisq, "/", meas.num - quant.num,
    ", CL= ", 100*(1-pchisq(chisq, df=meas.num-quant.num)), "%\n",sep="")
cat("##\n")
width.max = max(nchar(quant.names.true))
for (iquant in 1:quant.num.true) {
  cat(format(quant.names[iquant], width=width.max+2, justify = "right"), "=", quant[iquant], "+-", quant.err[iquant], "\n")
}
if (quant.num.true > 1) {
  cat("correlation\n")
  show(quant.corr[1:quant.num.true,1:quant.num.true])
}
cat("\n")

if (!flag.no.maxLik) {
##
## solve for quantities with iterative chi-square minimization
##

logLik.average = function(par) {
  chisq = t(meas - delta %*% par) %*% invcov %*% (meas - delta %*% par)
  return(-1/2*chisq)
}

fit = maxLik(logLik.average, start=quant*1.1, method="BFGS")

quant = coef(fit)
quant.cov = vcov(fit)
quant.err = sqrt(diag(quant.cov))
quant.corr = quant.cov / (quant.err %o% quant.err)

chisq = t(meas - delta %*% quant) %*% invcov %*% (meas - delta %*% quant)
chisq.fit = -2*logLik(fit)

cat("##\n")
cat("## numerical fit, chisq/d.o.f = ", chisq.fit, "/", meas.num - quant.num,
    ", CL= ", 100*(1-pchisq(chisq.fit, df=meas.num-quant.num)), "%\n",sep="")
cat("##\n")
width.max = max(nchar(quant.names.true))
for (iquant in 1:quant.num.true) {
  cat(format(quant.names[iquant], width=width.max+2, justify = "right"), "=", quant[iquant], "+-", quant.err[iquant], "\n")
}
if (quant.num.true > 1) {
  cat("correlation\n")
  show(quant.corr[1:quant.num.true,1:quant.num.true])
}
cat("\n")

cat("## begin fit summary\n")
show(fit)
cat("## end fit summary\n")
}

##++} ##-- end function alucomb

args <- commandArgs(TRUE)
if (length(args) > 0) alucomb(file = args[1]) 
