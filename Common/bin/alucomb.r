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
rc = try(library(maxLik), silent=TRUE)
if (inherits(rc, "try-error")) flag.no.maxLik = TRUE

source("../../../Common/bin/alucomb-read.r")

## ////////////////////////////////////////////////////////////////////////////
## definitions

## ////////////////////////////////////////////////////////////////////////////
## code

##
## alucomb
##

file = "average.input"

##++ alucomb = function(file = "") {

rc = alucomb.read(file)
measurements = rc$measurements
combination = rc$combination
rc = NULL

##
## build list of all measurements mentioned in the COMBINE section
## (right now we only handle COMBINE * * *)
##

##-- all measurements in the input cards
meas.names = names(measurements)
##-- quantities to be averaged
quant.names = combination$quantities
quant.num = length(quant.names)
##-- measurements of quantities that are linear combination of the quantities to be averaged
meas.names.comb = names(combination$meas.lin.combs[unlist(
  lapply(combination$meas.lin.combs, function(el) all(names(el) %in% combination$quantities)))])
##-- quantity measured per measurement
meas.quantities = unlist(lapply(measurements, function(x) names(x$value)))
names(meas.quantities) = meas.names
##-- unique quantities corresponding to linear combinations, other than averaged quantities
quant.names.comb = unique(meas.quantities[meas.names.comb])
quant.names.comb = setdiff(quant.names.comb, quant.names)
##-- quantities to be averaged plus their linear combinations
quant.names.include = c(quant.names, quant.names.comb)

if (length(quant.names.comb) >0) {
  cat("The following measurements are included as linear combinations:\n")
  cat(paste("  ", meas.names[meas.quantities %in% quant.names.comb], collapse="\n"), "\n");
}
meas.included.list = meas.quantities %in% quant.names.include
names(meas.included.list) = names(meas.quantities)[meas.included.list]
meas.names.discarded =  meas.names[!meas.included.list]
if (length(meas.names.discarded) >0) {
  cat("The following measurements are discarded:\n")
  cat(paste("  ", meas.names.discarded, collapse="\n"), "\n");
}

##-- update
measurements = measurements[meas.included.list]
meas.names = names(measurements)
meas.num = length(measurements[meas.included.list])
meas.quantities = meas.quantities[meas.included.list]

##-- checks on syst terms
for (meas in names(measurements)) {
  ##-- check that the sum of syst. terms does not exceed the syst. error
  syst.contribs = sqrt(sum(measurements[[meas]]$syst.terms^2))
  if (syst.contribs > (1+1e-3)*measurements[[meas]]$syst) {
    stop("error: syst. terms larger than syst. error\n    ",
         syst.contribs, " vs. ", measurements[[meas]]$syst, "\n    ",
         "measurement ", meas)
  } else if (syst.contribs > (1+1e-5)*measurements[[meas]]$syst) {
    cat("warning: syst. terms slightly larger than syst. error\n    ",
        syst.contribs, " vs. ", measurements[[meas]]$syst, "\n    ",
        "measurement ", meas, "\n", sep="")
  }
  ##-- warning if sum of syst terms different w.r.t. total syst
  if (abs(syst.contribs-measurements[[meas]]$syst) >
      1e-3 * sqrt(syst.contribs * measurements[[meas]]$syst)) {
    cat("warning: syst. terms do not match total syst. error, Combos requires that:\n  ",
            meas, ": ", measurements[[meas]]$syst, " vs. ", syst.contribs, "\n", sep="")
  }
}

##
## shift measurements according to updated external parameter dependencies
## update systematic terms according to updated external parameter errors
##
for (meas in names(measurements)) {
  value.delta = numeric()
  syst.term.deltasq = numeric()
  for (param.upd in names(combination$params)) {
    for (param.orig in names(measurements[[meas]]$params)) {
      if (param.orig == param.upd) {
        ##-- collect differences due to updated parameters
        value.delta = c(value.delta,
          ((combination$params[[param.upd]]["value"] - measurements[[meas]]$params[[param.orig]]["value"])
           * measurements[[meas]]$syst.terms[param.orig] / measurements[[meas]]$params[[param.orig]]["delta_pos"]))
        ##-- update systematic contribution according to updated parameter
        syst.term.orig = measurements[[meas]]$syst.terms[param.orig]
        syst.term.upd = syst.term.orig *
          combination$params[[param.upd]]["delta_pos"] / measurements[[meas]]$params[[param.orig]]["delta_pos"]
        measurements[[meas]]$syst.terms[param.orig] = syst.term.upd
        ##-- collect difference of syst term squares, to adjust the total systematic error as well
        syst.term.deltasq = c(syst.term.deltasq, (syst.term.upd^2 - syst.term.orig^2))
        if (FALSE) {
        cat(format(measurements[[meas]]$tag,width=30),
            format(param.orig,width=15),
            format(c(
            measurements[[meas]]$params[[param.orig]]["value"],
            combination$params[[param.upd]]["value"],
            (combination$params[[param.upd]]["value"] - measurements[[meas]]$params[[param.orig]]["value"])
            / measurements[[meas]]$params[[param.orig]]["delta_pos"],
            ## (value.upd - value.orig),
            measurements[[meas]]$params[[param.orig]]["delta_pos"],
            combination$params[[param.upd]]["delta_pos"],
            combination$params[[param.upd]]["delta_pos"] / measurements[[meas]]$params[[param.orig]]["delta_pos"]),
            width=10,digits=4,scientific=TRUE),
            "\n")
        }
      }
    }
  }
  ##-- update value
  measurements[[meas]]$value = measurements[[meas]]$value + sum(value.delta)
  ##-- update systematic error
  measurements[[meas]]$syst = sqrt(measurements[[meas]]$syst^2 + sum(syst.term.deltasq))
  ## cat(measurements[[meas]]$tag, measurements[[meas]]$value, measurements[[meas]]$stat, measurements[[meas]]$syst, "\n")
}

##-- for each syst. term external parameter, collect affected measurements
syst.terms.list = list()
for (meas in meas.names) {
  for (syst.term in names(measurements[[meas]]$syst.terms)) {
    ##-- add measurement to the list of the currect syst. term
    syst.terms.list[[syst.term]] = c(syst.terms.list[[syst.term]], meas)
  }
}

##-- retain just the syst. contributions that affect at least two measurements
syst.terms.corr = lapply(syst.terms.list, length) >= 2
if (length(syst.terms.corr) > 0) {
  syst.terms.corr = names(syst.terms.corr)[syst.terms.corr]
} else {
  syst.terms.corr = character()
}

##
## in the following, the covariance matrix for measurements is assembled
##

##-- get list of values, stat errors, syst errors
meas.value = unlist(lapply(measurements, function(x) x$value))
names(meas.value) = meas.names
meas.stat = unlist(lapply(measurements, function(x) x$stat))
names(meas.stat) = meas.names
meas.syst = unlist(lapply(measurements, function(x) x$syst))
names(meas.syst) = meas.names
meas.error = sqrt(meas.stat^2 + meas.syst^2)

##
## build correlation matrix
##
meas.corr = diag(rep(0,meas.num))
rownames(meas.corr) = meas.names
colnames(meas.corr) = meas.names
meas.corr.stat = meas.corr

##
## set off-diagonal correlation matrix coefficients from cards
## - meas.corr.stat means only stat. correlation, to be multiplied by stat. errors
## - meas.corr means total correlation, to be multiplied by total errors
##
##-- set off-diagonal statistical correlation matrix coefficients from cards
for (mi.name in meas.names) {
  for (mj.name in intersect(names(measurements[[mi.name]]$corr.terms), meas.names)) {
    meas.corr.stat[meas.names %in% mi.name, meas.names %in% mj.name] = measurements[[mi.name]]$corr.terms[[mj.name]]
  }
  for (mj.name in intersect(names(measurements[[mi.name]]$corr.terms.tot), meas.names)) {
    meas.corr[meas.names %in% mi.name, meas.names %in% mj.name] = measurements[[mi.name]]$corr.terms.tot[[mj.name]]
  }
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
## - stat. correlation is multiplied by stat. errors
## - total correlation is multiplied by total errors
##
meas.cov = meas.corr * (meas.error %o% meas.error)
meas.cov.stat = meas.corr.stat * (meas.stat %o% meas.stat)

##
## get syst. correlation corresponding to correlated syst. terms
##
meas.cov.syst = meas.corr * 0
for (meas.i in meas.names) {
  syst.i = measurements[[meas.i]]$syst.terms
  for (meas.j in meas.names) {
    ##-- no addition needed for on-diagonal terms
    if (meas.i == meas.j) next
    syst.j = measurements[[meas.j]]$syst.terms
    ##-- systematics common to the two measurements
    correl.i.j = intersect(names(syst.i), names(syst.j))
    ##-- remove syst. terms uncorrelated to two different measurements
    correl.i.j = intersect(correl.i.j, syst.terms.corr)
    if (length(correl.i.j) == 0) next
    meas.cov.syst[meas.i,meas.j] = sum(syst.i[correl.i.j] * syst.j[correl.i.j])
  }
}

##-- if total correlation specified, get stat. correlation by subtraction
meas.cov.stat = ifelse(meas.cov == 0, meas.cov.stat, meas.cov - meas.cov.syst)
meas.cov.stat = meas.cov.stat + diag(meas.stat^2)
##-- total covariance
meas.cov.syst = meas.cov.syst + diag(meas.syst^2)
meas.cov = meas.cov.stat + meas.cov.syst

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
delta = matrix(0, meas.num, quant.num)
colnames(delta) = quant.names
rownames(delta) = meas.names

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
for (meas in names(combination$meas.lin.combs)) {
  quants = names(combination$meas.lin.combs[[meas]])
  delta[meas,quants] = combination$meas.lin.combs[[meas]]
}

##-- print corrected measurements
if (FALSE) {
  show(meas.value)
  show(meas.stat)
  show(meas.syst)
  show(meas.error)
}

##-- simplify
meas = meas.value

if (FALSE && !flag.no.maxLik) {
##
## solve for quantities with iterative chi-square minimization
##

##-- compute weight matrix for computing chisq
invcov = solve(meas.cov)

logLik.average = function(par) {
  chisq = t(meas - delta %*% par) %*% invcov %*% (meas - delta %*% par)
  return(-1/2*chisq)
}

fit = maxLik(logLik.average, start=rep(0,quant.num), method="BFGS")

quant = coef(fit)
names(quant) = quant.names
quant.cov = vcov(fit)
rownames(quant.cov) = quant.names
colnames(quant.cov) = quant.names
quant.err = sqrt(diag(quant.cov))
quant.corr = quant.cov / (quant.err %o% quant.err)

chisq = drop(t(meas - delta %*% quant) %*% invcov %*% (meas - delta %*% quant))
chisq.fit = -2*logLik(fit)

cat("\n")
cat("## begin fit summary\n")
show(fit)
cat("## end fit summary\n")

cat("\n")
cat("##\n")
cat("## numerical fit, chisq/d.o.f = ", chisq.fit, "/", meas.num - quant.num,
    ", CL = ", (1-pchisq(chisq.fit, df=meas.num-quant.num)), "\n",sep="")
cat("##\n")
show(rbind(value=quant[1:quant.num], error=quant.err[1:quant.num]))
if (quant.num > 1) {
  cat("correlation\n")
  show(quant.corr[1:quant.num, 1:quant.num])
}
} ## !flag.no.maxLik && FALSE

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

chisq = drop(t(meas - delta %*% quant) %*% invcov %*% (meas - delta %*% quant))

cat("\n")
cat("##\n")
cat("## exact solution, chisq/d.o.f. = ",chisq, "/", meas.num - quant.num,
    ", CL = ", (1-pchisq(chisq, df=meas.num-quant.num)), "\n",sep="")
cat("##\n")
show(rbind(value=quant[1:quant.num], error=quant.err[1:quant.num]))
if (quant.num > 1) {
  cat("correlation\n")
  show(quant.corr[1:quant.num, 1:quant.num])
}

##
## each measurement is a linear combination of the quantities we fit
## here we collect all the unique linear combinations, named "types"
##
meas.types.id = unique(delta)
##++ get measurement names from "method" rather than MEASUREMENT card for safety
rownames(meas.types.id) = sub("[^.]*.([^.]*).[^.]*", "\\1", rownames(meas.types.id), perl=TRUE)

##-- for each "type", will set TRUE at the position of corresponding measurements
meas.types = matrix(FALSE, dim(meas.types.id)[1], meas.num)
meas.types.names = rownames(meas.types.id)
rownames(meas.types) = meas.types.names
colnames(meas.types) = meas.names
for (mt.name in meas.types.names) {
  for (m.name in meas.names) {
    meas.types[mt.name, m.name] = all(meas.types.id[mt.name,] == delta[m.name,])
  }
}

##
## measurement types that do not correspond to an averaged quantity
## - compute average and errors
##
meas.names.extra = meas.types.names[!(meas.types.names %in% quant.names)]
meas.extra.id = subset(meas.types.id, meas.types.names %in% meas.names.extra)
meas.extra = drop(meas.extra.id %*% quant[1:quant.num])
meas.extra.err = sqrt(diag(meas.extra.id %*% quant.cov %*% t(meas.extra.id)))

##-- chisq contribution and dof for each measurement type
chisq.types = rep(0, dim(meas.types.id)[1])
names(chisq.types) = meas.types.names
dof.types = chisq.types

##-- S-factor for each measurement type
sfact.types = rep(1, dim(meas.types.id)[1])
names(sfact.types) = meas.types.names

##-- S-factor for each true measurement
sfact = rep(1, meas.num)
names(sfact) = meas.names
keep.types.names = character(0)

##
## collect chisq contributions for each measurement type
## the following is filled:
## - matrix meas.types[meas.types,quant) with TRUE where a meas.type has a quantity as addendum
## - chisq.types: chisq per measurement type
## - dof.types: dof per measurement type
## - sfact.types: S-factor per measurement type
## - sfact: S-factor per measurement
##
for (mt.name in meas.types.names) {
  chisq.types.meas = numeric()
  ##-- linear comb. of averaged quantities corresponding to mt.name
  quant.comb = meas.types.id[mt.name,]
  mt.m.names = meas.names[meas.types[mt.name,]]
  for (m.name in mt.m.names) {
    ##-- chisq contribution of the single measurement
    tmp = (meas[m.name] - drop(quant.comb %*% quant)) / meas.error[m.name]
    chisq.types.meas = c(chisq.types.meas, tmp^2)
  }

  ##-- do not include chisq contribution of measurements with error >sqrt(3N) times average error (PDG recipe)
  num = length(chisq.types.meas)
  ##-- error on average as resulting from HFAG fit
  ## average.err = c(quant.err, meas.extra.err)[mt.name]
  ##-- error on average like for PDG, assuming no correlation
  average.err = 1/sqrt(sum(1/meas.error[mt.m.names]^2))
  error.max = 3*sqrt(num)*average.err
  keep = ifelse(meas.error[mt.m.names] <= error.max, 1, 0)
  ##++ keep = ifelse(meas.error[mt.m.names] <= (3*sqrt(num)*average.err), 1, 1)
  keep.types.names = c(keep.types.names, names(keep)[keep != 0])
  save = options()
  options(digits=4)
  if (any(keep == 0)) {
    cat("\nS-factor calculation, exclude because error > 3*sqrt(N)*av_err=", error.max, "\n")
    excl = meas.error[names(keep)[keep == 0]]
    cat(paste(names(excl), "error=", format(excl, digits=4), "nsigma=", format(excl/error.max*3, digits=4), sep=" "), sep="\n")
  }
  if (FALSE && any(keep != 0)) {
    cat("\nS-factor calculation, included measurements, 3*sqrt(N)*av_err=", error.max, "\n")
    excl = meas.error[names(keep)[keep != 0]]
    cat(paste(names(excl), "error=", format(excl, digits=4), "nsigma=", format(excl/error.max*3, digits=4), sep=" "), sep="\n")
  }
  options(save)
  ##-- chisq/dof for each type of measurement
  chisq.types[mt.name] = sum(keep*chisq.types.meas)
  ##-- number of degrees of freedom associated with measurements type
  dof.types[mt.name] = sum(keep != 0) - 1
  ##++ special treatment for when there is just 1 measurement of a specific type
  if (dof.types[mt.name] == 0) dof.types[mt.name] = 1
  ## cat(mt.name, "chisq=", chisq.types[mt.name], "dof= ", dof.types[mt.name], "\n")

  if (FALSE) {
    ##
    ## experimental code to compute chisq pertaining to meas. type by
    ## subsetting the chisq matrix formula
    ## - either subset to measurements of spec.type
    ## - or subset to all other meas. types and get difference from total chisq
    ##
    ## mt.select = !meas.types[mt.name,]
    mt.select = !meas.types[mt.name,]
    meas.cov.redu = subset(meas.cov, subset=mt.select, select=mt.select)
    invcov.redu = solve(meas.cov.redu)
    meas.delta = subset(meas - delta %*% quant, subset=mt.select)
    
    chisq.partial = drop(t(meas.delta) %*% invcov.redu %*% meas.delta)
    chisq.partial = chisq - chisq.partial
    ## dof.partial = length(meas.delta) -1
    dof.partial = (meas.num - quant.num) - (length(meas.delta) - (quant.num-1))
    
    chisq.types[mt.name] = chisq.partial
    dof.types[mt.name] = dof.partial
    ##++ special treatment for when there is just 1 measurement of a specific type
    if (dof.types[mt.name] == 0) dof.types[mt.name] = 1
    ##-- end
  }
  
  ##-- S-factor for each type of measurement and for each measurement
  tmp = sqrt(chisq.types[mt.name]/dof.types[mt.name])
  sfact.types[mt.name] = tmp
  sfact[meas.types[mt.name,]] = tmp
}

##-- recompute chisq and chisq/dof
dof = meas.num - quant.num
chisq = drop(t(meas - delta %*% quant) %*% invcov %*% (meas - delta %*% quant))

##-- update chisq using S-factors
meas.keep = ifelse(meas.names %in% keep.types.names, 1, 0)
names(meas.keep) = meas.names
dof.keep = length(keep.types.names) - quant.num
chisq.keep = drop(
  t(meas - delta %*% quant) %*% diag(meas.keep)
  %*% solve(diag(sfact) %*% meas.cov %*% diag(sfact))
  %*% diag(meas.keep) %*% (meas - delta %*% quant))

##-- adjust found S-factors to obtain chisq/dof = 1
sfact = sfact * sqrt(chisq.keep/dof.keep)
sfact.types = sfact.types * sqrt(chisq.keep/dof.keep)

sfact.types.orig = sfact.types
sfact.types = pmax(sfact.types, 1)
sfact.orig = sfact
sfact = pmax(sfact, 1)

##
## inflate the covariance matrix with the S-factors
##
meas2.cov = diag(sfact) %*% meas.cov %*% diag(sfact)
invcov2 = solve(meas2.cov)

##-- compute new inflated fitted quantities covariance matrix 
quant2.cov = solve(t(delta) %*% invcov2 %*% delta)
rownames(quant2.cov) = quant.names
colnames(quant2.cov) = quant.names

##-- updated errors and correlations
quant2.err = sqrt(diag(quant2.cov))
quant2.corr = quant2.cov / (quant2.err %o% quant2.err)

##-- all measurements chisq after S-factor inflation
chisq2 = drop(t(meas - delta %*% quant) %*% invcov2 %*% (meas - delta %*% quant))

##-- no-large-error measurements chisq after 2nd iteration S-factor inflation
chisq2.keep = drop(
  t(meas - delta %*% quant) %*% diag(meas.keep)
  %*% invcov2
  %*% diag(meas.keep) %*% (meas - delta %*% quant))

##-- no-large-error measurements chisq before any S-factor inflation
chisq.keep.0 = drop(
  t(meas - delta %*% quant) %*% diag(meas.keep)
  %*% invcov
  %*% diag(meas.keep) %*% (meas - delta %*% quant))

cat("\n")
cat("##\n")
cat("## S-factors accounting for larger than expected chi-quare\n")
cat("##\n")

show(rbind(original = c(
             chisq=chisq, dof=dof,
             "chisq/dof"=chisq/dof,
             CL=pchisq(chisq, dof, lower.tail=FALSE))
           ,"  after S-factor, 2nd stage" = c(
              chisq=chisq2, dof=dof,
              "chisq/dof"=chisq2/dof,
              CL=pchisq(chisq2, dof, lower.tail=FALSE))
           ,"no-large-error" = c(
                        chisq=chisq.keep.0, dof=dof.keep,
                        "chisq/dof"=chisq.keep.0/dof.keep,
                        CL=pchisq(chisq.keep.0, dof.keep, lower.tail=FALSE))
           ,"  after S-factor, 1st stage" = c(
                        chisq=chisq.keep, dof=dof.keep,
                        "chisq/dof"=chisq.keep/dof.keep,
                        CL=pchisq(chisq.keep, dof.keep, lower.tail=FALSE))
           ,"  after S-factor, 2nd stage" = c(
                        chisq=chisq2.keep, dof=dof.keep,
                        "chisq/dof"=chisq2.keep/dof.keep,
                        CL=pchisq(chisq2.keep, dof.keep, lower.tail=FALSE))
           ))

cat("Averaged quantities: value, error, error with S-factor, S-factor\n") 
if (quant.num > 1) {
  ##-- if multiple average, use dedicated S-factor computation
  sfact.row = quant2.err / quant.err
  sfact0.row = sfact.types[meas.types.names %in% quant.names]
  chisq.row = chisq.types[meas.types.names %in% quant.names]
  dof.row = dof.types[meas.types.names %in% quant.names]
} else {
  ##-- if averaging a single quantity, use global chisq to compute S-factor
  sfact0.row = sfact.types[meas.types.names %in% quant.names]
  sfact.row = sqrt(chisq/dof)
  quant2.err[1] = quant.err[1]*sfact.row[1]
  chisq.row = chisq
  dof.row = dof
}
show(rbind(value=quant,
           error=quant.err,
           upd.error=quant2.err,
           "S-factor"=sfact.row,
           "S-factor_0"=sfact0.row,
           chisq=chisq.row,
           dof=dof.row
           ))

if (quant.num > 1) {
  cat("correlation\n") 
  show(quant2.corr)
}

##-- update with S-factors measurement types that do not correspond to an averaged quantity
meas.extra.err.upd = sqrt(diag(meas.extra.id %*% quant2.cov %*% t(meas.extra.id)))

if (length(meas.extra) >0) {
  cat("Non-averaged measurement types: value, error, error with S-factor, S-factor\n") 
  show(rbind(value=meas.extra,
             error=meas.extra.err,
             upd.error=meas.extra.err.upd,
             "S-factor"=meas.extra.err.upd/meas.extra.err,
             "S-factor_0"=sfact.types[meas.types.names %in% meas.names.extra],
             chisq=chisq.types[meas.types.names %in% meas.names.extra],
             dof=dof.types[meas.types.names %in% meas.names.extra]
           ))
}

##++ } ##-- end function alucomb

args <- commandArgs(TRUE)
if (length(args) > 0 && exists("alucomb")) alucomb(file = args[1]) 
