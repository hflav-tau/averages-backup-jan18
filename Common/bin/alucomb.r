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

##-- get quantities measured by each experiment
meas.quantities = unlist(lapply(measurements, function(x) names(x$value)))

##-- build list of all measurements mentioned in the COMBINE section
meas.list = rep(FALSE, length(measurements))
names(meas.list) = names(measurements)
for (quant in combination$quantities) {
  meas.list = meas.list | (quant == meas.quantities)
}

##-- include measurements that correspond to combination of quantities
##++ should probably check also that _all_ quantities are in the combination
for (meas in names(combination$meas.lin.combs)) {
  if (sum(combination$quantities %in% names(combination$meas.lin.combs[[meas]])) != 0) {
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

##-- collect what measurements are affected by each syst. term
syst.terms.list = list()
for (meas in names(measurements)) {
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
## add correlated syst. terms to the measurements vector
## as fake measurements equal to zero with error one
##
meas.quantities.true = meas.quantities
meas.num.true = meas.num
meas.names.true = meas.names
meas.num.fake = length(syst.terms.corr)
meas.num = meas.num.true + meas.num.fake
meas.names.fake = character()
if (meas.num.fake > 0) {
  meas.names.fake = paste(syst.terms.corr, "m", sep=".")
}
meas.names = c(meas.names.true, meas.names.fake)

##
## in the following, the covariance matrix for measurements is assembled
##

##-- get list of stat and syst errors
meas.stat.true = unlist(lapply(measurements, function(x) x$stat))
names(meas.stat.true) = meas.names.true
meas.syst.true = unlist(lapply(measurements, function(x) x$syst))
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
    meas.corr.stat[meas$tag,other.tag] <<- other.corr
  }, names(meas$corr.terms), meas$corr.terms)
}

##
## set off-diagonal total correlation matrix coefficients from cards
## total correlation terms are to be multiplied by the total errors
##
for (meas in measurements) {
  mapply(function(other.tag, other.corr) {
    meas.corr[meas$tag,other.tag] <<- other.corr
  }, names(meas$corr.terms.tot), meas$corr.terms.tot)
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
  syst.i = measurements[[meas.i]]$syst.terms
  for (meas.j in meas.names.true) {
    if (meas.i != meas.j && meas.cov.tot[meas.i,meas.j] == 0) next
    syst.j = measurements[[meas.j]]$syst.terms
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
quant.names.true = combination$quantities
quant.num.true = length(quant.names.true)
quant.num.fake = length(syst.terms.corr)
quant.num = quant.num.true + quant.num.fake
quant.names.fake = character()
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
for (meas in names(combination$meas.lin.combs)) {
  quants = names(combination$meas.lin.combs[[meas]])
  delta[meas,quants] = combination$meas.lin.combs[[meas]]
}

##
## for each measurement, systematic term = fake quantity,
## the delta coefficient is the respective syst. term
##
for (syst.term.name in names(syst.terms.list[syst.terms.corr])) {
  quant.name.fake = paste(syst.term.name,"q",sep=".")
  for (meas in syst.terms.list[syst.terms.corr][[syst.term.name]]) {
    delta[meas, quant.name.fake] = measurements[[meas]]$syst.terms[syst.term.name]
  }
}

##-- get measurement values
meas.true = unlist(lapply(measurements, function(x) { x$value }))
names(meas.true) = meas.names.true
meas.fake = rep(0, meas.num.fake)
names(meas.fake) = meas.names.fake
names(meas.fake) = meas.names.fake
meas = c(meas.true, meas.fake)

##-- print corrected measurements
if (FALSE) {
  show(meas.true)
  show(meas.stat.true)
  show(meas.syst.true)
  show(meas.error.true)
}

if (!flag.no.maxLik) {
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
show(rbind(value=quant[1:quant.num.true], error=quant.err[1:quant.num.true]))
if (quant.num.true > 1) {
  cat("correlation\n")
  show(quant.corr[1:quant.num.true,1:quant.num.true])
}
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

chisq = drop(t(meas - delta %*% quant) %*% invcov %*% (meas - delta %*% quant))

cat("\n")
cat("##\n")
cat("## exact solution, chisq/d.o.f. = ",chisq, "/", meas.num - quant.num,
    ", CL = ", (1-pchisq(chisq, df=meas.num-quant.num)), "\n",sep="")
cat("##\n")
show(rbind(value=quant[1:quant.num.true], error=quant.err[1:quant.num.true]))
if (quant.num.true > 1) {
  cat("correlation\n")
  show(quant.corr[1:quant.num.true,1:quant.num.true])
}

##
## each measurement is a linear combination of the quantities we fit
## here we collect all the unique linear combinations, named "types"
##
meas.types.id = unique(delta[1:meas.num.true,1:quant.num.true,drop=FALSE])
rownames(meas.types.id) = sub("[^.]*.([^.]*).[^.]*", "\\1", rownames(meas.types.id), perl=TRUE)

##-- for each "type", will set TRUE at the position of corresponding measurements
meas.types = matrix(FALSE, dim(meas.types.id)[1], meas.num.true)
meas.types.names = rownames(meas.types.id)
rownames(meas.types) = meas.types.names
colnames(meas.types) = meas.names.true
         
##-- chisq contribution and dof for each measurement type
chisq.types = rep(0, dim(meas.types.id)[1])
names(chisq.types) = meas.types.names
dof.types = chisq.types

##-- S-factor for each measurement type
sfact.types = rep(1, dim(meas.types.id)[1])
names(sfact.types) = meas.types.names

##-- S-factor for each true measurement
sfact.true = rep(1, meas.num.true)
names(sfact.true) = meas.names.true

##
## collect chisq contributions for each measurement type
## the following is filled:
## - matrix meas.types[meas.types,quant) with TRUE where a meas.type has a quantity as addendum
## - chisq.types: chisq per measurement type
## - dof.types: dof per measurement type
## - sfact.types: S-factor per measurement type
## - sfact.true: S-factor per true measurement
##
for (mt.name in meas.types.names) {
  chisq.types.meas = numeric()
  for (m.name in meas.names.true) {
    quant.comb = quant.comb = (delta[1:meas.num.true,1:quant.num.true,drop=FALSE])[m.name,]
    if (all(quant.comb == meas.types.id[mt.name,])) {
      ##-- take note which measurements belong to each type
      meas.types[mt.name, m.name] = TRUE
      ##-- chisq contribution of the single measurement
      tmp = (meas[m.name] - drop(quant.comb %*% quant[1:quant.num.true])) / meas.error[m.name]
      chisq.types.meas = c(chisq.types.meas, tmp^2)
    }
  }
  ##-- chisq/dof for each type of measurement
  chisq.types[mt.name] = sum(chisq.types.meas)
  ##-- number of degrees of freedom associated with measurements type
  dof.types[mt.name] = length(chisq.types.meas)-1
  ##++ special treatment for when there is just 1 measurement of a specific type
  if (dof.types[mt.name] == 0) dof.types[mt.name] = 1

  if (FALSE) {
    ##
    ## experimental code to compute chisq pertaining to meas. type by
    ## subsetting the chisq matrix formula
    ## - either subset to measurements of spec.type
    ## - or subset to all other meas. types and get difference from total chisq
    ##
    ## mt.select = !meas.types[mt.name,]
    mt.select = c(!meas.types[mt.name,], rep(TRUE, meas.num.fake))
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
  sfact.true[meas.types[mt.name,]] = tmp
}

sfact.types.floored = pmax(sfact.types, 1)
sfact.true.floored = pmax(sfact.true, 1)

##
## S-factor per measurement
## - true measurements get the computed S-factors
## - fake measurements get S-factor=1
## if any S-factor is less than 1, it is set to 1
##
sfact = rep(1, meas.num)
names(sfact) = meas.names
sfact[1:meas.num.true] = sfact.true.floored

##
## inflate the true+fake (=whole) covariance matrix with the S-factors
##
meas2.cov = diag(sfact) %*% meas.cov %*% diag(sfact)
invcov2 = solve(meas2.cov)

##
## inflate the delta matrix coefficients corresponding to the correlated
## systematic errors with the proper per measurement S-factor
##
delta2 = delta
for (m.name in meas.names.true) {
  delta2[m.name, quant.names %in% quant.names.fake] = delta[m.name, quant.names %in% quant.names.fake] * sfact[m.name]
}

##-- compute new inflated fitted quantities covariance matrix 
quant2.cov = solve(t(delta2) %*% invcov2 %*% delta2)
rownames(quant2.cov) = quant.names
colnames(quant2.cov) = quant.names

##-- updated errors and correlations
quant2.err = sqrt(diag(quant2.cov))
quant2.corr = quant2.cov / (quant2.err %o% quant2.err)

##-- recompute chisq and chisq/dof
chisq = drop(t(meas - delta %*% quant) %*% invcov %*% (meas - delta %*% quant))

##
## compute new chi-square
## ++ just here we use re-minimized quantities with S-factor inflated errors
## ++ in this way the chisq ends up the same when replacing covariance with
## ++ dummy external parameters
##
quant.sfact = drop(quant2.cov %*% t(delta2) %*% (invcov2 %*% meas))
names(quant.sfact) = quant.names
chisq2 = drop(t(meas - delta2 %*% quant.sfact) %*% invcov2 %*% (meas - delta2 %*% quant.sfact))

cat("\n")
cat("##\n")
cat("## S-factors accounting for larger than expected chi-quare\n")
cat("##\n")
dof = (meas.num-quant.num)

##-- compute chisq for true measurements only (experimental)
meas.cov.true = meas.cov[1:meas.num.true,1:meas.num.true,drop=FALSE]
meas.delta = (meas - delta %*% quant)[1:meas.num.true, drop=FALSE]
invcov.true = solve(meas.cov.true)
chisq.true = drop(t(meas.delta) %*% invcov.true %*% meas.delta)
meas2.cov.true = meas2.cov[1:meas.num.true,1:meas.num.true,drop=FALSE]
meas2.delta = (meas - delta2 %*% quant)[1:meas.num.true, drop=FALSE]
invcov.true = solve(meas2.cov.true)
chisq2.true = drop(t(meas2.delta) %*% invcov.true %*% meas2.delta)

show(rbind(original   = c(chisq=chisq, dof=dof, "chisq/dof"=chisq/dof, CL=pchisq(chisq,dof,lower.tail=FALSE))
           ,updated    = c(chisq=chisq2, dof=dof, "chisq/dof"=chisq2/dof, CL=pchisq(chisq2,dof,lower.tail=FALSE))
           ##,original.t = c(chisq=chisq.true, dof=dof, "chisq/dof"=chisq.true/dof, CL=pchisq(chisq.true,dof,lower.tail=FALSE))
           ##,updated.t  = c(chisq=chisq2.true, dof=dof, "chisq/dof"=chisq2.true/dof, CL=pchisq(chisq2.true,dof,lower.tail=FALSE))
           ))

cat("Averaged quantities: value, error, error with S-factor, S-factor\n") 
if (quant.num.true > 1) {
  ##-- if multiple average, use dedicated S-factor computation
  sfact.row = quant2.err[1:quant.num.true] / quant.err[1:quant.num.true]
  sfact0.row = sfact.types[meas.types.names %in% quant.names.true]
  chisq.row = chisq.types[meas.types.names %in% quant.names.true]
  dof.row = dof.types[meas.types.names %in% quant.names.true]
} else {
  ##-- if averaging a single quantity, use global chisq to compute S-factor
  sfact0.row = sfact.types[meas.types.names %in% quant.names.true]
  sfact.row = sqrt(chisq/dof)
  quant2.err[1] = quant.err[1]*sfact.row[1]
  chisq.row = chisq
  dof.row = dof
}
show(rbind(value=quant[1:quant.num.true],
           error=quant.err[1:quant.num.true],
           upd.error=quant2.err[1:quant.num.true],
           "S-factor"=sfact.row,
           "S-factor_0"=sfact0.row,
           chisq=chisq.row,
           dof=dof.row
           ))

if (quant.num.true > 1) {
  cat("correlation\n") 
  show(quant2.corr[1:quant.num.true,1:quant.num.true])
}

##-- measurement types that do not correspond to an averaged quantity
meas.names.extra = meas.types.names[!(meas.types.names %in% quant.names)]
##meas.extra.id = meas.types.id[meas.types.names %in% meas.names.extra,]
meas.extra.id = subset(meas.types.id, meas.types.names %in% meas.names.extra)
meas.extra = drop(meas.extra.id %*% quant[1:quant.num.true])
meas.extra.err = sqrt(diag(meas.extra.id %*% quant.cov[1:quant.num.true,1:quant.num.true] %*% t(meas.extra.id)))
meas.extra.err.upd = sqrt(diag(meas.extra.id %*% quant2.cov[1:quant.num.true,1:quant.num.true] %*% t(meas.extra.id)))

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
