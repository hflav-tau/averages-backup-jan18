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

source("../../../Common/bin/alu-utils.r")

## ////////////////////////////////////////////////////////////////////////////
## definitions

## ////////////////////////////////////////////////////////////////////////////
## code

##
## alucomb
##

args <- commandArgs(TRUE)
if (length(args) > 0) {
  file = args[1]
} else {
  file = "average.input"
}

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

##
## transform COMBOFMEAS and SUMOF meas information
## (quantities that are combination of other quantities)
## into constraints consisting in combinaiton of values equal to constrats
##
tmp = mapply(function(name, value) {
  rc = c(value, -1)
  names(rc)[length(rc)] = name
  rc
},
  names(combination$meas.lin.combs), combination$meas.lin.combs,
  SIMPLIFY=FALSE)
if (length(tmp) > 0) {
  names(tmp) = paste(names(tmp), "l", sep=".")
  combination$constr.comb = c(combination$constr.comb, tmp)
  tmp2 = as.list(rep(0, length(tmp)))
  names(tmp2) = names(tmp)
  combination$constr.val = c(combination$constr.val, tmp2)
  rm(tmp2)
}
rm(tmp)

##-- retain only constraints whose terms are all included in the fitted quantities
constr.select = sapply(combination$constr.comb, function(x) all(names(x) %in% quant.names))
if (any(!constr.select)) {
  cat("The following constraints are dropped:\n")
  mapply(function(comb, val) {
    show(c(comb, constr=val))
  }, combination$constr.comb[!constr.select], combination$constr.val[!constr.select])
}
combination$constr.comb = combination$constr.comb[constr.select]
combination$constr.val = combination$constr.val[constr.select]

##-- quantity measured per measurement
meas.quantities = unlist(lapply(measurements, function(x) names(x$value)))
names(meas.quantities) = meas.names

meas.included.list = meas.quantities %in% quant.names
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
      if (param.orig == param.upd && param.orig %in% names(measurements[[meas]]$syst.terms)) {
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
meas.corr = diag.m(rep(0,meas.num))
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
meas.cov.stat = meas.cov.stat + diag.m(meas.stat^2)
##-- total covariance
meas.cov.syst = meas.cov.syst + diag.m(meas.syst^2)
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

##-- print corrected measurements
if (FALSE) {
  show(meas.value)
  show(meas.stat)
  show(meas.syst)
  show(meas.error)
}

##-- simplify
meas = meas.value

##-- set very large line width to print even large amount of averaged quantities on single line
options.save = options()
options(width=2000)

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
quant.err = sqrt(diag.m(quant.cov))
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
  ##++ show(quant.corr[1:quant.num, 1:quant.num])
}
} ## !flag.no.maxLik && FALSE

##
## analytical minimum chi-square solution for quantities
##

quant = rep(0, quant.num)
names(quant) = quant.names

meas.invcov = solve(meas.cov)

quant.invcov = t(delta) %*% meas.invcov %*% delta
##-- multiply constraint equation to avoid unbalanced singular values
quant.invcon.order = 10^round(log(det(quant.invcov)^(1/quant.num))/log(10))

##-- constraints
constr.num = length(quant.names.comb)
constr.names = names(combination$constr.comb)
constr.m =  do.call(rbind, lapply(combination$constr.comb, function(x) {tmp = quant; tmp[names(x)] = x; tmp}))
constr.m = quant.invcon.order * constr.m
constr.v = quant.invcon.order * unlist(combination$constr.val)

##-- build full matrix in front of c(quant vector, lagr.mult. vector)
full.m = rbind(
  cbind(quant.invcov, t(constr.m)),
  cbind(constr.m, matrix(0, constr.num, constr.num)))
##-- build full matrix in front of c(meas, constraint values)
full.v.m = rbind(
  cbind(t(delta) %*% meas.invcov, matrix(0, quant.num, constr.num, dimnames=list(NULL, constr.names))),
  cbind(matrix(0, constr.num, meas.num), diag.m(rep(1, constr.num))))
##-- build full vector c(measurements vector, constraint values)
full.v = c(meas, constr.v)

##-- matrix that applied to c(meas, constr. values) gives quant
solve.m = solve(full.m) %*% full.v.m

##-- solve for both quantities and lagrange multipliers
## quant.constr.val = solve(full.m, (full.v.m %*% full.v))
quant.constr.val = solve.m %*% full.v
quant = quant.constr.val[1:quant.num, drop=FALSE]
names(quant) = quant.names

##
## full covariance of measurements and constraint values
## - there is no error on the constraint values
##
full.v.cov = rbind(
  cbind(meas.cov, matrix(0, meas.num, constr.num, dimnames=list(NULL, constr.names))),
  cbind(matrix(0, constr.num, meas.num, dimnames=list(constr.names)), matrix(0, constr.num, constr.num)))

##-- covariance matrix of fitted values
quant.constr.cov = solve.m %*% full.v.cov %*% t(solve.m)
quant.cov = quant.constr.cov[1:quant.num, 1:quant.num, drop=FALSE]
## rownames(quant.cov) = quant.names
## colnames(quant.cov) = quant.names
quant.err = sqrt(diag.m(quant.cov))

quant = drop(quant.cov %*% t(delta) %*% (meas.invcov %*% meas))
names(quant) = quant.names

quant.corr = quant.cov / (quant.err %o% quant.err)

chisq = drop(t(meas - delta %*% quant) %*% meas.invcov %*% (meas - delta %*% quant))

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
cat("## end\n")

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

##-- chisq contribution, dof, S-factor for each measurement type
chisq.types = 0*meas.types.id[,1]
names(chisq.types) = meas.types.names
dof.types = chisq.types
num.types = chisq.types
num.types.keep = chisq.types
sfact.types = chisq.types + 1

##-- S-factor for each measurement
sfact = 0*meas + 1
meas.keep = meas & FALSE

##
## collect chisq contributions for each measurement type
## the following is filled:
## - matrix meas.types[meas.types,quant) with TRUE where a meas.type has a quantity as addendum
## - chisq.types: chisq per measurement type
## - dof.types: dof per measurement type
## - sfact.types: S-factor per measurement type
## - sfact: S-factor per measurement
##
chisq.out = numeric(0)
for (mt.name in meas.types.names) {
  ##-- linear comb. of averaged quantities corresponding to mt.name
  quant.comb = meas.types.id[mt.name,]

  ##-- selection of measurements of type mt.name
  meas.mt = meas.types[mt.name,]
  meas.mt.num = sum(meas.mt)
  num.types[mt.name] = meas.mt.num

  ##-- error on average like for PDG, assuming no correlation
  average.err = 1/sqrt(sum(1/meas.error[meas.mt]^2))
  error.max = 3*sqrt(meas.mt.num)*average.err

  ##-- keep only measurement with not too large errors as in PDG S-factor calculation
  meas.mt.keep = meas.mt & (meas.error <= error.max)
  meas.keep = meas.keep | meas.mt.keep

  ##
  ## chisq contribution from measurement of one type
  ##

  ##-- chisq when disregarding correlation between measurements of the same type
  meas.mt.chisq = sum(((meas[meas.mt.keep] - drop(quant.comb %*% quant)) / meas.error[meas.mt.keep])^2)

  ##-- chisq from measurements of the current type, including correlations
  meas.mt.chisq.corr = drop(
    t((meas - delta %*% quant)[meas.mt.keep])
    %*% solve(meas.cov[meas.mt.keep, meas.mt.keep])
    %*% (meas - delta %*% quant)[meas.mt.keep])
  
  ##-- collect chisq/dof for each type of measurement
  chisq.types[mt.name] = meas.mt.chisq.corr
  ##-- number of measurements of the current type that have not too-large errors
  num.types.keep[mt.name] = sum(meas.mt.keep)
  ##
  ## compute dof for kept measurements of the current type
  ##
  ## special treatment when there is just 1 kept measurement of a specific type
  ## when there is one measurement and one corresponding fitted quantity we expect
  ## a residual close to zero, because the fitted quantity converges to the measurement
  ## relatively undisturbed by chisq contributions due to correlation with other measurement
  ## however it can happen that even a single measurement of one type has a significant
  ## chisq contribution when it must match a combination of other measurements
  ## if there are many measurements of the type tau -> hhh nu with h = pi, K they will
  ## determine also the total tau -> hhh nu BR when h is not identified.  Even a single
  ## measurement of tau -> hhh undifferentiated can have a large chisq then.
  ## by setting dof = 1 rather than zero when a measurement type has a single measurement
  ## we acknowledge that because of correlations the chisq contribution can be different from zero
  ## and we even consider the possibility of applying an S-factor when the chisq contribution is large
  ## 
  dof.types[mt.name] = max(num.types[mt.name] -1, 1)

  ##-- print different chisq/dof
  if (FALSE) {
    cat(mt.name,
        format(meas.mt.chisq/dof.types[mt.name], width=15),
        format(meas.mt.chisq.corr/dof.types[mt.name], width=15),
        "\n")
  }

  save = options()
  options(digits=4)
  if (any(xor(meas.mt, meas.mt.keep))) {
    cat("\nS-factor calculation, exclude because error > 3*sqrt(N)*av_err=", error.max, "\n")
    excl = meas.error[xor(meas.mt, meas.mt.keep)]
    cat(paste(names(excl), "error=", format(excl, digits=4), "nsigma=", format(excl/error.max*3, digits=4), sep=" "), sep="\n")
  }
  if (FALSE && any(meas.mt.keep)) {
    cat("\nS-factor calculation, included measurements, 3*sqrt(N)*av_err=", error.max, "\n")
    excl = meas.error[meas.mt.keep]
    cat(paste(names(excl), "error=", format(excl, digits=4), "nsigma=", format(excl/error.max*3, digits=4), sep=" "), sep="\n")
  }
  options(save)

  tmp = sqrt(chisq.types[mt.name]/dof.types[mt.name])
  sfact.types[mt.name] = tmp
  sfact[meas.mt] = tmp
}

##-- show chisq differences due to correlation between same types measurements
if (FALSE && length(chisq.out)>0) {
  cat("Chisq differences due to correlation of measurements\n")
  show(chisq.out)
}
rm(chisq.out)

##-- recompute chisq and chisq/dof
dof = meas.num - quant.num + constr.num
chisq = drop(t(meas - delta %*% quant) %*% meas.invcov %*% (meas - delta %*% quant))

##
## PDG computes and uses S-factors only if chisq/dof > 1
## similarly, we compute and use S-factors only for the subset of measurements
## belonging to measurement types for which chisq/dof > 1
## 

##-- logical array of all measurements with S-factor > 1
meas.sfact = sfact > 1
##-- measurements with S>1 and not-too-large errors
meas.select = meas.keep & meas.sfact

##
## dof corresponding to selected measurement is the number of quantities
## that are related to such measurements. The delta matrix says which
## quantities relate to which measurements, use it to determine the
## quantities that are related to the selected measurements
##
quant.select.num = sum((rep(1, sum(meas.select)) %*% abs(delta[meas.select,])) > 0)
dof.select = sum(meas.select) - quant.select.num + constr.num

##++ prevent problems when all measurements have large errors but one
if (any(meas.select) && dof.select < 1) {
  cat("warning: too few measurements for S-factor calculation, setting dof=1\n")
  dof.select = 1
}

##-- save step 0
sfact.types.0 = sfact.types
sfact.0 = sfact

##-- chisq without S-factors, just for selected measurements
chisq.select.0 = drop(
  t((meas - delta %*% quant)[meas.select])
  %*% solve(meas.cov[meas.select, meas.select])
  %*% (meas - delta %*% quant)[meas.select])

##-- chisq with S-factors, just for selected measurements
chisq.select.1 = drop(
  t((meas - delta %*% quant)[meas.select])
  %*% solve(diag.m(sfact[meas.select]) %*% meas.cov[meas.select, meas.select] %*% diag.m(sfact[meas.select]))
  %*% (meas - delta %*% quant)[meas.select])

repeat {
  ##-- chisq with S-factors
  chisq.select = drop(
    t((meas - delta %*% quant)[meas.select])
    %*% solve(diag.m(sfact[meas.select]) %*% meas.cov[meas.select, meas.select] %*% diag.m(sfact[meas.select]))
    %*% (meas - delta %*% quant)[meas.select])
  
  ##-- adjust S-factors to obtain chisq/dof = 1, when restricted to selected measurements
  sfact[meas.select] = sfact[meas.select] * sqrt(chisq.select/dof.select)
  
  if (!any(meas.select)) break
  if (abs(chisq.select/dof.select -1) < 1e-6) break
}

##-- update sfact.types with adjusted sfact for selected meas
for (meas.s in names(meas.select[meas.select])) {
  sfact.types[names(which(meas.types[,meas.s]))] = sfact[meas.s]
}

##-- save step 1
sfact.types.1 = sfact.types
sfact.1 = sfact

##-- do not apply S-factors less than one
sfact.types = pmax(sfact.types, 1)
sfact = pmax(sfact, 1)

##
## inflate the covariance matrix with the S-factors
##
meas2.cov = diag.m(sfact) %*% meas.cov %*% diag.m(sfact)
rownames(meas2.cov) = meas.names
colnames(meas2.cov) = meas.names
meas2.invcov = solve(meas2.cov)

##-- compute new inflated fitted quantities covariance matrix 
##
## full covariance of measurements and constraint values
## - there is no error on the constraint values
##
full2.v.cov = rbind(
  cbind(meas2.cov, matrix(0, meas.num, constr.num, dimnames=list(NULL, constr.names))),
  cbind(matrix(0, constr.num, meas.num, dimnames=list(constr.names)), matrix(0, constr.num, constr.num)))

##-- covariance matrix of fitted values
quant2.constr.cov = solve.m %*% full2.v.cov %*% t(solve.m)
quant2.cov = quant2.constr.cov[1:quant.num, 1:quant.num, drop=FALSE]
## rownames(quant2.cov) = quant.names
## colnames(quant2.cov) = quant.names

##-- updated errors and correlations
quant2.err = sqrt(diag.m(quant2.cov))
quant2.corr = quant2.cov / (quant2.err %o% quant2.err)

##-- all measurements chisq after S-factor inflation
chisq2 = drop(t(meas - delta %*% quant) %*% meas2.invcov %*% (meas - delta %*% quant))

cat("\n")
cat("##\n")
cat("## S-factors accounting for larger than expected chi-square\n")
cat("##\n")

out = rbind(
  original = c(
    chisq=chisq, dof=dof,
    "chisq/dof"=chisq/dof,
    CL=pchisq(chisq, dof, lower.tail=FALSE)))
if (any(meas.select)) {
  out = rbind(out
    ,"  after S-factor, 2nd stage" = c(
       chisq=chisq2, dof=dof,
       "chisq/dof"=chisq2/dof,
       CL=pchisq(chisq2, dof, lower.tail=FALSE))
    ,"no-large-error" = c(
       chisq=chisq.select.0, dof=dof.select,
       "chisq/dof"=chisq.select.0/dof.select,
       CL=pchisq(chisq.select.0, dof.select, lower.tail=FALSE))
    ,"  after S-factor, 1st stage" = c(
       chisq=chisq.select.1, dof=dof.select,
       "chisq/dof"=chisq.select.1/dof.select,
       CL=pchisq(chisq.select.1, dof.select, lower.tail=FALSE))
    ,"  after S-factor, 2nd stage" = c(
       chisq=chisq.select, dof=dof.select,
       "chisq/dof"=chisq.select/dof.select,
       CL=pchisq(chisq.select, dof.select, lower.tail=FALSE)))
}
show(out)

cat("Averaged quantities: value, error, error with S-factor, S-factor\n") 
if (quant.num > 1) {
  ##-- if multiple average, use dedicated S-factor computation
  sfact.row = quant2.err / quant.err
  sfact0.row = sfact.types.0[quant.names]
  sfact1.row = sfact.types.1[quant.names]
  chisq.row = chisq.types[quant.names]
  dof.row = dof.types[quant.names]
} else {
  ##-- if averaging a single quantity, use global chisq to compute S-factor
  sfact0.row = sfact.types.0[quant.names]
  sfact1.row = sfact.types.1[quant.names]
  sfact.row = max(sqrt(chisq/dof), 1)
  ##-- also for single quantity, never deflate error
  quant2.err[1] = quant.err[1] * sfact.row
  chisq.row = chisq
  dof.row = dof
}
show(rbind(value=quant,
           error=quant.err,
           upd.error=quant2.err,
           "S-factor"=sfact.row,
           "S-factor_2"=sfact1.row,
           "S-factor_1"=sfact0.row,
           chisq=chisq.row,
           dof=dof.row
           ))

if (quant.num > 1) {
  cat("correlation\n") 
  show(quant.corr)
  cat("correlation, S-factor inflated\n") 
  show(quant2.corr)
}

##
## measurement types that do not correspond to an averaged quantity
## - compute average and errors
##
meas.names.extra = meas.types.names[!(meas.types.names %in% quant.names)]
if (length(meas.names.extra) >0) {
  meas.extra.id = subset(meas.types.id, meas.types.names %in% meas.names.extra)
  meas.extra = drop(meas.extra.id %*% quant[1:quant.num])
  meas.extra.err = sqrt(diag.m(meas.extra.id %*% quant.cov %*% t(meas.extra.id)))
  ##-- update with S-factors measurement types that do not correspond to an averaged quantity
  meas.extra.err.upd = sqrt(diag.m(meas.extra.id %*% quant2.cov %*% t(meas.extra.id)))
  
  cat("Non-averaged measurement types\n") 
  show(rbind(value=meas.extra,
             error=meas.extra.err,
             upd.error=meas.extra.err.upd,
             "S-factor"=meas.extra.err.upd/meas.extra.err,
             "S-factor_2"=sfact.types.1[meas.names.extra],
             "S-factor_1"=sfact.types.0[meas.names.extra],
             chisq=chisq.types[meas.names.extra],
             dof=dof.types[meas.names.extra]
           ))
}
cat("## end\n")

options(options.save)
##++ } ##-- end function alucomb

args <- commandArgs(TRUE)
if (length(args) > 0 && exists("alucomb")) alucomb(file = args[1]) 
