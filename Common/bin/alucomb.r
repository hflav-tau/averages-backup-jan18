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
  file.name = args[1]
} else {
  file.name = "average.input"
}

##++ alucomb = function(file.name = "") {

##-- set very large line width to print even large amount of averaged quantities on single line
options.save = options()
options(width=10000)

file.name.data = gsub("[.][^.]*$", "_alucomb.rdata", file.name)

rc = alucomb.read(file.name)
measurements = rc$measurements
combination = rc$combination
rm(rc)

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
## transform COMBOFQUANT and SUMOFQUANT cards
## (quantities that are combination of other quantities)
## into constraints consisting in requirements that
## a combination of values be equal to a constant numeric value
##
tmp = mapply(function(name, value) {
  rc = c(value, -1)
  names(rc)[length(rc)] = name
  rc
},
  names(combination$meas.lin.combs), combination$meas.lin.combs,
  SIMPLIFY=FALSE)
if (length(tmp) > 0) {
  names(tmp) = paste(names(tmp), "coq", sep=".")
  combination$constr.comb = c(combination$constr.comb, tmp)
  tmp2 = as.list(rep(0, length(tmp)))
  names(tmp2) = names(tmp)
  combination$constr.val = c(combination$constr.val, tmp2)
  rm(tmp2)
}
rm(tmp)

if (length(combination$constr.comb) > 1) {
  dupl.constr = NULL
  for (i in seq(1, length(combination$constr.comb))) {
    for (j in seq(i+1, length=length(combination$constr.comb)-i)) {
      if (length(setdiff(names(combination$constr.comb[[i]]), names(combination$constr.comb[[j]]))) == 0) {
        qn = names(combination$constr.comb[[i]])
        dupl.constr = rbind(dupl.constr,
          matrix(c(combination$constr.val[[i]], unlist(combination$constr.comb[[i]][qn])), nrow=1,
                 dimnames=list(names(combination$constr.comb[i]), c("val", qn))),
          matrix(c(combination$constr.val[[j]], unlist(combination$constr.comb[[j]][qn])), nrow=1,
                 dimnames=list(names(combination$constr.comb[j]), c("val", qn))))
      }
    }
  }
  if (length(dupl.constr) > 0) {
    cat("warning: duplicated constraints, listed in row pairs\n")
    show(dupl.constr)
  }
}

if (length(combination$constr.comb) > 0) {
  ##-- retain only constraints whose terms are all included in the fitted quantities
  constr.select = sapply(combination$constr.comb, function(x) all(names(x) %in% quant.names))
  if (any(!constr.select)) {
    cat("\nThe following constraints are dropped:\n")
    mapply(function(comb, val, val.name) {
      tmp = val
      names(tmp) = val.name
      show(c(comb, tmp))
    },
    combination$constr.comb[!constr.select],
    combination$constr.val[!constr.select],
    names(combination$constr.val[!constr.select]))
    cat("\nThe following measurement types are missing:\n")
    show(unique(unlist(lapply(combination$constr.comb, function(x) setdiff(names(x), quant.names)))))
  }
  combination$constr.comb = combination$constr.comb[constr.select]
  combination$constr.val = combination$constr.val[constr.select]
}

##-- quantity measured per measurement
meas.quantities = sapply(measurements, function(x) names(x$value))

meas.included.list = meas.quantities %in% quant.names
names(meas.included.list) = names(meas.quantities)[meas.included.list]
meas.names.discarded =  meas.names[!meas.included.list]
if (length(meas.names.discarded) >0) {
  cat("\nThe following measurements are discarded:\n")
  cat(paste("  ", meas.names.discarded, collapse="\n"), "\n");
}

##-- update
measurements = measurements[meas.included.list]
meas.names = names(measurements)
meas.num = length(measurements[meas.included.list])
meas.quantities = meas.quantities[meas.included.list]

larger = numeric()
slightly.larger = numeric()
not.matching = numeric()

##-- checks on syst terms
for (mn in meas.names) {
  ##-- check that the sum of syst. terms does not exceed the syst. error
  syst.contribs = sqrt(sum(measurements[[mn]]$syst.terms^2))
  if (syst.contribs > (1+1e-3)*measurements[[mn]]$syst) {
    larger = rbind(larger, matrix(c(measurements[[mn]]$syst, syst.contribs), 1, 2, dimnames=list(mn)))
  } else if (syst.contribs > (1+1e-5)*measurements[[mn]]$syst) {
    slightly.larger = rbind(slightly.larger, matrix(c(measurements[[mn]]$syst, syst.contribs), 1, 2, dimnames=list(mn)))
  }
  ##-- warning if sum of syst terms different w.r.t. total syst
  if (abs(syst.contribs-measurements[[mn]]$syst) >
      1e-3 * sqrt(syst.contribs * measurements[[mn]]$syst)) {
    not.matching = rbind(not.matching, matrix(c(measurements[[mn]]$syst, syst.contribs), 1, 2, dimnames=list(mn)))
  }
}

if (length(slightly.larger) > 0) {
  cat("\nwarning: syst. terms sum slightly larger than total syst. error\n")
  colnames(slightly.larger) = c("total", "sum of terms")
  show(slightly.larger)
}
if (length(larger) > 0) {
  cat("\nerror: syst. terms sum larger than total syst. error\n")
  colnames(larger) = c("total", "sum of terms")
  show(larger)
  stop("aborting")
}
if (length(not.matching) > 0) {
  cat("\nwarning: syst. terms do not match total syst. error, Combos requires that\n")
  colnames(not.matching) = c("total", "sum of terms")
  show(not.matching)
}

##
## shift measurements according to updated external parameter dependencies
## update systematic terms according to updated external parameter errors
##

for (mn in meas.names) {
  value.delta = numeric()
  syst.term.deltasq = numeric()
  for (param.upd in names(combination$params)) {
    for (param.orig in names(measurements[[mn]]$params)) {
      if (param.orig == param.upd && param.orig %in% names(measurements[[mn]]$syst.terms)) {
        ##-- collect differences due to updated parameters
        value.delta = c(value.delta,
          ((combination$params[[param.upd]]["value"] - measurements[[mn]]$params[[param.orig]]["value"])
           * measurements[[mn]]$syst.terms[param.orig] / measurements[[mn]]$params[[param.orig]]["delta_pos"]))
        ##-- update systematic contribution according to updated parameter
        syst.term.orig = measurements[[mn]]$syst.terms[param.orig]
        syst.term.upd = syst.term.orig *
          combination$params[[param.upd]]["delta_pos"] / measurements[[mn]]$params[[param.orig]]["delta_pos"]
        measurements[[mn]]$syst.terms[param.orig] = syst.term.upd
        ##-- collect difference of syst term squares, to adjust the total systematic error as well
        syst.term.deltasq = c(syst.term.deltasq, (syst.term.upd^2 - syst.term.orig^2))
        if (FALSE) {
        cat(format(measurements[[mn]]$tag,width=30),
            format(param.orig,width=15),
            format(c(
            measurements[[mn]]$params[[param.orig]]["value"],
            combination$params[[param.upd]]["value"],
            (combination$params[[param.upd]]["value"] - measurements[[mn]]$params[[param.orig]]["value"])
            / measurements[[mn]]$params[[param.orig]]["delta_pos"],
            ## (value.upd - value.orig),
            measurements[[mn]]$params[[param.orig]]["delta_pos"],
            combination$params[[param.upd]]["delta_pos"],
            combination$params[[param.upd]]["delta_pos"] / measurements[[mn]]$params[[param.orig]]["delta_pos"]),
            width=10,digits=4,scientific=TRUE),
            "\n")
        }
      }
    }
  }
  ##-- update value
  measurements[[mn]]$value.orig = measurements[[mn]]$value
  measurements[[mn]]$value = measurements[[mn]]$value + sum(value.delta)
  ##-- update systematic error
  measurements[[mn]]$syst.orig = measurements[[mn]]$syst
  measurements[[mn]]$syst = sqrt(measurements[[mn]]$syst^2 + sum(syst.term.deltasq))
}

##-- get list of values, stat errors, syst errors
meas.val = sapply(measurements, function(x) {tmp=x$value; names(tmp)=NULL; tmp})
meas.stat = sapply(measurements, function(x) {tmp=x$stat; names(tmp)=NULL; tmp})
meas.syst = sapply(measurements, function(x) {tmp=x$syst; names(tmp)=NULL; tmp})
meas.err = sqrt(meas.stat^2 + meas.syst^2)

##--- unshifted values
meas.val.orig = sapply(measurements, function(x) {tmp=x$value.orig; names(tmp)=NULL; tmp})
meas.syst.orig = sapply(measurements, function(x) {tmp=x$syst.orig; names(tmp)=NULL; tmp})

##-- which measurements got shifted in value or syst. error
meas.shifted = (meas.val.orig != meas.val) | (meas.syst.orig != meas.syst)

if (any(meas.shifted)) {
  cat("\nThe following measurements were shifted from updated external parameters\n")
  show(rbind(orig=meas.val.orig[meas.shifted],
             value=meas.val[meas.shifted],
             stat=meas.stat[meas.shifted],
             orig=meas.syst.orig[meas.shifted],
             syst=meas.syst[meas.shifted]))
}
  
##-- for each syst. term external parameter, collect affected measurements
syst.terms.list = list()
for (mn in meas.names) {
  for (syst.term in names(measurements[[mn]]$syst.terms)) {
    ##-- add measurement to the list of the currect syst. term
    syst.terms.list[[syst.term]] = c(syst.terms.list[[syst.term]], mn)
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

flag = FALSE

##-- check that the STAT_CORRELATED_WITH terms are symmetric
if (any(meas.corr.stat != t(meas.corr.stat))) {
  errors = character(0)
  for (mi in 1:meas.num) {
    for (mj in seq(mi+1, length=meas.num-mi)) {
      if (meas.corr.stat[mi, mj] != meas.corr.stat[mj, mi]) {
        errors = c(errors, paste(meas.names[mi], " - ", meas.names[mj],
          " : ", meas.corr.stat[mi, mj], " , ",  meas.corr.stat[mj, mi], sep=""))
        flag = TRUE
      }
    }
  }
  cat("error: asymmetric statistical correlation\n  ", paste(errors, collapse="\n  "), "\n", sep="")
}

##-- check that the TOTAL_CORRELATED_WITH terms are symmetric
if (any(meas.corr != t(meas.corr))) {
  errors = character(0)
  for (mi in 1:meas.num) {
    for (mj in seq(mi+1, length=meas.num-mi)) {
      if (meas.corr[mi, mj] != meas.corr[mj, mi]) {
        errors = c(errors, paste(meas.names[mi], " - ", meas.names[mj],
          " : ", meas.corr[mi, mj], " , ",  meas.corr[mj, mi], sep=""))
        flag = TRUE
      }
    }
  }
  cat("error: asymmetric total correlation\n  ", paste(errors, collapse="\n  "), "\n", sep="")
}

if (flag) stop("quitting")

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
meas.cov = meas.corr * (meas.err %o% meas.err)
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
## - quantities are the results of the fit
## measurements are linear combinations of quantities, meas_i = delta_ij * quant_j
## if a measurement i corresponds to a quantity j then delta_ij = 1
## all remaining delta matrix terms are zero
## one can generalize the above concepts to measurements that are linear combinations
## of quantities by using proper coefficients different from 1
##
delta = sapply(quant.names, function(x) as.numeric(x == meas.quantities))
rownames(delta) = meas.names

##-- print corrected measurements
if (TRUE) {
  cat("\n##\n")
  cat("## Using the following measurements\n")
  cat("##\n")
  show(rbind(value=meas.val,
             stat=meas.stat,
             syst=meas.syst,
             error=meas.err))
}

if (FALSE && !flag.no.maxLik) {
##
## solve for quantities with iterative chi-square minimization
## numerical minimization needs to be updated to handle constraints
##

##-- compute weight matrix for computing chisq
meas.invcov = solve(meas.cov)

logLik.average = function(par) {
  chisq = t(meas.val - delta %*% par) %*% meas.invcov %*% (meas.val - delta %*% par)
  return(-1/2*chisq)
}

fit = maxLik(logLik.average, start=rep(0,quant.num), method="BFGS")

quant.val = coef(fit)
names(quant.val) = quant.names
quant.cov = vcov(fit)
rownames(quant.cov) = quant.names
colnames(quant.cov) = quant.names
quant.err = sqrt(diag.m(quant.cov))
quant.corr = quant.cov / (quant.err %o% quant.err)

chisq = drop(t(meas.val - delta %*% quant.val) %*% meas.invcov %*% (meas.val - delta %*% quant.val))
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
show(rbind(value=quant.val[1:quant.num], error=quant.err[1:quant.num]))
if (quant.num > 1) {
  cat("correlation\n")
  ##++ show(quant.corr[1:quant.num, 1:quant.num])
}
} ## !flag.no.maxLik && FALSE

##
## analytical minimum chi-square solution for quantities
##
quant.val = rep(0, quant.num)
names(quant.val) = quant.names
meas.invcov = solve(meas.cov)
quant.invcov = t(delta) %*% meas.invcov %*% delta

##-- constraints
constr.num = length(combination$constr.comb)
constr.names = names(combination$constr.comb)
constr.m =  do.call(rbind, lapply(combination$constr.comb, function(x) {tmp = quant.val; tmp[names(x)] = x; tmp}))
constr.v = unlist(combination$constr.val)

if (constr.num > 0) {
  ##-- determine the typical size of quant.invcov elements
  sv = svd(quant.invcov)$d
  sv.central = round(quant.num*1/3):round(quant.num*2/3)
  sv.log.mean = mean(log(sv[sv.central]))
  quant.invcov.order = 10^round(sv.log.mean/log(10))

  ##-- to avoid computationally singular matrix, apply proper factor to constraint equations
  constr.m.order = 10^round(log(mean(abs(constr.m[constr.m!=0])))/log(10))
  constr.m = constr.m * quant.invcov.order/constr.m.order
  constr.v = constr.v * quant.invcov.order/constr.m.order

  cat("\n## Constraint equations\n")
  tmp = mapply(function(name, val, comb) {names(val) = name; cat("\n"); show(c(val, unlist(comb)))},
    names(constr.v), constr.v/quant.invcov.order*constr.m.order,
    apply(constr.m/quant.invcov.order*constr.m.order, 1, function(x) list(x[x!=0])))

  ##-- build full matrix in front of c(quant.val vector, lagr.mult. vector)
  full.m = rbind(
    cbind(quant.invcov, t(constr.m)),
    cbind(constr.m, matrix(0, constr.num, constr.num)))
  ##-- build full matrix in front of c(meas.val, constraint values)
  full.v.m = rbind(
    cbind(t(delta) %*% meas.invcov, matrix(0, quant.num, constr.num, dimnames=list(NULL, constr.names))),
    cbind(matrix(0, constr.num, meas.num, dimnames=list(constr.names)), diag.m(rep(1, constr.num))))
  ##-- build full vector c(measurements vector, constraint values)
  full.v = c(meas.val, constr.v)
} else {
  full.m = quant.invcov
  full.v.m = t(delta) %*% meas.invcov
  full.v = meas.val
}

##-- matrix that applied to c(meas.val, constr. values) gives c(quant.val, lagr.mult.)
solve.m = solve(full.m) %*% full.v.m

##-- solve for both quantities and lagrange multipliers
## quant.constr.val = solve(full.m, (full.v.m %*% full.v))
quant.constr.val = solve.m %*% full.v
quant.val = quant.constr.val[1:quant.num, drop=FALSE]
names(quant.val) = quant.names

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
quant.corr = quant.cov / (quant.err %o% quant.err)

chisq = drop(t(meas.val - delta %*% quant.val) %*% meas.invcov %*% (meas.val - delta %*% quant.val))

cat("\n")
cat("##\n")
cat("## exact solution, chisq/d.o.f. = ",chisq, "/", meas.num - quant.num,
    ", CL = ", (1-pchisq(chisq, df=meas.num-quant.num)), "\n",sep="")
cat("##\n")
show(rbind(value=quant.val[1:quant.num], error=quant.err[1:quant.num]))
if (quant.num > 1) {
  cat("correlation\n")
  show(quant.corr[1:quant.num, 1:quant.num])
}
cat("## end\n")

##-- chisq contribution, dof, S-factor for each measurement type (now same as quant.val)
quant.chisq = quant.val*0
quant.nums = quant.val*0
quant.nums.keep = quant.val*0
quant.sfact = quant.chisq + 1

##-- S-factor for each measurement
meas.sfact = meas.val*0 + 1
meas.keep = meas.val & FALSE

##
## compute S-factors
## the following is filled:
## - quant.chisq: chisq per measurement type
## - quant.sfact: S-factor per measurement type
## - meas.sfact: S-factor per measurement
##
for (mt.name in quant.names) {
  ##-- selection of measurements of type mt.name
  meas.mt = delta[, mt.name] != 0
  meas.mt.num = sum(meas.mt)
  quant.nums[mt.name] = meas.mt.num

  ##-- compute max error to retain measurement for S-factor calculation
  error.max = 3*sqrt(meas.mt.num)*quant.err[mt.name]

  ##-- keep only measurement with not too large errors as in PDG S-factor calculation
  meas.mt.keep = meas.mt & (meas.err < error.max)
  meas.mt.keep.num = sum(meas.mt.keep)
  quant.nums.keep[mt.name] = meas.mt.keep.num
  meas.keep = meas.keep | meas.mt.keep

  ##
  ## chisq contribution from measurement of one type
  ##

  ##-- chisq when disregarding correlation between measurements of the same type
  meas.mt.chisq = sum(((meas.val[meas.mt] - quant.val[mt.name]) / meas.err[meas.mt])^2)

  ##-- chisq from measurements of the current type, including correlations
  if (meas.mt.num > 0) {
    meas.mt.chisq.corr = drop(
      t((meas.val - delta %*% quant.val)[meas.mt])
      %*% solve(meas.cov[meas.mt, meas.mt])
      %*% (meas.val - delta %*% quant.val)[meas.mt])
  } else {
    meas.mt.chisq.corr = 0
  }

  ##-- collect chisq contribution for each type of measurement
  quant.chisq[mt.name] = meas.mt.chisq.corr

  ##
  ## print measurement with too large errors
  ##
  save = options()
  options(digits=4)
  if (any(xor(meas.mt, meas.mt.keep))) {
    cat("\nS-factor calculation, exclude because error > 3*sqrt(N)*av_err=", error.max, "\n")
    excl = meas.err[xor(meas.mt, meas.mt.keep)]
    cat(paste(names(excl), "error=", format(excl, digits=4), "nsigma=", format(excl/error.max*3, digits=4), sep=" "), sep="\n")
  }
  if (FALSE && any(meas.mt.keep)) {
    cat("\nS-factor calculation, included measurements, 3*sqrt(N)*av_err=", error.max, "\n")
    excl = meas.err[meas.mt.keep]
    cat(paste(names(excl), "error=", format(excl, digits=4), "nsigma=", format(excl/error.max*3, digits=4), sep=" "), sep="\n")
  }
  options(save)

  ##
  ## PDG S-factor for fits
  ## p.13 of http://pdg.lbl.gov/2009/reviews/rpp2009-rev-rpp-intro.pdf
  ##
  if (meas.mt.keep.num >0) {
    res = meas.val[meas.mt.keep] - quant.val[mt.name]
    dr = abs(res) / sqrt(meas.err[meas.mt.keep] * quant.err[mt.name])
    res.sq.exp = meas.err[meas.mt.keep]^2 - quant.err[mt.name]^2
    if (any(res.sq.exp < 0)) {
      show(mt.name)
      show(res.sq.exp[res.sq.exp < 0])
      cat("warning: measurement on fitted quantity is larger than error on one of its measurements\n")
      ## stop("measurement on fitted quantity is larger than error on one of its measurements")
    }
    ## tmp = sqrt(mean(ifelse(dr > 1e-6, res^2 / res.sq.exp, 1)))
    ##--- protect against negative expected error
    tmp = sqrt(mean(ifelse(dr > 1e-6 & res.sq.exp > 0, res^2 / res.sq.exp, 1)))
  } else {
    tmp = 1
  }

  quant.sfact[mt.name] = tmp
  meas.sfact[meas.mt] = tmp
}

##-- recompute chisq and chisq/dof
dof = meas.num - quant.num + constr.num
chisq = drop(t(meas.val - delta %*% quant.val) %*% meas.invcov %*% (meas.val - delta %*% quant.val))

##-- do not apply S-factors less than one
quant.sfact = pmax(quant.sfact, 1)
meas.sfact = pmax(meas.sfact, 1)

##-- chisq with S-factors, just for not-too-large error measurements
chisq.keep = drop(
  t((meas.val - delta %*% quant.val)[meas.keep])
  %*% solve(diag.m(meas.sfact[meas.keep]) %*% meas.cov[meas.keep, meas.keep] %*% diag.m(meas.sfact[meas.keep]))
  %*% (meas.val - delta %*% quant.val)[meas.keep])

##++ quant.keep.num = sum((rep(1, sum(meas.keep)) %*% abs(delta[meas.keep,])) > 0)
##++ assume there are no measurement types left without measurements
quant.keep.num = quant.num - constr.num
dof.keep = sum(meas.keep) - quant.keep.num

##
## inflate the covariance matrix with the S-factors
##
meas2.cov = diag.m(meas.sfact) %*% meas.cov %*% diag.m(meas.sfact)
rownames(meas2.cov) = meas.names
colnames(meas2.cov) = meas.names
meas2.invcov = solve(meas2.cov)

##-- compute new inflated fitted quantities covariance matrix 
full2.v.cov = rbind(
  cbind(meas2.cov, matrix(0, meas.num, constr.num, dimnames=list(NULL, constr.names))),
  cbind(matrix(0, constr.num, meas.num, dimnames=list(constr.names)), matrix(0, constr.num, constr.num)))

##-- covariance matrix of fitted values
quant2.constr.cov = solve.m %*% full2.v.cov %*% t(solve.m)
quant2.cov = quant2.constr.cov[1:quant.num, 1:quant.num, drop=FALSE]

##-- updated errors and correlations
quant2.err = sqrt(diag.m(quant2.cov))
quant2.corr = quant2.cov / (quant2.err %o% quant2.err)

##-- all measurements chisq after S-factor inflation
chisq2 = drop(t(meas.val - delta %*% quant.val) %*% meas2.invcov %*% (meas.val - delta %*% quant.val))

##
## new S-factors computed after 1st S-factor inflation
##
if (quant.num > 1) {
  ##-- if multiple average, use dedicated S-factor computation
  quant2.sfact = quant2.err / quant.err
} else {
  ##-- if averaging a single quantity, use global chisq to compute S-factor
  quant2.sfact = max(sqrt(chisq/dof), 1)
}

##
## final S-factors
## - use 1st stage S-factors to inflate errors
## - recompute errors and S-factors
##
quant2.sfact = pmax(quant2.sfact, 1)
meas2.sfact = drop(delta %*% quant2.sfact)

##-- chisq with S-factors, just for not-too-large error measurements
chisq2.keep = drop(
  t((meas.val - delta %*% quant.val)[meas.keep])
  %*% solve(diag.m(meas2.sfact[meas.keep]) %*% meas.cov[meas.keep, meas.keep] %*% diag.m(meas2.sfact[meas.keep]))
  %*% (meas.val - delta %*% quant.val)[meas.keep])

dof2.keep = dof.keep

if (quant.num > 1) {
  ##
  ## inflate the covariance matrix with the S-factors
  ##
  meas3.cov = diag.m(meas2.sfact) %*% meas.cov %*% diag.m(meas2.sfact)
  rownames(meas3.cov) = meas.names
  colnames(meas3.cov) = meas.names
  meas3.invcov = solve(meas3.cov)

  ##-- compute new inflated fitted quantities covariance matrix 
  full3.v.cov = rbind(
    cbind(meas3.cov, matrix(0, meas.num, constr.num, dimnames=list(NULL, constr.names))),
    cbind(matrix(0, constr.num, meas.num, dimnames=list(constr.names)), matrix(0, constr.num, constr.num)))
  
  ##-- covariance matrix of fitted values
  quant3.constr.cov = solve.m %*% full3.v.cov %*% t(solve.m)
  quant3.cov = quant3.constr.cov[1:quant.num, 1:quant.num, drop=FALSE]
  
  ##-- updated errors and correlations
  quant3.err = sqrt(diag.m(quant3.cov))
  quant3.corr = quant3.cov / (quant3.err %o% quant3.err)
  
  ##-- all measurements chisq after S-factor inflation
  chisq3 = drop(t(meas.val - delta %*% quant.val) %*% meas3.invcov %*% (meas.val - delta %*% quant.val))
} else {
  quant3.cov = quant2.cov
  quant3.err = quant2.err
  quant3.corr = quant2.corr
  chisq3 = chisq2
}

##
## final S-factors computed after 1st S-factor inflation
##
if (quant.num > 1) {
  ##-- if multiple average, use dedicated S-factor computation
  quant3.sfact = quant3.err / quant.err
} else {
  ##-- if averaging a single quantity, use global chisq to compute S-factor
  quant3.sfact = quant2.sfact
}

quant.sf.err = quant3.err
quant.sf.cov = quant3.cov
quant.sf.corr = quant3.corr
quant.sf.sfact = quant3.sfact
rm(quant3.err)
rm(quant3.cov)
rm(quant3.corr)
rm(quant3.sfact)

cat("\n")
cat("##\n")
cat("## S-factors accounting for larger than expected chi-square\n")
cat("##\n")

out = rbind(
  "fit.1" = c(
    chisq=chisq, dof=dof,
    "chisq/dof"=chisq/dof,
    CL=pchisq(chisq, dof, lower.tail=FALSE)))
if (any(meas.keep)) {
  out = rbind(out
    ,"fit.sf.1" = c(
       chisq=chisq.keep, dof=dof.keep,
       "chisq/dof"=chisq.keep/dof.keep,
       CL=pchisq(chisq.keep, dof.keep, lower.tail=FALSE))
    ,"fit.sf" = c(
       chisq=chisq2.keep, dof=dof2.keep,
       "chisq/dof"=chisq2.keep/dof2.keep,
       CL=pchisq(chisq2.keep, dof2.keep, lower.tail=FALSE)))
}
show(out)

cat("Averaged quantities: value, error, error with S-factor, S-factor\n") 
show(rbind(value=quant.val,
           error=quant.err,
           error.sf.2=quant2.err,
           error.sf=quant.sf.err,
           "S-factor.1"=quant.sfact,
           "S-factor.2"=quant2.sfact,
           "S-factor"=quant.sf.sfact
           ))

if (quant.num > 1) {
  cat("correlation\n") 
  show(quant.corr)
  cat("correlation, S-factor inflated\n") 
  show(quant.sf.corr)
}

##--- save data and results
save(file=file.name.data,
     measurements, combination,
     chisq, chisq.keep, chisq2.keep,
     dof, dof.keep, dof2.keep,
     meas.val,   meas.err,   meas.cov,  meas.cov.stat, meas.cov.syst, meas.corr,
     quant.val,  quant.err,  quant.cov,  quant.corr,  quant.sfact, 
                 quant.sf.err, quant.sf.cov, quant.sf.corr, quant.sf.sfact,
     constr.m, constr.v)

cat("\n")
cat(paste("file", file.name.data, "produced\n"))
cat("\n")
cat("## end\n")

options(options.save)
##++ } ##-- end function alucomb

args <- commandArgs(TRUE)
if (length(args) > 0 && exists("alucomb")) alucomb(file = args[1]) 
