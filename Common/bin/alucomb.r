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

##+++ alucomb = function(file.name = "") {

##--- set very large line width to print even large amount of averaged quantities on single line
options.save = options()
## options(width=10000)

file.name.data = gsub("[.][^.]*$", "_alucomb.rdata", file.name)

rc = alucomb.read(file.name)
measurements = rc$measurements
combination = rc$combination
rm(rc)

##
## build list of all measurements mentioned in the COMBINE section
## (right now we only handle COMBINE * * *)
##

##--- all measurements in the input cards
meas.names = names(measurements)
##--- quantities to be averaged
quant.names = combination$quantities

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
  ##--- retain only constraints whose terms are all included in the fitted quantities
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

##--- quantity measured per measurement
meas.quantities = sapply(measurements, function(x) names(x$value))

##--- discard measurements that are not associated to a declared fitted quantity
meas.included.list = meas.quantities %in% quant.names
names(meas.included.list) = names(meas.quantities)[meas.included.list]
meas.names.discarded =  meas.names[!meas.included.list]
if (length(meas.names.discarded) >0) {
  cat("\nwarning: the following measurements are discarded:\n")
  cat(paste("  ", meas.names.discarded, collapse="\n"), "\n");
}

##--- quantities involved in constraints
constr.quantities = unique(unlist(lapply(combination$constr.comb, function(x) names(x)), use.names=FALSE))

##--- discard fitted quantities that are not defined by measurements and constraints
quant.discarded = !(quant.names %in% c(meas.quantities, constr.quantities))
if (any(quant.discarded)) {
  cat("\nwarning: the following fitted quantities are discarded:\n")
  cat(paste("  ", quant.names[quant.discarded], collapse="\n"), "\n");
  combination$quantities = combination$quantities[!quant.discarded]
  quant.names = combination$quantities
}

##--- update
measurements = measurements[meas.included.list]
meas.names = names(measurements)
meas.num = length(measurements[meas.included.list])
meas.quantities = meas.quantities[meas.included.list]
quant.num = length(quant.names)

larger = numeric()
slightly.larger = numeric()
not.matching = numeric()

##--- checks on syst terms
for (mn in meas.names) {
  ##--- check that the sum of syst. terms does not exceed the syst. error
  syst.contribs = sqrt(sum(measurements[[mn]]$syst.terms^2))
  if (syst.contribs > (1+1e-3)*measurements[[mn]]$syst) {
    larger = rbind(larger, matrix(c(measurements[[mn]]$syst, syst.contribs), 1, 2, dimnames=list(mn)))
  } else if (syst.contribs > (1+1e-5)*measurements[[mn]]$syst) {
    slightly.larger = rbind(slightly.larger, matrix(c(measurements[[mn]]$syst, syst.contribs), 1, 2, dimnames=list(mn)))
  }
  ##--- warning if sum of syst terms different w.r.t. total syst
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
## scale all measurements of the specified type according to the scale parameter in the cards
##
quant.sfact.list = lapply(combination$quantities.options, function(el) { el["scale"] })
meas.sfact.cards = rep(1, meas.num)
if (length(quant.sfact.list) > 0) {
  names(meas.sfact.cards) = meas.names
  cat("\n")
  rc = lapply(names(quant.sfact.list), function(el) {
    sel = which(meas.quantities == el)
    rc = lapply(sel, function(i) {
      measurements[[i]]$stat <<- measurements[[i]]$stat * quant.sfact.list[[el]]
      measurements[[i]]$syst <<- measurements[[i]]$syst * quant.sfact.list[[el]]
      measurements[[i]]$syst.terms <<- measurements[[i]]$syst.terms * quant.sfact.list[[el]]
    })
    meas.sfact.cards[sel] <<- quant.sfact.list[[el]]
    cat("applying s-factor = ", quant.sfact.list[[el]], " for quantity ", el, "for measurements:\n")
    show(names(meas.sfact.cards[sel]))
  })
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
        ##--- collect differences due to updated parameters
        value.delta = c(value.delta,
          ((combination$params[[param.upd]]["value"] - measurements[[mn]]$params[[param.orig]]["value"])
           * measurements[[mn]]$syst.terms[param.orig] / measurements[[mn]]$params[[param.orig]]["delta_pos"]))
        ##--- update systematic contribution according to updated parameter
        syst.term.orig = measurements[[mn]]$syst.terms[param.orig]
        syst.term.upd = syst.term.orig *
          combination$params[[param.upd]]["delta_pos"] / measurements[[mn]]$params[[param.orig]]["delta_pos"]
        measurements[[mn]]$syst.terms[param.orig] = syst.term.upd
        ##--- collect difference of syst term squares, to adjust the total systematic error as well
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
  ##--- update value
  measurements[[mn]]$value.orig = measurements[[mn]]$value
  measurements[[mn]]$value = measurements[[mn]]$value + sum(value.delta)
  ##--- update systematic error
  measurements[[mn]]$syst.orig = measurements[[mn]]$syst
  measurements[[mn]]$syst = sqrt(measurements[[mn]]$syst^2 + sum(syst.term.deltasq))
}

##--- get list of values, stat errors, syst errors
meas.val = sapply(measurements, function(x) {tmp=x$value; names(tmp)=NULL; tmp})
meas.stat = sapply(measurements, function(x) {tmp=x$stat; names(tmp)=NULL; tmp})
meas.syst = sapply(measurements, function(x) {tmp=x$syst; names(tmp)=NULL; tmp})
meas.err = sqrt(meas.stat^2 + meas.syst^2)

##--- unshifted values
meas.val.orig = sapply(measurements, function(x) {tmp=x$value.orig; names(tmp)=NULL; tmp})
meas.syst.orig = sapply(measurements, function(x) {tmp=x$syst.orig; names(tmp)=NULL; tmp})

##--- which measurements got shifted in value or syst. error
meas.shifted = (meas.val.orig != meas.val) | (meas.syst.orig != meas.syst)

if (any(meas.shifted)) {
  cat("\nThe following measurements were shifted from updated external parameters\n")
  show(rbind(orig=meas.val.orig[meas.shifted],
             value=meas.val[meas.shifted],
             stat=meas.stat[meas.shifted],
             orig=meas.syst.orig[meas.shifted],
             syst=meas.syst[meas.shifted]))
}
  
##--- for each syst. term external parameter, collect affected measurements
syst.terms.list = list()
for (mn in meas.names) {
  for (syst.term in names(measurements[[mn]]$syst.terms)) {
    ##--- add measurement to the list of the currect syst. term
    syst.terms.list[[syst.term]] = c(syst.terms.list[[syst.term]], mn)
  }
}

##--- retain just the syst. contributions that affect at least two measurements
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
##--- set off-diagonal statistical correlation matrix coefficients from cards
for (mi.name in meas.names) {
  for (mj.name in intersect(names(measurements[[mi.name]]$corr.terms), meas.names)) {
    meas.corr.stat[meas.names %in% mi.name, meas.names %in% mj.name] = measurements[[mi.name]]$corr.terms[[mj.name]]
  }
  for (mj.name in intersect(names(measurements[[mi.name]]$corr.terms.tot), meas.names)) {
    meas.corr[meas.names %in% mi.name, meas.names %in% mj.name] = measurements[[mi.name]]$corr.terms.tot[[mj.name]]
  }
}

flag = FALSE

##--- check that the STAT_CORRELATED_WITH terms are symmetric
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

##--- check that the TOTAL_CORRELATED_WITH terms are symmetric
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

##--- not handled and forbidden to enter both total and stat. only correlations
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
    ##--- no addition needed for on-diagonal terms
    if (meas.i == meas.j) next
    syst.j = measurements[[meas.j]]$syst.terms
    ##--- systematics common to the two measurements
    correl.i.j = intersect(names(syst.i), names(syst.j))
    ##--- remove syst. terms uncorrelated to two different measurements
    correl.i.j = intersect(correl.i.j, syst.terms.corr)
    if (length(correl.i.j) == 0) next
    meas.cov.syst[meas.i,meas.j] = sum(syst.i[correl.i.j] * syst.j[correl.i.j])
  }
}

##--- if total correlation specified, get stat. correlation by subtraction
meas.cov.stat = ifelse(meas.cov == 0, meas.cov.stat, meas.cov - meas.cov.syst)
meas.cov.stat = meas.cov.stat + diag.m(meas.stat^2)
##--- total covariance
meas.cov.syst = meas.cov.syst + diag.m(meas.syst^2)
meas.cov = meas.cov.stat + meas.cov.syst

##--- total correlation
meas.corr = meas.cov / (meas.err %o% meas.err)

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

##--- print corrected measurements
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

##--- compute weight matrix for computing chisq
meas.invcov = solve(meas.cov)
meas.invcov = (meas.invcov + t(meas.invcov))/2

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
  show(quant.corr[1:quant.num, 1:quant.num])
}
} ## !flag.no.maxLik && FALSE

##
## analytical minimum chi-square solution for quantities
##
quant.val = rep(0, quant.num)
names(quant.val) = quant.names
meas.invcov = solve(meas.cov)
meas.invcov = (meas.invcov + t(meas.invcov))/2
##--- useful in calculations but not actual quant.invcov if there are constraints
quant.invcov = t(delta) %*% meas.invcov %*% delta
quant.invcov = (quant.invcov + t(quant.invcov))/2

##
## obtain constraint equations
##
constr.num = length(combination$constr.comb)
constr.names = names(combination$constr.comb)
constr.m =  do.call(rbind, lapply(combination$constr.comb, function(x) {tmp = quant.val; tmp[names(x)] = x; tmp}))
constr.v = unlist(combination$constr.val)

if (constr.num > 0) {
  ##--- determine the typical size of quant.invcov elements
  sv = svd(quant.invcov)$d
  sv.central = round(quant.num*1/3):round(quant.num*2/3)
  sv.log.mean = mean(log(sv[sv.central]))
  quant.invcov.order = 10^round(sv.log.mean/log(10))

  ##--- to avoid computationally singular matrix, apply proper factor to constraint equations
  constr.m.order = 10^round(log(mean(abs(constr.m[constr.m!=0])))/log(10))
  constr.m = constr.m * quant.invcov.order/constr.m.order
  constr.v = constr.v * quant.invcov.order/constr.m.order

  cat("\n## Constraint equations\n")
  tmp = mapply(function(name, val, comb) {names(val) = name; cat("\n"); show(c(val, unlist(comb)))},
    names(constr.v), constr.v/quant.invcov.order*constr.m.order,
    apply(constr.m/quant.invcov.order*constr.m.order, 1, function(x) list(x[x!=0])))
}

##
## if there are constraints, assemble full matrix equation
##
if (constr.num > 0) {
  ##--- build full matrix in front of c(quant.val vector, lagr.mult. vector)
  full.m = rbind(
    cbind(quant.invcov, t(constr.m)),
    cbind(constr.m, matrix(0, constr.num, constr.num)))
  ##--- build full matrix in front of c(meas.val, constraint values)
  full.v.m = rbind(
    cbind(t(delta) %*% meas.invcov, matrix(0, quant.num, constr.num, dimnames=list(NULL, constr.names))),
    cbind(matrix(0, constr.num, meas.num, dimnames=list(constr.names)), diag(constr.num)))
  ##--- build full vector c(measurements vector, constraint values)
  full.v = c(meas.val, constr.v)
} else {
  full.m = quant.invcov
  full.v.m = t(delta) %*% meas.invcov
  full.v = meas.val
}

##--- matrix that applied to c(meas.val, constr. values) gives c(quant.val, lagr.mult.)
solve.m = solve(full.m) %*% full.v.m

##--- solve for both quantities and lagrange multipliers
## quant.constr.val = solve(full.m, (full.v.m %*% full.v))
quant.constr.val = solve.m %*% full.v
quant.val = drop(quant.constr.val)[1:quant.num]
quant.cov = solve.m[1:quant.num,1:meas.num, drop=FALSE] %*% meas.cov %*% t(solve.m[1:quant.num,1:meas.num, drop=FALSE])
quant.cov = (quant.cov + t(quant.cov))/2
quant.err = sqrt(diag(quant.cov))
quant.corr = quant.cov / (quant.err %o% quant.err)

chisq = drop(t(meas.val - delta %*% quant.val) %*% meas.invcov %*% (meas.val - delta %*% quant.val))
dof = meas.num - quant.num + constr.num

cat("\n")
cat("##\n")
cat("## exact solution, chisq/d.o.f. = ",chisq, "/", dof, ", CL = ", (1-pchisq(chisq, df=dof)), "\n",sep="")
cat("##\n")
show(rbind(value=quant.val[1:quant.num], error=quant.err[1:quant.num]))
if (FALSE && quant.num > 1) {
  cat("correlation\n")
  show(quant.corr[1:quant.num, 1:quant.num])
}
cat("## end\n")

##
## use S-factor inflated measurement covariance to update fitted quantities covariance
##
get.fit.errors = function(meas.cov, meas.keep) {
  rc = list()
  rc$meas.cov = meas.cov
  rc$meas.keep = meas.keep

  meas.invcov = solve(meas.cov)
  meas.invcov = (meas.invcov + t(meas.invcov))/2
  
  rc$quant.cov = solve.m[1:quant.num,1:meas.num, drop=FALSE] %*% meas.cov %*% t(solve.m[1:quant.num,1:meas.num, drop=FALSE])
  rc$quant.cov = (rc$quant.cov + t(rc$quant.cov))/2

  rc$quant.err = sqrt(diag(rc$quant.cov))
  rc$quant.corr = rc$quant.cov / (rc$quant.err %o% rc$quant.err)
  rc$chisq = drop(t(meas.val - delta %*% quant.val) %*% meas.invcov %*% (meas.val - delta %*% quant.val))
  rc$dof = dof
  if (any(meas.keep)) {
    rc$chisq.keep = drop(
      t((meas.val - delta %*% quant.val)[meas.keep])
      %*% meas.invcov[meas.keep, meas.keep]
      %*% (meas.val - delta %*% quant.val)[meas.keep])
    rc$dof.keep = sum(meas.keep) - (quant.num - constr.num)
  }
  rc$quant.sfact = rc$quant.err / quant.err
  
  return(rc)
}

##
## S-factors, no grouping (i.e. all measurements together)
##

##--- computational tolerance
tol = sqrt(.Machine$double.eps)

##--- compute covariance matrix for pull averages
## pull.cov = meas.cov + (delta %*% quant.cov %*% t(delta)) -
##   delta %*% solve.m[1:quant.num,1:meas.num] %*% meas.cov -
##   meas.cov %*% t(solve.m[1:quant.num,1:meas.num]) %*% t(delta)
pull.cov = meas.cov - delta %*% quant.cov %*% t(delta)

##--- square root of pull weight matrix, to get pull direction
pull.invcov.sqrt = alu.matr.inv.sqrt.symm.semipos.norm(pull.cov, meas.err)
dof.pull = attr(pull.invcov.sqrt, "pos.eigen.num")
##--- pull
meas.pull = drop(pull.invcov.sqrt %*% (meas.val - (delta %*% quant.val)))
##--- force to zero very small pulls
meas.pull[abs(meas.pull) <= tol] = 0

##--- alu prescription using full pull vector term-by-term normalized with effective dof
meas.dof.eff = attr(pull.invcov.sqrt, "dof.eff")
meas.pull.weight = sapply(meas.dof.eff, function(x) ifelse(x!=0, 1/sqrt(x), 0))
meas.pull.dof = abs(meas.pull) * meas.pull.weight
meas.sfact.alu.full = ifelse(meas.pull.dof > 1, meas.pull.dof, 1)

##--- Orin Dahl S-factors prescription, all measurements together
meas.pull.sq = sum(meas.pull^2)
if (meas.pull.sq != 0) {
  meas.pull.versor = meas.pull / sqrt(meas.pull.sq)
} else {
  meas.pull.versor = meas.pull * 0
}
if (meas.pull.sq > dof.pull && meas.pull.sq > 1) {
  ##--- add to measurements correlation a matrix proportional to pull vector outer product
  meas.corr.pull = drop(t(meas.pull.versor) %*% meas.corr %*% meas.pull.versor)
  meas.pull.sf.corr = (meas.pull.sq - 1) * meas.corr.pull
  meas.corr.orin.full = meas.corr + (meas.pull.versor %o% meas.pull.versor) * meas.pull.sf.corr
} else {
  meas.corr.orin.full = meas.corr
}

##--- Orin Dahl S-factor inflated measurement covariance
meas.cov.orin.full = meas.corr.orin.full * (meas.err %o% meas.err)
orin.full = get.fit.errors(meas.cov.orin.full, (meas.val | TRUE))
##--- alu S-factor inflated measurement covariance
meas.cov.alu.full = meas.cov * (meas.sfact.alu.full %o% meas.sfact.alu.full)
alu.full = get.fit.errors(meas.cov.alu.full, (meas.val | TRUE))

##
## S-factors, grouping by relation to fitted quantity
##
meas.sfact.alu.fq = meas.val * 0 + 1
quant.sfact.alu.fq = quant.val * 0 + 1
meas.keep.fq = meas.val & FALSE
meas.corr.orin.fq = meas.corr
for (mt.name in quant.names) {
  ##--- selection of measurements connected to the specified fitted quantity
  mm.list = (delta[, mt.name] == 1)
  mm.num = sum(mm.list)
  if (mm.num == 0) next

  ##--- compute max error to retain measurement for S-factor calculation
  error.max = 3*sqrt(mm.num)*quant.err[mt.name]
  ##--- keep only measurement with not too large errors as in PDG S-factor calculation
  mm.list.keep = mm.list & (meas.err < error.max)
  meas.keep.fq = meas.keep.fq | mm.list.keep

  ##
  ## print measurement with too large errors
  ##
  save = options()
  options(digits=4)
  if (any(xor(mm.list, mm.list.keep))) {
    cat("\nS-factor calculation, exclude because error > 3*sqrt(N)*av_err=", error.max, "\n")
    excl = meas.err[xor(mm.list, mm.list.keep)]
    cat(paste(names(excl), "error=", format(excl, digits=4), "nsigma=", format(excl/error.max*3, digits=4), sep=" "), sep="\n")
  }
  if (FALSE && any(mm.list.keep)) {
    cat("\nS-factor calculation, included measurements, 3*sqrt(N)*av_err=", error.max, "\n")
    excl = meas.err[mm.list.keep]
    cat(paste(names(excl), "error=", format(excl, digits=4), "nsigma=", format(excl/error.max*3, digits=4), sep=" "), sep="\n")
  }
  options(save)

  mm.list.all = mm.list
  mm.list = mm.list.keep
  if (sum(mm.list) == 0) next

  ##--- get proper section of full pull covariance
  pull.mm.cov = pull.cov[mm.list, mm.list, drop=FALSE]
  ##--- get just diagonal part of pull matrix (equivalent to old implementation)
  ## pull.mm.cov = diag.m(diag(pull.mm.cov))

  pull.mm.invcov.sqrt = alu.matr.inv.sqrt.symm.semipos(pull.mm.cov)
  dof.mm = attr(pull.mm.invcov.sqrt, "pos.eigen.num")
  ##--- go to next measurement group if the pull matrix is entirely singular
  if (dof.mm == 0) next

  pull = drop(pull.mm.invcov.sqrt %*% (meas.val[mm.list] - (delta %*% quant.val)[mm.list]))
  ##--- force to zero very small pulls
  pull[abs(pull) <= tol] = 0
  pull.sq = sum(pull^2)
  
  ##--- one single S-factor per measurement group
  ## sfact = sqrt(pull.sq/dof.mm)
  ##--- compute average chisq term
  sfact = sqrt(sum(meas.pull.dof[mm.list]^2 * meas.dof.eff[mm.list])/sum(meas.dof.eff[mm.list]))

  if (sfact < 1) sfact = 1
  quant.sfact.alu.fq[mt.name] = sfact
  meas.sfact.alu.fq[mm.list.all] = sfact

  ##--- use pull direction as in Orin Dahl S-factors prescription for grouping by statistical correlation
  if (pull.sq != 0) {
    pull.versor = pull / sqrt(pull.sq)
  } else {
    pull.versor = pull * 0
  }
  if (pull.sq > dof.mm && pull.sq > 1) {
    ##--- add to measurements correlation a matrix proportional to pull vector outer product
    meas.corr.pull = drop(t(pull.versor) %*% meas.corr[mm.list, mm.list] %*% pull.versor)
    pull.sf.corr = (pull.sq - 1) * meas.corr.pull
    meas.corr.orin.fq[mm.list, mm.list] = meas.corr[mm.list, mm.list, drop=FALSE] + (pull.versor %o% pull.versor) * pull.sf.corr
  }
}
##--- Orin Dahl S-factor inflated measurement covariance, alternative grouping by fitted quantity
meas.cov.orin.fq = meas.corr.orin.fq * (meas.err %o% meas.err)
orin.fq = get.fit.errors(meas.cov.orin.fq, meas.keep.fq)
##--- alu S-factor inflated measurement covariance, grouping by fitted quantity
meas.cov.alu.fq = meas.cov * (meas.sfact.alu.fq %o% meas.sfact.alu.fq)
alu.fq = get.fit.errors(meas.cov.alu.fq, meas.keep.fq)

##
## S-factors, grouping by statistical correlation (Orin Dahl)
##
## - group measurements
##   - first among statistically correlated ones
##   - second by their connection to a fitted quantity
## - compute pulls for each measurement group
##

##--- fix diagonal of statistical correlation matrix
diag(meas.corr.stat) = 1
meas.correlated = colSums(meas.corr.stat != 0) > 1
meas.correlated.which = which(meas.correlated)
meas.correlated.list = list()

##--- group statistically correlated measurements
repeat {
  if (length(meas.correlated.which) == 0) break
  mm = meas.correlated.which[1]
  mm.list = which(meas.corr.stat[mm, ] != 0)
  meas.correlated.list = c(meas.correlated.list, list(mm.list))
  meas.correlated.which = setdiff(meas.correlated.which, mm.list)
}

##--- group statistically uncorrelated measurements per fitted quantity
meas.uncorrelated.list = list()
for (mt.name in quant.names) {
  ##--- selection of measurements of type mt.name but not stat correlated
  mm.list = (delta[, mt.name] == 1) & !meas.correlated
  mm.num = sum(mm.list)
  if (mm.num == 0) next

  ##--- compute max error to retain measurement for S-factor calculation
  error.max = 3*sqrt(mm.num)*quant.err[mt.name]
  ##--- keep only measurement with not too large errors as in PDG S-factor calculation
  mm.list.all = mm.list
  mm.list = mm.list & (meas.err < error.max)
  if (sum(mm.list) == 0) next

  mm.list.which = which(mm.list)
  attr(mm.list.which, "fitted.quantity") = mt.name
  attr(mm.list.which, "not.kept") = which(mm.list.all & ! mm.list)
  meas.uncorrelated.list[[mt.name]] = mm.list.which
}

##--- compute S-factors
meas.corr = meas.cov / (meas.err %o% meas.err)
meas.corr.orin.sc = meas.corr
meas.sfact.alu.sc = meas.val * 0 + 1
meas.sfact.orin.sc = meas.val * 0 + 1
meas.keep.sc = meas.val & FALSE
for (mm.list in c(meas.correlated.list, meas.uncorrelated.list)) {
  ##--- assemble list of kept measurements for S-factors determination
  meas.keep.sc[mm.list] = TRUE

  ## pull.mm.cov = meas.cov[mm.list, mm.list, drop=FALSE] - (delta %*% quant.cov %*% t(delta))[mm.list, mm.list, drop=FALSE]
  pull.mm.cov = pull.cov[mm.list, mm.list, drop=FALSE]

  pull.mm.invcov.sqrt = alu.matr.inv.sqrt.symm.semipos(pull.mm.cov)
  dof.mm = attr(pull.mm.invcov.sqrt, "pos.eigen.num")

  ##--- go to next measurement group if the pull matrix is entirely singular
  if (dof.mm == 0) next

  pull = drop(pull.mm.invcov.sqrt %*% (meas.val[mm.list] - (delta %*% quant.val)[mm.list]))
  ##--- force to zero very small pulls
  pull[abs(pull) <= tol] = 0
  pull.sq = sum(pull^2)

  mt.name = attr(mm.list, "fitted.quantity")
  if (!is.null(mt.name)) {
    ##--- measurements connected to a single fit quantity
    ## sfact = sqrt(pull.sq/dof.mm) 
    ##--- compute average chisq term
    sfact = sqrt(sum(meas.pull.dof[mm.list]^2 * meas.dof.eff[mm.list])/sum(meas.dof.eff[mm.list]))

    if (sfact < 1) sfact = 1
    meas.sfact.alu.sc[c(mm.list, attr(mm.list, "not.kept"))] = sfact
    meas.sfact.orin.sc[c(mm.list, attr(mm.list, "not.kept"))] = sfact
  } else {
    ##--- measurements grouped by statistical correlation

    ##--- Orin Dahl S-factors prescription, grouping by statistical correlation
    if (pull.sq != 0) {
      pull.versor = pull / sqrt(pull.sq)
    } else {
      pull.versor = pull * 0
    }
    if (pull.sq > dof.mm && pull.sq > 1) {
      ##--- add to measurements correlation a matrix proportional to pull vector outer product
      meas.corr.pull = drop(t(pull.versor) %*% meas.corr[mm.list, mm.list] %*% pull.versor)
      pull.sf.corr = (pull.sq - 1) * meas.corr.pull
      meas.corr.orin.sc[mm.list, mm.list] = meas.corr[mm.list, mm.list, drop=FALSE] + (pull.versor %o% pull.versor) * pull.sf.corr
    }

    ##--- A.Lusiani S-factors prescription #2, grouping by statistical correlation
    ## pull.dof = abs(pull) * sapply(attr(pull.mm.invcov.sqrt, "dof.eff"), function(x) ifelse(x!=0, 1/sqrt(x), 0))
    pull.dof = meas.pull.dof[mm.list]
    meas.sfact.alu.sc[mm.list] = ifelse(pull.dof > 1, pull.dof, 1)
  }
}
##--- Orin Dahl S-factor inflated measurement covariance, grouping by statistical correlation
meas.cov.orin.sc = meas.corr.orin.sc * (meas.err %o% meas.err)
meas.cov.orin.sc = meas.cov.orin.sc * (meas.sfact.orin.sc %o% meas.sfact.orin.sc)
orin.sc = get.fit.errors(meas.cov.orin.sc, meas.keep.sc)
##--- alu S-factor inflated measurement covariance, grouping by statistical correlation
meas.cov.alu.sc = meas.cov * (meas.sfact.alu.sc %o% meas.sfact.alu.sc)
alu.sc = get.fit.errors(meas.cov.alu.sc, meas.keep.sc)

if (any(meas.keep.fq)) {
  ##--- chisq/dof for kept measurements grouped by fitted quantity
  chisq.keep.fq = drop(
    t((meas.val - delta %*% quant.val)[meas.keep.fq])
    %*% meas.invcov[meas.keep.fq, meas.keep.fq]
    %*% (meas.val - delta %*% quant.val)[meas.keep.fq])
  ##--- assume no quantity left without a kept measurement
  dof.keep.fq = sum(meas.keep.fq) - (quant.num - constr.num)
}

if (any(meas.keep.sc)) {
  ##--- chisq/dof for kept measurements grouped by statistical correlation
  chisq.keep.sc = drop(
    t((meas.val - delta %*% quant.val)[meas.keep.sc])
    %*% meas.invcov[meas.keep.sc, meas.keep.sc]
    %*% (meas.val - delta %*% quant.val)[meas.keep.sc])
  ##--- assume no quantity left without a kept measurement
  dof.keep.sc = sum(meas.keep.sc) - (quant.num - constr.num)
}

quant.sfact = quant.sfact.alu.fq
quant.sf.err = alu.full$quant.err
quant.sf.cov = alu.full$quant.cov
quant.sf.corr = alu.full$quant.corr
quant.sf.sfact = alu.full$quant.sfact

cat("\n")
cat("##\n")
cat("## S-factors accounting for larger than expected chi-square\n")
cat("##\n")

get.chisq.dof.pchisq = function(fit.errors) {
  c("chisq"     = fit.errors$chisq,
    "dof"       = fit.errors$dof,
    "chisq/dof" = fit.errors$chisq/dof,
    "CL"        = pchisq(fit.errors$chisq, fit.errors$dof, lower.tail=FALSE))
}

get.chisq.dof.pchisq.keep = function(fit.errors) {
  if (any(fit.errors$meas.keep)) {
    c("chisq"     = fit.errors$chisq.keep,
      "dof"       = fit.errors$dof.keep,
      "chisq/dof" = fit.errors$chisq.keep/fit.errors$dof.keep,
      "CL"        = pchisq(fit.errors$chisq.keep, fit.errors$dof.keep, lower.tail=FALSE))
  } else {
    NULL
  }
}

out = rbind(
  "fit" = c(
    chisq=chisq, dof=dof,
    "chisq/dof"=chisq/dof,
    CL=pchisq(chisq, dof, lower.tail=FALSE))
  , "  orin.full" = get.chisq.dof.pchisq(orin.full)
  , "  alu.full" = get.chisq.dof.pchisq(alu.full)
  , "  orin.sc" = get.chisq.dof.pchisq(orin.sc)
  , "  alu.sc" = get.chisq.dof.pchisq(alu.sc)
  , "  orin.fq" = get.chisq.dof.pchisq(orin.fq)
  , "  alu.fq" = get.chisq.dof.pchisq(alu.fq)
  , "keep.sc" =
  if (any(meas.keep.sc)) {
    c(chisq=chisq.keep.sc, dof=dof.keep.sc,
      "chisq/dof"=chisq.keep.sc/dof.keep.sc,
      CL=pchisq(chisq.keep.sc, dof.keep.sc, lower.tail=FALSE))
  } else NULL
  , "  alu.sc" = get.chisq.dof.pchisq.keep(alu.sc)
  , "  orin.sc" = get.chisq.dof.pchisq.keep(orin.sc)
  , "keep.fq" =
  if (any(meas.keep.fq)) {
    c("chisq"     = chisq.keep.fq,
      "dof"       = dof.keep.fq,
      "chisq/dof" = chisq.keep.fq/dof.keep.fq,
      "CL"        = pchisq(chisq.keep.fq, dof.keep.fq, lower.tail=FALSE))
  } else NULL
  , "  alu.fq" = get.chisq.dof.pchisq.keep(alu.fq)
  , "  orin.fq" = get.chisq.dof.pchisq.keep(orin.fq)
  )
show(out)

cat("Averaged quantities: value, error, error with S-factor, S-factor\n") 
show(rbind("value"=quant.val,
           "error"=quant.err,
           "  orin.full"=orin.full$quant.err,
           "  alu.full"=alu.full$quant.err,
           "  orin.sc"=orin.sc$quant.err,
           "  alu.sc"=alu.sc$quant.err,
           "  orin.fq"=orin.fq$quant.err,
           "  alu.fq"=alu.fq$quant.err,
           "S-factor orin.full"=orin.full$quant.sfact,
           "S-factor alu.full"=alu.full$quant.sfact,
           "S-factor orin.sc"=orin.sc$quant.sfact,
           "S-factor alu.sc"=alu.sc$quant.sfact,
           "S-factor orin.fq"=orin.fq$quant.sfact,
           "S-factor alu.fq"=alu.fq$quant.sfact,
           "  quant"=quant.sfact.alu.fq
           ))

if (FALSE && quant.num > 1) {
  ## cat("correlation\n") 
  ## show(quant.corr)
  cat("correlation, S-factor inflated\n") 
  show(orin$quant.corr)
}

##--- save data and results
rc = save(file=file.name.data,
  measurements, combination, delta,
  chisq, dof,
  meas.val,   meas.err,   meas.cov,  meas.cov.stat, meas.cov.syst, meas.corr, meas.sfact.cards,
  quant.val,  quant.err,    quant.cov,    quant.corr,    quant.sfact,
              quant.sf.err, quant.sf.cov, quant.sf.corr, quant.sf.sfact,
  alu.full, orin.full,
  alu.fq, alu.sc, orin.fq, orin.sc,
  constr.m, constr.v)

cat("\n")
cat(paste("file", file.name.data, "produced\n"))
cat("\n")
cat("## end\n")

options(options.save)
##+++ } ##--- end function alucomb

args <- commandArgs(TRUE)
if (length(args) > 0 && exists("alucomb")) alucomb(file = args[1]) 
