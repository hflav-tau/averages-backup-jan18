#!/usr/bin/env Rscript

##
## alucomb.r
##
## Copyright Alberto Lusiani 2010. All rights reserved.
## an open source license will be set once tested and working
##
## - averages measurements in a way that is mostly compatible with Combos
## - can read Combos input files
## - can average multiple quantities related to multiple statistically correlated measurements
##

library(methods)
source("../../../Common/bin/alu-utils2.r")

## ////////////////////////////////////////////////////////////////////////////
## definitions

method="solnp"
method="alucomb2"

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

file.name.data = gsub("[.][^.]*$", ".rdata", file.name)

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
quant.names = names(combination$quantities)

##--- check duplicate linear constraints
if (length(combination$constr.lin.comb) > 1) {
  dupl.constr = NULL
  for (i in seq(1, length(combination$constr.lin.comb))) {
    for (j in seq(i+1, length=length(combination$constr.lin.comb)-i)) {
      if (length(setdiff(names(combination$constr.lin.comb[[i]]), names(combination$constr.lin.comb[[j]]))) == 0) {
        qn = names(combination$constr.lin.comb[[i]])
        dupl.constr = rbind(dupl.constr,
          matrix(c(combination$constr.lin.val[[i]], unlist(combination$constr.lin.comb[[i]][qn])), nrow=1,
                 dimnames=list(names(combination$constr.lin.comb[i]), c("val", qn))),
          matrix(c(combination$constr.lin.val[[j]], unlist(combination$constr.lin.comb[[j]][qn])), nrow=1,
                 dimnames=list(names(combination$constr.lin.comb[j]), c("val", qn))))
      }
    }
  }
  if (length(dupl.constr) > 0) {
    cat("warning: duplicated constraints, listed in row pairs\n")
    show(dupl.constr)
  }
}

##--- retain only constraints whose terms are all included in the fitted quantities
if (length(combination$constr.lin.comb) > 0) {
  constr.select = sapply(combination$constr.lin.comb, function(x) all(names(x) %in% quant.names))
  if (any(!constr.select)) {
    cat("\nThe following constraints are dropped:\n")
    mapply(function(comb, val, val.name) {
      tmp = val
      names(tmp) = val.name
      show(c(comb, tmp))
    },
    combination$constr.lin.comb[!constr.select],
    combination$constr.lin.val[!constr.select],
    names(combination$constr.lin.val[!constr.select]))
    cat("\nThe following measurement types are missing:\n")
    show(unique(unlist(lapply(combination$constr.lin.comb, function(x) setdiff(names(x), quant.names)))))
  }
  combination$constr.lin.comb = combination$constr.lin.comb[constr.select]
  combination$constr.lin.val = combination$constr.lin.val[constr.select]
}

##--- quantity measured per measurement
meas.quantities = sapply(measurements, function(x) x$quant)

##--- discard measurements that are not associated to a declared fitted quantity
meas.included.list = meas.quantities %in% quant.names
names(meas.included.list) = names(meas.quantities)[meas.included.list]
meas.names.discarded =  meas.names[!meas.included.list]
if (length(meas.names.discarded) >0) {
  cat("\nwarning: the following measurements are discarded:\n")
  cat(paste("  ", meas.names.discarded, collapse="\n"), "\n");
}

##--- quantities involved in constraints
constr.quantities = unique(unlist(lapply(combination$constr.lin.comb, function(x) names(x)), use.names=FALSE))

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
meas.num = length(measurements)
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
quant.sfact.list = lapply(combination$quantities.options, function(el) { unname(el["scale"]) })
quant.sfact.list = quant.sfact.list[!is.na(quant.sfact.list)]
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
          ##--- log updates to measurements values and uncertainties
          rc = cat(format(measurements[[mn]]$tag,width=30),
            format(param.orig,width=15),
            format(c(measurements[[mn]]$params[[param.orig]]["value"],
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
meas.val = sapply(measurements, function(x) {unname(x$value)})
meas.stat = sapply(measurements, function(x) {unname(x$stat)})
meas.syst = sapply(measurements, function(x) {unname(x$syst)})
meas.err = sqrt(meas.stat^2 + meas.syst^2)

##--- unshifted values
meas.val.orig = sapply(measurements, function(x) {unname(x$value.orig)})
meas.syst.orig = sapply(measurements, function(x) {unname(x$syst.orig)})

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
## set seed values for quantities in quant.seed.val
## - this is needed to linearize the constraint equations
## - where measurements are not available, the seed values must be set in the MEASUREMENT cards in the COMBINE block
## - elsewhere seed values are computed using the measurements without using the constraints
##

##--- get quantities that have at least one measurement
quant.measured.bool = quant.names %in% meas.quantities
##--- assemble delta matrix for just those quantities
delta.measured = delta[,quant.measured.bool]
##--- fit for quantities without constraints
quant.seed.val = quant.val
quant.seed.val[quant.measured.bool] =
  solve(t(delta.measured) %*% meas.invcov %*% delta.measured) %*%
  (t(delta.measured) %*% meas.invcov %*% meas.val)

##--- get seed values for quantities without measurements
quant.cards.seed.val = unlist(lapply(combination$quantities, function(el) { unname(el["seed"]) }))
if (is.null(quant.cards.seed.val)) {
  quant.cards.seed.val = numeric(0)
}
quant.cards.seed.val = quant.cards.seed.val[!is.na(quant.cards.seed.val)]

##--- set seed values for quantities that have no measurements
seed.needed = setdiff(quant.names, quant.names[quant.measured.bool])
seed.needed.notincards = setdiff(seed.needed, names(quant.cards.seed.val))
if (length(seed.needed.notincards)>0) {
  cat("\nwarning: no seed value provided for the following quantities\n")
  show(seed.needed.notincards)
  cat("warning: please set them in the cards, (default values of zero used)\n")
}
quant.seed.val[seed.needed] = quant.cards.seed.val[seed.needed]

## 
## prepare constraints for minimization constraints
## - transform linear constraints combinations into expressions
## - transform non-linear constraint strings into expressions
## - join linear and non-linear expressions in a single list
## - join linear and non-linear constraint constants in a single list
## - save all lists in "combination" list
##

##--- get expressions corresponding to linear constraints combinations
constr.lin.expr = mapply(function(comb) {
  constr = paste(mapply(function(x,y) paste(x, "*", y, sep=""), comb, names(comb)), collapse = "+")
  parse(text=constr)
}, combination$constr.lin.comb)

##--- get expressions corresponding to non-linear constraints combinations
constr.nl.expr = parse(text=combination$constr.nl.expr)
names(constr.nl.expr) = names(combination$constr.nl.expr)

##--- join linear and non-linear constrainst expressions
combination$constr.all.expr = c(constr.lin.expr, constr.nl.expr)

##--- join linear and non-linear constraint values
combination$constr.all.val = unlist(c(combination$constr.lin.val, combination$constr.nl.val))

##--- flag which constraints are non-linear
combination$constr.all.nl = c(rep(FALSE, length(combination$constr.lin.val)), rep(TRUE, length(combination$constr.nl.val)))

##--- print constraint equations
if (length(combination$constr.all.val) > 0) {
  cat("\n## Constraint equations begin\n\n")
  cat(paste(combination$constr.all.val, as.character(combination$constr.all.expr), sep=" = "), sep="\n")
  cat("\n## Constraint equations end\n")
}

##
## ////////////////////////////////////////////////////////////////////////////
##

if (method == "alucomb") {

##
## obtain constraint equations
##
constr.num = length(combination$constr.lin.comb)
constr.names = names(combination$constr.lin.comb)
constr.m =  do.call(rbind, lapply(combination$constr.lin.comb, function(x) {tmp = quant.val; tmp[names(x)] = x; tmp}))
constr.v = unlist(combination$constr.lin.val)

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

  if (FALSE) {
    cat("\n## Constraint equations begin\n\n")
    tmp = mapply(function(name, val, comb) {names(val) = name; show(c(val, unlist(comb)))},
      names(constr.v), constr.v/quant.invcov.order*constr.m.order,
      apply(constr.m/quant.invcov.order*constr.m.order, 1, function(x) list(x[x!=0])))
    cat("\n## Constraint equations end\n")
  }
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

##--- save data and results
rc = save(file=file.name.data,
  measurements, combination, delta,
  chisq, dof,
  meas.val,   meas.err,   meas.cov,  meas.cov.stat, meas.cov.syst, meas.corr, meas.sfact.cards,
  quant.val,  quant.err,  quant.cov, quant.corr,
  constr.m, constr.v)

} # end if method alucomb

##
## ////////////////////////////////////////////////////////////////////////////
##

if (method == "alucomb2") {

##--- init quant.val
quant.val = quant.seed.val

first.iteration = TRUE
constr.num = length(combination$constr.all.expr)
constr.names = names(combination$constr.all.expr)
constr.nl = combination$constr.all.nl

##--- derivative and gradient of constraint equation expressions
constr.expr = lapply(combination$constr.all.expr, function(x) deriv(x, all.vars(x)))

##--- print linearized constraints
print.linearized.constraints = function(val, comb) {
  tmp = mapply(function(val, comb) {
    cat(val, "=", paste(mapply(function(name, val) {
      paste(val, "*", name, sep="")
    }, names(comb), comb), collapse=" + "), "\n")
  }, val, comb)
}

repeat {
  ##
  ## linearize just the non-linear constraint equations
  ## f_j(x_i) = c_j becomes
  ## ... f(x^0_i) + f_j'_i*(x_i - x^0_i) = c_j
  ## ... f_j'_i*x_i = c_j - f(x^0_i) + f_j'_i*x^0_i
  ##
  constr.expr.val = lapply(constr.expr[constr.nl], function(x) eval(x, as.list(quant.val)))
  constr.grad.comb = lapply(constr.expr.val, function(x) drop(attr(x, "gradient")))
  constr.grad.val = combination$constr.all.val[constr.nl] - as.vector(unlist(constr.expr.val))
  if (sum(constr.nl) > 0) {
    constr.grad.val = constr.grad.val + sapply(constr.grad.comb, function(x) drop(x %*% quant.val[names(x)]))
  }

  ##--- join linear constraints
  constr.grad.comb = c(combination$constr.lin.comb, constr.grad.comb)
  constr.grad.val = c(unlist(combination$constr.lin.val), constr.grad.val)
  
  ##--- obtain constraint equations
  constr.m = do.call(rbind, lapply(constr.grad.comb, function(x) {tmp = quant.val*0; tmp[names(x)] = x; tmp}))
  constr.v = constr.grad.val
  
  if (constr.num > 0) {
    if (first.iteration) {
      ##--- determine the typical size of quant.invcov elements
      sv = svd(quant.invcov)$d
      sv.central = round(quant.num*1/3):round(quant.num*2/3)
      sv.log.mean = mean(log(sv[sv.central]))
      quant.invcov.order = 10^round(sv.log.mean/log(10))
      
      ##--- determine the typical size of the constraint equation terms
      constr.m.order = 10^round(log(mean(abs(constr.m[constr.m!=0])))/log(10))

      if (TRUE) {
        cat("\n## Begin of linearized constraint equations (1st iteration)\n\n")
        print.linearized.constraints(constr.grad.val[constr.nl], constr.grad.comb[constr.nl])
        cat("\n## End of linearized constraint equations (1st iteration)\n")
      }

      if (FALSE) {
        cat("\n## Begin of linearized constraint equations (1st iteration)\n\n")
        tmp = mapply(function(name, val, comb) {names(val) = name; show(c(val, unlist(comb)))},
          names(constr.v), constr.v,
          apply(constr.m, 1, function(x) list(x[x!=0])))
        cat("\n## End of linearized constraint equations (1st iteration)\n")
      }
      
      cat("\n## Begin of constraint percent change summaries\n\n")
    }
    ##--- to avoid computationally singular matrix, apply proper factor to constraint equations
    constr.m = constr.m * quant.invcov.order/constr.m.order
    constr.v = constr.v * quant.invcov.order/constr.m.order
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
  
  if (length(constr.grad.comb[constr.nl]) == 0) {
    ##--- there is no non-linear constraint whose linearization must iteratively converge
    break
  }

  if (!first.iteration) {
    ##--- compute sum of squares of differences between previous and current iteration linearized constraints
    constr.diff = mapply(function(c2, v2, c1, v1) {
      x2 = c(c2, v2)
      x1 = c(c1, v1)
      norm = pmax(x1,x2,x2+x1)
      norm = mean(norm)
      ifelse(norm==0, 0, sum(((x2-x1)/norm)^2))
    }, constr.grad.comb[constr.nl], constr.grad.val[constr.nl], constr.grad.prev.comb[constr.nl], constr.grad.prev.val[constr.nl])
    show(constr.diff)
    ##--- end if the average percent change of constraint equation coefficients is small enough
    if (sum(constr.diff) / constr.num < 1e-10) break
  }

  ##--- save before calculation of next updated values
  constr.grad.prev.comb = constr.grad.comb
  constr.grad.prev.val = constr.grad.val
  first.iteration = FALSE
}
if (constr.num > 0) {
  cat("\n## End of constraint percent change summaries\n")
}
rm(first.iteration)
rm(constr.diff)

if (TRUE) {
  cat("\n## Begin of linearized constraint equations (after convergence)\n\n")
  print.linearized.constraints(constr.grad.val[constr.nl], constr.grad.comb[constr.nl])
  cat("\n## End of linearized constraint equations (after convergence)\n")
}

##
## compute errors and chi square
##
quant.cov = solve.m[1:quant.num,1:meas.num, drop=FALSE] %*% meas.cov %*% t(solve.m[1:quant.num,1:meas.num, drop=FALSE])
quant.cov = (quant.cov + t(quant.cov))/2
quant.err = sqrt(diag(quant.cov))
quant.corr = quant.cov / (quant.err %o% quant.err)

chisq = drop(t(meas.val - delta %*% quant.val) %*% meas.invcov %*% (meas.val - delta %*% quant.val))
dof = meas.num - quant.num + constr.num

cat("\n##\n")
cat("## alucomb2 solution, chisq/d.o.f. = ",chisq, "/", dof, ", CL = ", (1-pchisq(chisq, df=dof)), "\n",sep="")
cat("##\n")
show(rbind(value=quant.val[1:quant.num], error=quant.err[1:quant.num]))
if (FALSE && quant.num > 1) {
  cat("correlation\n")
  show(quant.corr[1:quant.num, 1:quant.num])
}
cat("## end\n")

##--- save data and results
rc = save(file=file.name.data,
  measurements, combination, delta,
  chisq, dof,
  meas.val,   meas.err,   meas.cov,  meas.cov.stat, meas.cov.syst, meas.corr, meas.sfact.cards,
  quant.val,  quant.err,  quant.cov, quant.corr,
  constr.m, constr.v)

} # end if method alucomb2

##
## ////////////////////////////////////////////////////////////////////////////
##

if (method == "solnp") {

library(truncnorm)
library(Rsolnp)

##--- degrees of freedom
dof = meas.num - quant.num + length(combination$constr.all.val)

##--- function returns vector of non-linear constraint functions
constr.fun <- function (qv)  {
  env.list = as.list(qv)
  rc = unlist(lapply(combination$constr.all.expr, function(expr) {eval(expr, env.list)}))
  return(rc)
}

##--- function returns 1/2 of chi square
half.chisq.fun <- function(qv) {
  drop(t(meas.val - delta %*% qv) %*% meas.invcov %*% (meas.val - delta %*% qv))/2
}

##--- minimize half chi square to get the correct covariance for the fitted parameters
cat("\nBegin of solnp minimization\n")
rc.solnp <- solnp(quant.seed.val, fun = half.chisq.fun, eqfun = constr.fun, eqB = combination$constr.all.val)
cat("\nEnd of solnp minimization\n")

chisq = 2*tail(rc.solnp$values, 1)
quant.val = rc.solnp$pars
quant.cov = solve(rc.solnp$hessian)
rownames(quant.cov) = quant.names
colnames(quant.cov) = quant.names
quant.err = sqrt(diag(quant.cov))
quant.corr = quant.cov / (quant.err %o% quant.err)

cat("\n")
cat("##\n")
cat("## solnp solution, chisq/d.o.f. = ", chisq, "/", dof, ", CL = ", (1-pchisq(chisq, df=dof)), "\n",sep="")
cat("##\n")
show(rbind(value=quant.val, error=quant.err))
if (FALSE && quant.num > 1) {
  cat("correlation\n")
  show(quant.corr)
}
cat("## end\n")

##--- cleanup
if (FALSE) {
  rm(constr.nl.expr, rc.solnp)
}

##--- save data and results
rc = save(file=file.name.data,
  measurements, combination, delta,
  chisq, dof,
  meas.val,   meas.err,   meas.cov,  meas.cov.stat, meas.cov.syst, meas.corr, meas.sfact.cards,
  quant.val,  quant.err,  quant.cov, quant.corr)

} # end if method solnp

##
## ////////////////////////////////////////////////////////////////////////////
##

cat("\n")
cat(paste("file", file.name.data, "produced\n"))
cat("\n")
cat("## end\n")

options(options.save)
##+++ } ##--- end function alucomb

args <- commandArgs(TRUE)
if (length(args) > 0 && exists("alucomb")) alucomb(file = args[1]) 
