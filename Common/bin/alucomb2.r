#!/usr/bin/env Rscript

##
## alucomb2.r
##
## Copyright Alberto Lusiani 2010. All rights reserved.
## an open source license will be set once tested and working
## usage of this code is allowed provided it is properly referenced
##
## - averages measurements by minimizing the chi square (Gaussian errors assumed)
## - reads combos-like cards with:
##   - measurements and their statistical and systematic errors
##   - systematic error terms due to external parameters
##   - statistical correlations among different measurements
##   - linear and non-linear constraints on the fitted quantities
## - non-linear constraints are linearized and the chi square is minimized analytically
##   by solving a set of linear equations; the constraint linearization is iteratively
##   optimized by repeating the fit until convergence is reached
##

source("../../../Common/bin/alucomb2-utils.r")

## ////////////////////////////////////////////////////////////////////////////
## definitions

method="solnp"
method="alabama"
method="alucomb2"

## ////////////////////////////////////////////////////////////////////////////
## functions

##--- utility function to only have an empty output line when needed
print.empty.line.first.time = function() {
  if (flag.empty.line.first.time) {
    cat("\n")
    flag.empty.line.first.time <<- FALSE
  }
}

##--- print parameters
alucomb2.print.params = function(params) {
  if (length(params) > 0) {
    cat("\nPARAMETERS\n")
    mapply(function(label, value) {
      cat(" ", format(label, width=16), sprintf("%-12g", value[1]))
      if (value[2] == -value[3]) {
        cat(" +-", value[2], "\n", sep="")
      } else {
        cat(" +", sprintf("%-12g", value[2]), " -", value[3], "\n", sep="")
      }
    }, names(params), params)
  }
}

##--- print systematic terms, the original input string and the resulting absolute value
alucomb2.print.meas.syst.terms = function(syst.label, syst.terms, mask) {
  syst.terms.input = attr(syst.terms, "input")[mask]
  syst.terms = syst.terms[mask]
  if (length(syst.terms) > 0) {
    cat("\n", syst.label, "\n", sep="")
    mapply(function(label, val.abs, val.input) {
      cat("  ", format(val.input, width=16), " ", format(label, width=16), sep="")
      if (val.input != val.abs) cat(" # ", val.abs, sep="")
      cat("\n")
    }, names(syst.terms), sprintf("%+g", syst.terms), syst.terms.input)
  }
}

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

alucomb = function(file.name = "") {

##--- set very large line width to print even large amount of averaged quantities on single line
options.save = options()
## options(width=10000)

file.name.data = gsub("[.][^.]*$", ".rdata", file.name)

rc = alucomb.read(file.name)
measurements = rc$measurements[sort(names(rc$measurements))]
combination = rc$combination
rm(rc)

##
## build list of all measurements mentioned in the COMBINE section
## (right now we only handle COMBINE * * *)
##

##--- quantities to be averaged
quant.names = combination$combine

cat("\n##\n")
cat("## averaging the following quantities\n")
cat("##\n")
cat("COMBINE\n")
quant.fmt = format(quant.names)
maxlen = max(nchar(quant.fmt))
items.per.row = floor((79-2)/maxlen)
for (i.first in seq(1, length(quant.fmt), by=items.per.row)) {
  i.last = min(i.first + items.per.row - 1, length(quant.fmt) )
  cat("  ", paste(quant.fmt[i.first:i.last]), "\n", sep="")
}

flag.empty.line.first.time <<- TRUE
fields.descr.list = c("node", "descr")
for (quant.name in quant.names) {
  quant = combination$quantities[[quant.name]]
  fields = names(quant)
  fields.descr = intersect(fields.descr.list, fields)
  if (length(fields.descr) > 0 && any(quant[fields.descr] != "")) {
    print.empty.line.first.time()
    cat("QUANTITY", quant.name)
    for (field in fields.descr.list) {
      val = quant[[field]]
      if (!is.null(val) && val != "") cat("", field, val)
    }
    cat("\n")
  }
}

flag.empty.line.first.time <<- TRUE
for (quant.name in quant.names) {
  quant = combination$quantities[[quant.name]]
  fields = names(quant)
  fields.nodescr = setdiff(fields, fields.descr.list)
  if (length(fields.nodescr) > 0) {
    print.empty.line.first.time()
    cat("QUANTITY ", quant.name, " ", paste(fields.nodescr, quant[fields.nodescr]), "\n", sep="")
  }
}

cat("\n##\n")
cat("## quantities with no description\n")
cat("##\n")
for (quant.name in quant.names) {
  quant = combination$quantities[[quant.name]]
  if (is.null(quant$descr) || quant$descr == "") {
    cat("QUANTITY", quant.name, "node \"\" descr \"\"\n")
  }
}

##--- print measurements as read from the cards (short form)
if (FALSE) {
  cat("\n##\n")
  cat("## measurements as read from the cards\n")
  cat("##\n")
  rc = alu.rbind.print(
    rbind(value=meas.val,
          stat=meas.stat,
          syst=meas.syst,
          error=meas.err),
    num.columns=1)
}

##
## print measurements as read from the cards (all fields)
##
if (TRUE) {
  cat("\n##\n")
  cat("## measurements as read from the cards\n")
  cat("##\n")
  
  for (meas in measurements) {
    cat("\n")
    cat("# ", meas$quant, " (PDG node = ", combination$quantities[[meas$quant]]$node, ")\n", sep="")
    descr = combination$quantities[[meas$quant]]$descr
    if (!is.null(descr)) cat("# ", descr, "\n", sep="")
    cat("BEGIN MEASUREMENT ", paste(meas$tags, collapse=" "), "\n", sep="")

    cat("\n")
    cat("VALUE", meas$value)
    
    if (meas$stat.p == -meas$stat.n) {
      if (meas$stat != meas$stat.p) {
        stop(names(meas),": stat != stat.p, stat, stat.p, stat.n:", meas$stat, meas$stat.p, meas$stat.n)
      }
      cat(" +-", meas$stat, sep="")
    } else {
      cat(" (+", meas$stat.p, " -", -meas$stat.n, ")", sep="")
    }
    
    if (meas$syst.p == -meas$syst.n) {
      if (meas$syst != meas$syst.p) {
        stop(names(meas),": syst != syst.p, syst, syst.p, syst.n:", meas$syst, meas$syst.p, meas$syst.n)
      }
      cat(" +-", meas$syst, sep="")
    } else {
      cat(" (+", meas$syst.p, " -", -meas$syst.n, ")", sep="")
    }
    cat("\n")

    if (length(meas$corr.terms.stat) > 0) {
      cat("\n")
      mapply(function(label, value, input) {
        cat("STAT_CORR_WITH ", format(input, width=12), " ", gsub(".", " ", label, fixed=TRUE), sep="")
        if (value != input) cat(" # ", value, sep="")
        cat("\n")
      }, names(meas$corr.terms.stat), sprintf("%+-g", meas$corr.terms.stat), attr(meas$corr.terms.stat, "input"))
    }

    if (length(meas$corr.terms.tot) > 0) {
      cat("\n")
      mapply(function(label, value, input) {
        cat("TOT_CORR_WITH ", format(input, width=12), " ", gsub(".", " ", label, fixed=TRUE), sep="")
        if (value != input) cat(" # ", value, sep="")
        cat("\n")
      }, names(meas$corr.terms.tot), sprintf("%+-g", meas$corr.terms.tot), attr(meas$corr.terms.tot, "input"))
    }

    meas.name = paste(meas$tags, collapse=".")
    paper.name = paste(meas$tags[-2], collapse=".")
    syst.local.mask = substr(names(meas$syst.terms), 1, length(meas.name)) == meas.name
    syst.paper.mask = substr(names(meas$syst.terms), 1, length(meas.name)) == paper.name

    alucomb2.print.meas.syst.terms("SYSTEMATICS", meas$syst.terms, !(syst.local.mask | syst.paper.mask))
    alucomb2.print.meas.syst.terms("SYSTPAPER", meas$syst.terms, syst.paper.mask)
    alucomb2.print.meas.syst.terms("SYSTLOCAL", meas$syst.terms, syst.local.mask)
    
    alucomb2.print.params(meas$params)

    cat("\n")
    cat("END\n")
  }
}

##--- discard measurements according to COMBINATION cards
meas.names = names(measurements)
meas.drop.cards = combination$meas.drop.cards
if (is.null(meas.drop.cards)) meas.drop.cards = character(0)
meas.drop.cards.existing = meas.drop.cards %in% meas.names
if (any(!meas.drop.cards.existing)) {
  cat("\nwarning, cards require dropping non-existing measurements\n")
  cat(paste("  ", meas.drop.cards[!meas.drop.cards.existing], collapse="\n"), "\n")
}
if (any(meas.drop.cards.existing)) {
  cat("\n##\n")
  cat("## dropping the following measurements according to cards\n")
  cat("##\n")
  cat(paste("  ", meas.drop.cards[meas.drop.cards.existing], collapse="\n"), "\n")
  meas.names = setdiff(meas.names, meas.drop.cards[meas.drop.cards.existing])
}

##--- discard measurements that are not associated to a declared fitted quantity
meas.quantities = sapply(measurements[meas.names], function(x) x$quant)
meas.included.list = meas.quantities %in% quant.names
meas.names.discarded =  meas.names[!meas.included.list]
combination$meas.drop.unass = meas.names.discarded
if (length(meas.names.discarded) >0) {
  cat("\nwarning: measurements discarded because not in combined quantities:\n")
  cat(paste("  ", meas.names.discarded, collapse="\n"), "\n")
  meas.names = setdiff(meas.names, meas.names.discarded)
}

combination$measurements = meas.names
meas.used = names(measurements) %in% meas.names
meas.quantities = sapply(measurements[meas.used], function(x) x$quant)
meas.num = length(meas.names)

##--- get list of values, stat errors, syst errors
meas.val = sapply(measurements[meas.used], function(x) {unname(x$value)})
meas.stat = sapply(measurements[meas.used], function(x) {unname(x$stat)})
meas.syst = sapply(measurements[meas.used], function(x) {unname(x$syst)})
meas.err = sqrt(meas.stat^2 + meas.syst^2)

cat("\n##\n")
cat("## using the following updated global parameters\n")
cat("##\n")
alucomb2.print.params(combination$params)

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
    cat(dupl.constr, "\n")
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
      print(c(comb, tmp))
    },
    combination$constr.lin.comb[!constr.select],
    combination$constr.lin.val[!constr.select],
    names(combination$constr.lin.val[!constr.select]))
    cat("\nThe following quantities are missing:\n")
    print(unique(unlist(lapply(combination$constr.lin.comb, function(x) setdiff(names(x), quant.names)))))
  }
  combination$constr.lin.comb = combination$constr.lin.comb[constr.select]
  combination$constr.lin.val = combination$constr.lin.val[constr.select]
}

##--- quantities involved in constraints
constr.lin.quantities = unique(unlist(lapply(combination$constr.lin.comb, function(x) names(x)), use.names=FALSE))
constr.nl.quantities = all.vars(parse(text=combination$constr.nl.expr))
involved.quantities = unique(c(meas.quantities, constr.lin.quantities, constr.nl.quantities))

##--- discard fitted quantities that are not defined by measurements and constraints
quant.discarded = setdiff(quant.names, involved.quantities)
if (length(quant.discarded) > 0) {
  cat("\nwarning: following quantities discarded (not in measurements & constraints):\n")
  cat(paste("  ", quant.discarded, collapse="\n"), "\n");
  quant.names = setdiff(quant.names, quant.discarded)
  combination$combine.old = combination$combine
  combination$combine = quant.names
}

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
  print(slightly.larger)
}
if (length(larger) > 0) {
  cat("\nerror: syst. terms sum larger than total syst. error\n")
  colnames(larger) = c("total", "sum of terms")
  print(larger)
  stop("aborting")
}
if (length(not.matching) > 0) {
  cat("\nwarning: syst. terms do not match total syst. error, Combos requires that\n")
  colnames(not.matching) = c("total", "sum of terms")
  print(not.matching)
}

##
## shift measurements according to updated external parameter dependencies
## update systematic terms according to updated external parameter errors
##
flag.header.printed = FALSE
for (mn in meas.names) {
  ##--- save original values
  measurements[[mn]]$value.orig = measurements[[mn]]$value
  measurements[[mn]]$syst.orig = measurements[[mn]]$syst

  syst.upd.names = intersect(names(measurements[[mn]]$syst.terms), names(measurements[[mn]]$params))
  syst.upd.names = intersect(syst.upd.names, names(combination$params))
  if (length(syst.upd.names) == 0) next

  syst.term.orig = measurements[[mn]]$syst.terms[syst.upd.names]

  params.old       = sapply(measurements[[mn]]$params[syst.upd.names], function(x) unname(x["value"]))
  params.old.delta = sapply(measurements[[mn]]$params[syst.upd.names], function(x) unname(x["delta_pos"]))
  params.new       = sapply(combination$params[syst.upd.names], function(x) unname(x["value"]))
  params.new.delta = sapply(combination$params[syst.upd.names], function(x) unname(x["delta_pos"]))

  ##--- shift measurement values according to updated external parameters
  value.delta = (params.new - params.old) * (syst.term.orig / params.old.delta)

  ##--- differences of squares of uncertainties due to updated parameters
  syst.term.upd = syst.term.orig * (params.new.delta / params.old.delta)
  syst.term.deltasq = syst.term.upd^2 - syst.term.orig^2

  need.update = (value.delta != 0) || (syst.term.deltasq != 0)
  if (all(!need.update)) next

  ##--- update value
  measurements[[mn]]$value = measurements[[mn]]$value + sum(value.delta)
  ##--- update systematic error
  measurements[[mn]]$syst = sqrt(measurements[[mn]]$syst^2 + sum(syst.term.deltasq))
  ##--- update systematic terms
  measurements[[mn]]$syst.terms[syst.upd.names] = syst.term.upd
  
  if (TRUE) {
    params.old = params.old[need.update]
    params.old.delta = params.old.delta[need.update]
    params.new = params.new[need.update]
    params.new.delta = params.new.delta[need.update]
    syst.term.orig = syst.term.orig[need.update]
    syst.term.upd = syst.term.upd[need.update]
    
    ##--- log updates to measurements values and uncertainties
    if (!flag.header.printed) {
      cat("\n##\n")
      cat("## measurements updates due to updated parameters\n")
      cat("##\n")
      flag.header.printed = TRUE
    }
    cat("\n", mn, "\n", sep="")
    rc = print(rbind(
      old=c(value=measurements[[mn]]$value.orig, syst=measurements[[mn]]$syst.orig),
      new=c(value=measurements[[mn]]$value, syst=measurements[[mn]]$syst)))
    print(rbind(params.old, delta.old=params.old.delta, params.new, delta.new=params.new.delta,
                  syst.old=syst.term.orig, syst.new=syst.term.upd))
  }
}
suppressWarnings(rm(flag.header.printed, need.update, syst.term.orig, syst.term.upd, syst.term.deltasq))
suppressWarnings(rm(params.old, params.old.delta, params.new, params.new.delta))

##--- get list of updated values, stat errors, syst errors
meas.val = sapply(measurements[meas.used], function(x) {unname(x$value)})
meas.stat = sapply(measurements[meas.used], function(x) {unname(x$stat)})
meas.syst = sapply(measurements[meas.used], function(x) {unname(x$syst)})
meas.err = sqrt(meas.stat^2 + meas.syst^2)

##--- unshifted values
meas.val.orig = sapply(measurements[meas.used], function(x) {unname(x$value.orig)})
meas.syst.orig = sapply(measurements[meas.used], function(x) {unname(x$syst.orig)})

##--- which measurements got shifted in value or syst. error
meas.shifted = (meas.val.orig != meas.val) | (meas.syst.orig != meas.syst)

if (FALSE && any(meas.shifted)) {
  cat("\n##\n")
  cat("## following measurements were shifted from updated external parameters\n")
  cat("##\n")
  rc = alu.rbind.print(
    rbind(orig=meas.val.orig[meas.shifted],
          value=meas.val[meas.shifted],
          ## stat=meas.stat[meas.shifted],
          orig=meas.syst.orig[meas.shifted],
          syst=meas.syst[meas.shifted]),
    num.columns=1)
}
  
##--- print updated measurements
if (FALSE) {
  cat("\n##\n")
  cat("## measurements after update due to updated parameters\n")
  cat("##\n")
  rc = alu.rbind.print(
    rbind(value=meas.val,
          stat=meas.stat,
          syst=meas.syst,
          error=meas.err),
    num.columns=1)
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
## warn about correlations with unknown measurements
##
## if strict=FALSE, deal with alucomb.r - style correlations terms,
## which have just 3 fields to identify the measurement
## and do not match the complete measurement name
##

##--- utility function to replace incomplete measurement names
alu.find.matching.meas = function(name, strict=FALSE) {
  regexp = paste("^", gsub("[.]", "[.]", name, perl=TRUE), ".*", ifelse(strict, "$", ""), sep="")
  mapply(function(name, x) {
    rc = grep(x, meas.names, perl=TRUE, value=TRUE)
    if (length(rc)>1) {
      stop("measurement", mi.name, "\n  correlation", name, "\n  matches multiple measurements\n", rc)
    }
    if (length(rc) == 0) {
      cat("warning, measurement", mi.name, "\n  correlation", name, "\n  does not match any measurement\n")
      return(name)
    }
    cat("warning, measurement", mi.name, "\n  correlation", name, "\n  updated with matching measurement", rc[1], "\n")
    return(rc[1])
  },
  name, regexp)
}

flag.empty.line.first.time <<- TRUE
strict.correlation.matching = FALSE
##--- replace incomplete measurement names in correlations
for (mi.name in meas.names) {
  corr.names = names(measurements[[mi.name]]$corr.terms.stat)
  corr.names.missing.mask = !(corr.names %in% meas.names)
  if (any(corr.names.missing.mask)) {
    print.empty.line.first.time()
    corr.names.missing = corr.names[corr.names.missing.mask]
    corr.names.missing.upd = alu.find.matching.meas(corr.names.missing, strict.correlation.matching)
    names(measurements[[mi.name]]$corr.terms.stat)[corr.names.missing.mask] = corr.names.missing.upd
  }
  
  corr.names = names(measurements[[mi.name]]$corr.terms.tot)
  corr.names.missing.mask = !(corr.names %in% meas.names)
  if (any(corr.names.missing.mask)) {
    print.empty.line.first.time()
    corr.names.missing = corr.names[corr.names.missing.mask]
    corr.names.missing.upd = alu.find.matching.meas(corr.names.missing, strict.correlation.matching)
    names(measurements[[mi.name]]$corr.terms.tot)[corr.names.missing.mask] = corr.names.missing.upd
  }
}

##
## set off-diagonal correlation matrix coefficients from cards
## - meas.corr.stat means only stat. correlation, to be multiplied by stat. errors
## - meas.corr means total correlation, to be multiplied by total errors
##

##--- set off-diagonal statistical correlation matrix coefficients from cards
for (mi.name in meas.names) {
  for (mj.name in intersect(names(measurements[[mi.name]]$corr.terms.stat), meas.names)) {
    meas.corr.stat[meas.names %in% mi.name, meas.names %in% mj.name] = measurements[[mi.name]]$corr.terms.stat[[mj.name]]
  }
  for (mj.name in intersect(names(measurements[[mi.name]]$corr.terms.tot), meas.names)) {
    meas.corr[meas.names %in% mi.name, meas.names %in% mj.name] = measurements[[mi.name]]$corr.terms.tot[[mj.name]]
  }
}

flag = FALSE

##--- check that the STAT_CORR_WITH terms are symmetric
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

##--- check that the ERROR_CORR_WITH terms are symmetric
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
flag = FALSE
for (i in 1:meas.num) {
  for (j in i:meas.num) {
    if (meas.corr[i,j] != 0 && meas.corr.stat[i,j] != 0) {
      flag = TRUE
      cat(paste("error: both total and statistical correlation specified for measurements:\n  ",
                meas.names[i], ", ", meas.names[j], "\n", collapse=""))
    }
  }
}
if (flag) stop("aborted because of above errors\n")

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

##--- if total correlation specified, get stat. correlation by subtraction (term-by-term)
meas.cov.stat = ifelse(meas.cov == 0, meas.cov.stat, meas.cov - meas.cov.syst)
meas.cov.stat = meas.cov.stat + diag.m(meas.stat^2)
##--- total covariance
meas.cov.syst = meas.cov.syst + diag.m(meas.syst^2)
meas.cov = meas.cov.stat + meas.cov.syst

##--- total correlation
meas.corr = meas.cov / (meas.err %o% meas.err)

##
## error scaling using the "scale" parameter for fitted quantities
##
quant.cards.sfact = unlist(lapply(combination$quantities[quant.names], function(el) { unname(el["scale"]) }))
if (!is.null(quant.cards.sfact)) {
  meas.scale.names = names(meas.quantities[meas.quantities %in% names(quant.cards.sfact)])
  meas.scale.stat = sapply(measurements[meas.used][meas.scale.names], function(x) {unname(x$stat)})
  meas.scale.syst = sapply(measurements[meas.used][meas.scale.names], function(x) {unname(x$syst.orig)})
  
  ##--- compute additional syst. error to scale the total original error as requested
  meas.scale.systsq = (quant.cards.sfact^2 -1) * (meas.scale.stat^2 + meas.scale.syst^2)
  ##--- add additional syst. contribution to the diagonal elements of the covariance
  diag(meas.cov)[meas.scale.names] = diag(meas.cov)[meas.scale.names] + meas.scale.systsq
  
  cat("\n")
  for (quant.i in names(quant.cards.sfact)) {
    cat("applying s-factor =", quant.cards.sfact[quant.i], "for quantity", quant.i, "in measurements:\n")
    cat(names(meas.quantities[meas.quantities %in% quant.i]), sep="\n")
  }
}

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
delta = as.matrix(sapply(quant.names, function(x) as.numeric(x == meas.quantities)))
rownames(delta) = meas.names
##--- needed when length(quant.names) == 1
colnames(delta) = quant.names

##--- print corrected measurements
if (FALSE) {
  cat("\n##\n")
  cat("## measurements after updating systematic parameters\n")
  cat("##\n")
  rc = alu.rbind.print(
    rbind(value=meas.val,
          stat=meas.stat,
          syst=meas.syst,
          error=meas.err),
    num.columns=1)
}

##
## preliminary computations for analytical minimum chi-square solution
##
quant.num = length(quant.names)
quant.val = rep(0, quant.num)
names(quant.val) = quant.names
meas.invcov = solve(meas.cov)
meas.invcov = (meas.invcov + t(meas.invcov))/2
## meas.invcov = chol2inv(chol(meas.cov))

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
delta.measured = delta[, quant.measured.bool]
##--- fit for quantities without constraints
quant.seed.val = quant.val
quant.seed.val[quant.measured.bool] =
  solve(t(delta.measured) %*% meas.invcov %*% delta.measured) %*%
  (t(delta.measured) %*% meas.invcov %*% meas.val)

##--- get seed values for quantities without measurements
quant.cards.seed.val = unlist(lapply(combination$quantities[quant.names], function(el) { unname(el["seed"]) }))
if (is.null(quant.cards.seed.val)) {
  quant.cards.seed.val = numeric(0)
}

##--- set seed values for quantities that have no measurements
seed.needed = setdiff(quant.names, quant.names[quant.measured.bool])
seed.needed.notincards = setdiff(seed.needed, names(quant.cards.seed.val))
if (length(seed.needed.notincards)>0) {
  cat("\nwarning: no seed value provided for the following quantities\n")
  print(seed.needed.notincards)
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
    tmp = mapply(function(name, val, comb) {names(val) = name; print(c(val, unlist(comb)))},
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
cat("##\n\n")
alu.rbind.print(rbind(value=quant.val[1:quant.num], error=quant.err[1:quant.num]))
if (FALSE && quant.num > 1) {
  cat("\ncorrelation\n\n")
  print(quant.corr[1:quant.num, 1:quant.num])
}
cat("\n## end\n")

##--- save data and results
rc = save(file=file.name.data,
  measurements, combination, delta,
  chisq, dof,
  meas.val,   meas.err,   meas.cov,  meas.cov.stat, meas.cov.syst, meas.corr,
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

##--- derivative and gradient of non-linear constraint equation expressions
constr.nl.expr = lapply(combination$constr.all.expr[combination$constr.all.nl], function(x) deriv(x, all.vars(x)))

##--- print linearized constraints
print.linearized.constraints = function(val, comb) {
  tmp = mapply(function(val, comb) {
    cat(sprintf("%.5g", val), "=", paste(mapply(function(name, val) {
      paste(sprintf("%.5g", val), name, sep="*")
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
  constr.nl.val = lapply(constr.nl.expr, function(x) eval(x, as.list(quant.val)))
  constr.nl.comb = lapply(constr.nl.val, function(x) drop(attr(x, "gradient")))
  constr.nl.val = unlist(combination$constr.nl.val) - unlist(constr.nl.val)
  if (any(combination$constr.all.nl)) {
    constr.nl.val = constr.nl.val + sapply(constr.nl.comb, function(x) drop(x %*% quant.val[names(x)]))
  }

  ##--- join linear constraints
  constr.all.comb = c(combination$constr.lin.comb, constr.nl.comb)
  constr.all.val = c(unlist(combination$constr.lin.val), constr.nl.val)
  ##--- convert linearized constraint equations into linear matricial form
  constr.m = do.call(rbind, lapply(constr.all.comb, function(x) {tmp = quant.val*0; tmp[names(x)] = x; tmp}))
  constr.v = constr.all.val
  
  if (constr.num > 0) {
    if (first.iteration) {
      ##--- determine the typical size of quant.invcov elements
      sv = svd(quant.invcov)$d
      sv.central = round(quant.num*1/3):round(quant.num*2/3)
      sv.log.mean = mean(log(sv[sv.central]))
      quant.invcov.order = 10^round(sv.log.mean/log(10))
      ##--- determine the typical size of the constraint equation terms
      constr.m.order = apply(constr.m, 1, function(x) 10^round(log(mean(abs(x[x!=0])))/log(10)))

      if (TRUE) {
        cat("\n## Begin of linearized constraint equations (1st iteration)\n\n")
        print.linearized.constraints(constr.nl.val, constr.nl.comb)
        cat("\n## End of linearized constraint equations (1st iteration)\n")
      }

      cat("\n## Begin of constraint percent change summaries\n")
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
  
  if (all(!combination$constr.all.nl)) {
    ##--- there is no non-linear constraint whose linearization must iteratively converge
    break
  }

  if (!first.iteration) {
    ##--- compute sum of squares of differences between previous and current iteration linearized constraints
    constr.diff = mapply(function(c2, v2, c1, v1) {
      x2 = c(c2, v2)
      x1 = c(c1, v1)
      norm = pmax(abs(x1),abs(x2),abs(x2+x1))
      norm = mean(norm)
      ifelse(norm==0, 0, sum(((x2-x1)/norm)^2))
    }, constr.nl.comb, constr.nl.val, constr.nl.prev.comb, constr.nl.prev.val)
    cat("\n")
    print(constr.diff, digits=3)
    ##--- end if the average percent change of constraint equation coefficients is small enough
    if (sum(constr.diff) / length(constr.nl.val) < 1e-10) break
  }

  ##--- save before calculation of next updated values
  constr.nl.prev.comb = constr.nl.comb
  constr.nl.prev.val = constr.nl.val
  first.iteration = FALSE
}
if (constr.num > 0) {
  cat("\n## End of constraint percent change summaries\n")
}
rm(first.iteration)
suppressWarnings(rm(constr.diff))

if (TRUE) {
  cat("\n## Begin of linearized constraint equations (after convergence)\n\n")
  print.linearized.constraints(constr.nl.val, constr.nl.comb)
  cat("\n## End of linearized constraint equations (after convergence)\n")
}

##
## compute errors and chi square
##
solve.cov.m = solve.m[1:quant.num,1:meas.num, drop=FALSE]
quant.cov = solve.cov.m %*% meas.cov %*% t(solve.cov.m)
quant.cov = (quant.cov + t(quant.cov))/2
quant.err = sqrt(diag(quant.cov))
quant.corr = quant.cov / (quant.err %o% quant.err)

chisq = drop(t(meas.val - delta %*% quant.val) %*% meas.invcov %*% (meas.val - delta %*% quant.val))
dof = meas.num - quant.num + constr.num

cat("\n##\n")
cat("## alucomb2 solution, chisq/d.o.f. = ",chisq, "/", dof, ", CL = ", (1-pchisq(chisq, df=dof)), "\n",sep="")
cat("##\n\n")
rc = alu.rbind.print(rbind(value=quant.val, error=quant.err))
if (FALSE && quant.num > 1) {
  cat("\ncorrelation\n\n")
  print(quant.corr)
}
cat("\n## end\n")

##--- save data and results
rc = save(file=file.name.data,
  measurements, combination, delta,
  chisq, dof,
  meas.val,   meas.err,   meas.cov,  meas.cov.stat, meas.cov.syst, meas.corr,
  quant.val,  quant.err,  quant.cov, quant.corr,
  constr.m, constr.v,
  solve.cov.m)

} # end if method alucomb2

##
## ////////////////////////////////////////////////////////////////////////////
##

if (method == "solnp") {

require(truncnorm, quietly=TRUE)
require(Rsolnp, quietly=TRUE)

##--- degrees of freedom
dof = meas.num - quant.num + length(combination$constr.all.val)

##--- function returns vector of non-linear constraint functions
constr.fun <- function(qv)  {
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
quant.hessian = rc.solnp$hessian
quant.cov = solve(quant.hessian)
rm(quant.hessian)
rownames(quant.cov) = quant.names
colnames(quant.cov) = quant.names
quant.err = sqrt(diag(quant.cov))
quant.corr = quant.cov / (quant.err %o% quant.err)

cat("\n")
cat("##\n")
cat("## solnp solution, chisq/d.o.f. = ", chisq, "/", dof, ", CL = ", (1-pchisq(chisq, df=dof)), "\n",sep="")
cat("##\n\n")
rc = alu.rbind.print(rbind(value=quant.val, error=quant.err))

if (FALSE && quant.num > 1) {
  cat("\ncorrelation\n\n")
  print(quant.corr)
}
cat("\n## end\n")

##--- cleanup
if (FALSE) {
  rm(constr.nl.expr, rc.solnp)
}

##--- save data and results
rc = save(file=file.name.data,
  measurements, combination, delta,
  chisq, dof,
  meas.val,   meas.err,   meas.cov,  meas.cov.stat, meas.cov.syst, meas.corr,
  quant.val,  quant.err,  quant.cov, quant.corr)

} # end if method solnp

##
## ////////////////////////////////////////////////////////////////////////////
##

if (method == "alabama") {

require(alabama, quietly=TRUE)
##
## note that we can give this the symbolic derivative of half.chisq.fun 
## and the jacobian of the constraints if we want to; will go faster 
##

##--- degrees of freedom
dof = meas.num - quant.num + length(combination$constr.all.val)

##
## function returns vector of non-linear constraint functions
## here we include the combination$constr.all.val for alabama
## meaning: [ f(x) = const ] goes to [ f(x) - const = 0 ]
##
constr.fun <- function (qv)  {
  env.list = as.list(qv)
  rc = unlist(lapply(combination$constr.all.expr, function(expr) {eval(expr, env.list)}))
  rc = rc - combination$constr.all.val
  return(rc)
}

##--- function returns 1/2 of chi square
half.chisq.fun <- function(qv) {
  drop(t(meas.val - delta %*% qv) %*% meas.invcov %*% (meas.val - delta %*% qv))/2
}

##--- minimize half chi square to get the correct covariance for the fitted parameters
cat("\nBegin of alabama minimization\n")
rc.alabama = auglag(par = quant.seed.val, fn = half.chisq.fun, heq = constr.fun)
cat("\nEnd of alabama minimization\n")

chisq = 2*tail(rc.alabama$value, 1)
quant.val = rc.alabama$par
quant.cov = solve(rc.alabama$hessian)
rownames(quant.cov) = quant.names
colnames(quant.cov) = quant.names
quant.err = sqrt(diag(quant.cov))
quant.corr = quant.cov / (quant.err %o% quant.err)

cat("\n")
cat("##\n")
cat("## alabama solution, chisq/d.o.f. = ", chisq, "/", dof, ", CL = ", (1-pchisq(chisq, df=dof)), "\n",sep="")
cat("##\n\n")
rc = alu.rbind.print(rbind(value=quant.val, error=quant.err))
if (FALSE && quant.num > 1) {
  cat("\ncorrelation\n\n")
  rc = alu.rbind.print(quant.corr)
}
cat("\n## end\n")

##--- cleanup
if (FALSE) {
  rm(rc.alabama)
}

##--- save data and results
rc = save(file=file.name.data,
  measurements, combination, delta,
  chisq, dof,
  meas.val,   meas.err,   meas.cov,  meas.cov.stat, meas.cov.syst, meas.corr,
  quant.val,  quant.err,  quant.cov, quant.corr)

}  # end of if method alabama

##
## ////////////////////////////////////////////////////////////////////////////
##

cat("\n")
cat(paste("file", file.name.data, "produced\n"))
cat("\n## end\n")

options(options.save)
} ##--- end function alucomb

args <- commandArgs(TRUE)
if (length(args) > 0 && exists("alucomb")) alucomb(file = args[1]) 
