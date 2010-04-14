#!/usr/bin/env Rscript

## ////////////////////////////////////////////////////////////////////////////
##
## aluplot, create plot.input files
##
## ////////////////////////////////////////////////////////////////////////////

library(methods)

source("../../../Common/bin/alucomb-read.r")

## ////////////////////////////////////////////////////////////////////////////
## definitions

##
## transform the name of a measurement in the HFAG-tau format into a format usable by Root
##
label.root = function(str) {
  str.orig = str

  str = gsub("(^|[^A-Z])([A-Z][a-y]*)([z+-]*)b", "\\1#bar{\\2}\\3", str, perl=TRUE)
  str = gsub("F1", "f_{1}",str)
  str = gsub("Pi", "#pi",str)
  str = gsub("Nu", "#nu",str)
  str = gsub("M", "#mu",str)
  str = gsub("H", "h",str)
  str = gsub("m($|[#}A-Zh])", "^{-}\\1", str, perl=TRUE)
  str = gsub("p($|[#}A-Zh])", "^{+}\\1", str, perl=TRUE)
  str = gsub("z($|[#}A-Zh])", "^{0}\\1", str, perl=TRUE)

  if (str.orig %in% c("HmHmHpNu", "PimKmPipNu", "PimPimPipNu")) {
    str = paste(str, "(ex.K^{0})")
  }

  str = paste("B(#tau^{-} #rightarrow ", str, ")", sep="")
  return(str)
}

##
## read log file pdg_average.log containing the PDG-like average
##
get.pdg.average = function(file) {
  rc = list()
  fh = pipe(paste(
    "awk 2>/dev/null '/Total Chi2 for/,EOF'",
    file, sep=" "))
  average = readLines(fh)
  close(fh)
  if (length(average) <= 0) return(list())
  average = paste(average, collapse=";")
  patt = paste(".*Scale factor = ([.0-9eE+-]+)",
    ".*.*Weighted average .with scale factor. = ([.0-9eE+-]+)[\\s+-]+([.0-9eE+-]+)",
    sep="")
  if (regexpr(patt, average, perl=TRUE) == -1) {
    warning("Cannot read PDG average from ", file)
    return(list())
  }
  fields = unlist(strsplit(sub(patt, "\\1;\\2;\\3", average, perl=TRUE), ";"))
  rc$scale = fields[1]
  rc$value = as.numeric(fields[2:3])
  return(rc)
}

##
## read log file of Combos using chi2_sym procedure
##
get.combos.chi2sym = function(file) {
  rc = list()
  fh = pipe(paste(
    "awk 2>/dev/null '/CHI2_SYM: CHI2\\/NDOF, SCALE FAC, CL/ {flag++}; flag>1'",
    file, sep=" "))
  average = readLines(fh)
  close(fh)
  if (length(average) <= 0) return(list())
  average = paste(average, collapse=";")
  patt = paste(".*CHI2_SYM\\s+:\\s+([^=\\s]+)\\s*=\\s*([.0-9eE+-]+)",
    "\\s*\\+-\\s*([.0-9eE+-]+)\\s*CL\\s*=\\s*([.0-9eE+-]+).*", sep="")
  if (regexpr(patt, average, perl=TRUE) == -1) {
    return(list())
  }
  fields = unlist(strsplit(sub(patt, "\\1;\\2;\\3;\\4", average, perl=TRUE), ";"))
  rc$CL = as.numeric(fields[4])
  ##++ need to recover also sqrt(chisq/dof)
  rc$scale = -rc$CL
  label = sub("^M_", "", fields[1])
  val = as.numeric(fields[2])
  err = as.numeric(fields[3])
  names(val) = label
  names(err) = label
  rc$val = val
  rc$err = err
  return(rc)
}

##
## read log file of Combos using chi2_nsym procedure
##
get.combos.chi2nsym = function(file) {
  rc = list()
  fh = pipe(paste(
    "awk 2>/dev/null '/CHI2_N_SYM: CHI2, NMEFF, NQUAN, NDOF/ {flag++}; flag>1'",
    file, sep=" "))
  average = readLines(fh)
  close(fh)
  if (length(average) <= 0) return(list())
  patt0 = paste(".*CHI2_N_SYM:\\s+CHI2,\\s+NMEFF,\\s+NQUAN,\\s+NDOF\\s*=",
    "\\s*([.0-9eE+-]+)\\s+([.0-9eE+-]+)\\s+([.0-9eE+-]+)\\s+([.0-9eE+-]+).*", sep="")
  if (regexpr(patt0, average[1]) == -1) return(list())
  fields = unlist(strsplit(sub(patt0, "\\1;\\2;\\3;\\4", average[1], perl=TRUE), ";"))
  rc$chisq = as.numeric(fields[1])
  rc$dof = as.numeric(fields[4])
  rc$CL = 1 - pchisq(rc[["chisq"]], rc[["dof"]])
  rc$scale = sqrt(rc$chisq/rc$dof)

  ##-- first line containig combos averages
  patt = paste(".*CHI2_N_SYM\\s+:\\s+([^=\\s]+)\\s*=\\s*([.0-9eE+-]+)",
    "\\s*\\+-\\s*([.0-9eE+-]+)\\s*CL\\s*=\\s*([.0-9eE+-]+).*", sep="")
  ##-- following lines
  patt2 = paste(".*\\s+([^=\\s]+)\\s*=\\s*([.0-9eE+-]+)",
    "\\s*\\+-\\s*([.0-9eE+-]+)\\s*.*", sep="")

  flag = FALSE  
  comb.names = character()
  val = numeric()
  err = numeric()
  for (line in average) {
    if (regexpr(patt, line) != -1) {
      fields = unlist(strsplit(sub(patt, "\\1;\\2;\\3;\\4", line, perl=TRUE), ";"))
      flag = TRUE
    } else if (flag) {
      if (regexpr(patt2, line) != -1) {
        fields = unlist(strsplit(sub(patt2, "\\1;\\2;\\3", line, perl=TRUE), ";"))
      } else {
        flag = FALSE
      }
    } else {
      flag = FALSE
    }
    if (flag) {
      quant = fields[1]
      quant = sub("^M_", "", quant)
      comb.names = c(comb.names, quant)
      val = c(val, as.numeric(fields[2]))
      err = c(err, as.numeric(fields[3]))
    }
  }
  names(val) = comb.names
  names(err) = comb.names
  rc$val = val
  rc$err = err
  return(rc)
}

##
## get Combos results in specified file
##
get.combos.results = function(file) {
  lines = suppressWarnings(try(get.file.lines(file), silent=TRUE))
  if (inherits(lines, "try-error")) {
    warning("Cannot open / read file ", file)
    return(list())
  }
  rc.chi2sym = get.combos.chi2sym(file)
  rc.chi2nsym = get.combos.chi2nsym(file)
  if (length(rc.chi2nsym) >0) {
    if (length(rc.chi2sym) >0) {
      warning("both chi2_sym and chi2_nsym results in ", file)
    }
    return(rc.chi2nsym)
  } else if (length(rc.chi2sym) >0) {
    return(rc.chi2sym)
  }
  warning("Cannot get Combos information from", file)
  return(list())
}

##
## get numeric data from aluxomb log
##
get.alucomb.data = function(lines) {
  data = numeric()
  li = 1
  if (is.na(lines[li])) {
    return(list(lines=li-1, data=data))
  }
  fields = unlist(strsplit(lines[li], "\\s+", perl=TRUE))
  if (fields[1] != "") {
    warning("Cannot find data header line in ", lines[li])
  }
  cols = fields[-1]
  rows = character()
  repeat {
    li = li+1
    if (is.na(lines[li])) {
      break
    }
    fields = unlist(strsplit(lines[li], "\\s+", perl=TRUE))
    prev.options = options()
    options(warn=-1)
    if (length(fields) <= length(cols)) break
    values = as.numeric(fields[-(1:(length(fields)-length(cols)))])
    options(prev.options)
    if (length(values) == 0 || any(is.na(values))) break
    data = rbind(data, values)
    rows = c(rows, paste(fields[1:(length(fields)-length(cols))], collapse=" "))
  }
  if (length(data) > 0) {
    if (is.null(dim(data))) {
      names(data) = cols
    } else {
      colnames(data) = cols
      rownames(data) = rows
    }
  }
  return(list(lines=li-1, data=data))
}

##
## get a section of numeric data from aluxomb log
##
get.alucomb.section = function(lines, pattern, file, offset=0) {
  liv = grep(pattern, lines)
  if (length(liv) == 0) {
    cat("warning: cannot find '", pattern, "' in\n  ", file, "\n", sep="")
    return(NULL)
  }
  li = liv[1]
  li = li+1+offset
  
  rc = get.alucomb.data(lines[-(1:(li-1))])
  if (length(data) == 0) {
    warning("could not read chisq data in ", file)
    return(NULL)
  }
  return(rc)
}

##
## get alucomb log data
##
get.alucomb = function(file) {
  ol = list()
  
  lines = suppressWarnings(try(get.file.lines(file), silent=TRUE))
  if (inherits(lines, "try-error")) {
    warning("Cannot open / read file ", file)
    return(ol)
  }
  
  rc = get.alucomb.section(lines, "^#+\\s+S-factors accounting", file, 1)
  if (is.null(rc)) {
    return(ol)
  }
  li = 1 + rc$lines
  ol$chisq = rc$data["no-large-error", "chisq"]
  ol$dof = rc$data["no-large-error", "dof"]
  
  rc.av = get.alucomb.section(lines[-(1:(li-1))], "^Averaged quantities", file)
  if (is.null(rc)) {
    return(ol)
  }
  li = li + rc$lines
  
  rc.corr = get.alucomb.section(lines[-(1:(li-1))], "^correlation", file)
  if (is.null(rc)) {
    return(ol)
  }
  li = li + rc$lines
  
  rc.extra = get.alucomb.section(lines[-(1:(li-1))], "^Non-averaged measurement types", file)
  if (is.null(rc)) {
    return(ol)
  }
  li = li + rc$lines

  data = cbind(rc.av$data, rc.extra$data)
  
  ol$val = data["value",]
  names(ol$val) = colnames(data)
  ol$err = data["error",]
  names(ol$err) = colnames(data)
  ol$sfact = data["S-factor",]
  names(ol$sfact) = colnames(data)
  ol$chisq.all = data["chisq",]
  names(ol$chisq.all) = colnames(data)
  ol$dof.all = data["dof",]
  names(ol$dof.all) = colnames(data)

  return(ol)
}

## ////////////////////////////////////////////////////////////////////////////
## code

file = "average.input"

##++ aluplot = function(file = "") {

rc = alucomb.read(file)
measurements = rc$measurements
combination = rc$combination
plot.data = list()

log.dir = file.path("../../../Data", sub("^.*/([^/]*/[^/]*/[^/]*)$", "\\1", getwd()))
log.file = function(fname) {file.path(log.dir, fname)}

##-- get quantities measured by each experiment
meas.names = names(measurements)
meas.num = length(meas.names)

##++ gets names in Combos MEASUREMENT card, sometimes incorrect
meas.quantities.raw = unlist(lapply(measurements, function(x) names(x$value)))
##++ get accurate names for quantities measured by experiments
meas.quantities = sub("[.][^.]+$", "", sub("^[^.]+[.]", "", names(measurements)))

##-- quantities we want to determine from measurements
quant.names = combination$quantities

##-- measurements of quantities that are linear combination of the quantities to be averaged
meas.names.comb = names(combination$meas.lin.combs[unlist(
  lapply(combination$meas.lin.combs, function(el) all(names(el) %in% combination$quantities)))])
##-- quantity measured by each measurement
meas.quantities = unlist(lapply(measurements, function(x) names(x$value)))
names(meas.quantities) = names(measurements)
##-- unique quantities corresponding to linear combinations of averaged quantities
quant.names.comb = unique(meas.quantities[meas.names.comb])
quant.names.comb = setdiff(quant.names.comb, quant.names)
##-- will plot all averaged quantities and all their measured combinations
quant.names = c(quant.names, quant.names.comb)

##-- get list of stat and syst errors
meas.value = unlist(lapply(measurements, function(x) x$value))
names(meas.value) = meas.names
meas.stat = unlist(lapply(measurements, function(x) x$stat))
names(meas.stat) = meas.names
meas.syst = unlist(lapply(measurements, function(x) x$syst))
names(meas.syst) = meas.names
##-- combine true and fake info, compute total error
meas.error = sqrt(meas.stat^2 + meas.syst^2)

##-- init meas.names.true to use same code as in alucomb
meas.names.true = meas.names

##++ begin of code in common with alucomb.r

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

##++ end of code in common with alucomb.r

##
## get syst. correlation corresponding to correlated syst. terms
##
meas.cov.syst = meas.corr * 0
for (meas.i in meas.names.true) {
  syst.i = measurements[[meas.i]]$syst.terms
  for (meas.j in meas.names.true) {
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
## get PDG average produced from the current directory
## (not used any more, a dedicated file is used for this information)
##
rc = get.pdg.average(log.file("pdg_average.log"))
if (length(rc) > 0) {
   plot.data$pdg = rc
}

##-- read dedicated file with PDG 2009 averages
lines = get.file.lines("../Common/pdg_averages.input")
pdg.averages = list()
for (line in lines) {
  ##-- remove comments up to end of line, and preceding space
  line = gsub("\\s*[*#;].*$", "", line, perl=TRUE)
  if (line == "") {
    next
  }
  ##-- remove leading space
  line = gsub("^\\s+", "", line, perl=TRUE)
  fields = unlist(strsplit(line, "\\s+", perl=TRUE))
  pdg.averages[[fields[1]]] = as.numeric(fields[-1])
}

##
## in the Combos cards the measurement results have symmetrized errors
## override symmetrized errors with the original ones usind a dedicated file
##
lines = get.file.lines("../Common/results_asymm_errors.input")
meas.asymm = list()
for (line in lines) {
  if (regexpr("^\\s*$", line, perl=TRUE) != -1 ||
      regexpr("^[*#;]", line, perl=TRUE) != -1) {
    next
  }
  line = gsub("\\s*!.*", "", line, perl=TRUE)
  fields = unlist(strsplit(line, "\\s+", perl=TRUE))
  tag = paste(fields[1:3], collapse=".")
  meas.asymm[[tag]] = as.numeric(fields[-(1:3)])
  names(meas.asymm[[tag]]) = c("meas", "stat_pos", "stat_neg", "syst_pos", "syst_neg")
}

for (quant in quant.names) {
  plot.data[[quant]] = list()
}

##-- get alucomb results
rc = get.alucomb(log.file("average_alucomb.log"))
if (length(rc) > 0) {
  plot.data$hfag = rc
} else {
  ##-- get Combos results
  rc = get.combos.results(log.file("average.log"))
  if (length(rc) > 0) {
    plot.data$hfag = rc
  }
}

##
## for each quantity to be plotted, collect in plot.data[[quant]]
## - minimum and maximum x coordinate for plotting
## - the best units/ precision for plotting/ printing
## - the list of the experiments that measured it
##
for (quant in quant.names) {
  x.mins = numeric()
  x.maxs = numeric()
  x.values = numeric()
  all.exp.meas = list()
  meas.list.tmp = list()

  ##-- collect value, stat, syst, bibitem of all relevant measurements
  for (meas in meas.names[quant == meas.quantities]) {
    meas.tmp = list()
    meas.tmp[[meas]] =
      list(value=measurements[[meas]]$value,
           stat=measurements[[meas]]$stat,
           syst=measurements[[meas]]$syst,
           bibitem=measurements[[meas]]$bibitem,
           index=which(meas.names == meas))
    meas.list.tmp = c(meas.list.tmp, meas.tmp)
  }

  if (quant == "HmHmHpNu") {
    ##-- special case, add BaBar and Belle hhh measurements as sum of exclusive modes
    combnames = c("PimPimPipNu", "PimKmPipNu", "PimKmKpNu", "KmKmKpNu")
    explist = list(
      list("BaBar", "published"),
      list("Belle", "submitted"))
    reslist = unlist(lapply(explist, function(exp) {
      rc = list()
      name = paste(exp[[1]], quant, exp[[2]], sep=".")
      rc[[name]] = names(meas.quantities) %in% paste(exp[[1]], combnames, exp[[2]], sep=".")
      return(rc)
    }), recursive=FALSE)
    
    reslist.sorted.names = character(0)
    reslist.which.meas = do.call(rbind, lapply(reslist, function(x) which(x)))
    if (length(reslist.which.meas) >0) {
      reslist.sorted.names = names(sort(reslist.which.meas[,1]))
      ##-- collect value, stat, syst, bibitem of the combination
      for (combres in reslist.sorted.names) {
        combres.vec = as.numeric(reslist[[combres]])
        value = drop(combres.vec %*% meas.value)
        stat = sqrt(drop(combres.vec %*% meas.cov.stat %*% combres.vec))
        syst = sqrt(drop(combres.vec %*% meas.cov.syst %*% combres.vec))
        ##-- order number of 1st measurement used in the combination
        index = which(combres.vec !=0)[1]
        bibitem = measurements[[index]]$bibitem
        ##-- substitute quantity with the combination we calculated
        bibitem[2] = quant
        ## cat(unlist(strsplit(combres, ".", fixed=TRUE)), value, stat, syst, "\n")
        meas.tmp = list()
        meas.tmp[[combres]] = list(value=value, stat=stat, syst=syst, bibitem=bibitem, index=index)
        meas.list.tmp = c(meas.list.tmp, meas.tmp)
      }
    }
  }
  
  ##-- sort all measurements according to their ordering in the Combos cards
  for (meas in names(sort(unlist(lapply(meas.list.tmp, function(x) x$index))))) {
    exp.meas = list()
    value = meas.list.tmp[[meas]]$value
    stat = meas.list.tmp[[meas]]$stat
    syst = meas.list.tmp[[meas]]$syst
    bibitem = meas.list.tmp[[meas]]$bibitem

    error = sqrt(sum(c(stat), syst)^2)
    x.mins = c(x.mins, value - error)
    x.maxs = c(x.maxs, value + error)
    x.values = c(x.values, value)

    exp.meas$value = c(value, +stat, -stat, +syst, -syst)
    exp.meas$bibitem = paste(c(bibitem[1], bibitem[-(1:3)]), collapse=" ")

    if (!is.null(meas.asymm[[meas]])) {
      ##-- override measurement with info on asymmetric errors
      cat("ovverride", meas, "\n")
      cat("  from", sprintf("%12.6g ", exp.meas$value), "\n")
      cat("    to", sprintf("%12.6g ", meas.asymm[[meas]]), "\n")
      exp.meas$value = meas.asymm[[meas]]
    }

    all.exp.meas[[meas]] = exp.meas
  }

  ##
  ## get average info either from alucomb or Combos
  ##
  if (is.null(plot.data$hfag$chisq.all)) {
    ##-- alucomb info not available, use combos if existing
    value = plot.data$hfag$val[toupper(quant)]
    if (is.null(value)) {
      stop("Neither alucomb nor Combos log files were found")
    }
    error = plot.data$hfag$err[toupper(quant)]
    conf.lev = plot.data$hfag$CL
    scale.factor = plot.data$hfag$scale
  } else {
    ##-- alucomb info available
    value = plot.data$hfag$val[quant]
    error = plot.data$hfag$err[quant]
    chisq = plot.data$hfag$chisq.all[quant]
    dof = plot.data$hfag$dof.all[quant]
    scale.factor = plot.data$hfag$sfact[quant]
    ##-- in case of multiple quantities average, consider global CL
    ## conf.lev = pchisq(chisq, df=dof, lower.tail=FALSE)
    conf.lev = pchisq(plot.data$hfag$chisq, df=plot.data$hfag$dof, lower.tail=FALSE)
  }
  
  ##
  ## HFAG does not usually use S-factors but rather quotes CL
  ## however when CL<CL_min we quote S-factors instead
  ## CL_min is set to the one-sided Gaussian fraction beyond 3 sigma
  ##
  conf.lev.plot = conf.lev
  if (conf.lev < pnorm(3, mean=0, sd=1, lower.tail=FALSE)) {
    ##++ remove following "if" when scale.factor for chi_sym fixed
    if (conf.lev != -scale.factor) {
      error = error*scale.factor
    }
    conf.lev.plot = -scale.factor
  }
  
  if (is.null(value) || is.na(value)) {
    warning("could not find HFAG average data for ", quant)
  } else {
    plot.data[[quant]]$hfag$avg = c(value, error)
    plot.data[[quant]]$hfag$CL = conf.lev
    plot.data[[quant]]$hfag$CL.plot = conf.lev.plot
    plot.data[[quant]]$hfag$scale = scale.factor
    x.mins = c(x.mins, value - error)
    x.maxs = c(x.maxs, value + error)
  }
  
  ##-- get xmin/xmax of PGD average if existing
  tmp = pdg.averages[[quant]]
  if (is.null(tmp)) {
    warning("could not find PDG average for ", quant)
  } else {
    x.mins = c(x.mins, tmp[1] - tmp[2])
    x.maxs = c(x.maxs, tmp[1] + tmp[2])
    plot.data[[quant]]$pdg$avg = tmp[1:2]
    plot.data[[quant]]$pdg$scale = tmp[3]
    ##-- in plot.input set scale = 0 if PDG used no scale factor
    if (tmp[3] == 1) {
      plot.data[[quant]]$pdg$scale = 0
    }
  }

  ##-- compute plot size, units, precision
  x.min = min(x.mins)
  x.max = max(x.maxs)
  x.mean = (x.min+x.max)/2
  not.higher.order = floor(log(max(x.values)*1.01)/log(10))
  if (not.higher.order == -1) {
    order = -2
    precision = 5.2
  } else if (not.higher.order == -2) {
    order = -2
    precision = 5.3
  } else if (not.higher.order == -3) {
    order = -2
    precision = 5.3
  } else if (not.higher.order == 0) {
    order = 0
    precision = 5.3
  } else if (not.higher.order == 1) {
    order = 0
    precision = 5.2
  } else if (not.higher.order == 2) {
    order = 0
    precision = 5.1
  } else if (not.higher.order == 3) {
    order = 0
    precision = 5.0
  } else {
    order = not.higher.order
    precision = 5.3
  }
  ## x.fraction.result = 2.5/9.5
  ## x.fraction.result = 1/3
  x.fraction.result = 0.3
  x.padding.left = 0.05 + (1/x.fraction.result -1)
  x.padding.right = 0.05
  x.min = x.min - x.padding.left*(x.max-x.min)
  x.max = x.max + x.padding.right*(x.max-x.min)
  plot.data[[quant]]$xmin = x.min
  plot.data[[quant]]$xmax = x.max
  plot.data[[quant]]$expts = all.exp.meas
  plot.data[[quant]]$order = order
  plot.data[[quant]]$units = exp(order*log(10))
  plot.data[[quant]]$precision = precision
}

##
## create plot.input
##
for (quant in quant.names) {
  fname = "plot.input"
  if (length(quant.names) > 1) {
    fname = paste("plot-", quant, ".input", sep="")
  }
  fh = file(fname, "w")
  quant.data = plot.data[[quant]]
  x.units = quant.data$units
  cat("# first line is xmin, xmax, precision, title\n", file=fh)
  units.label = sprintf(" [#times10^{%d}]", quant.data$order)
  if (quant.data$order == -2) {
    units.label = " [%]"
  } else if (quant.data$order == 0) {
    units.label = ""
  }

  range = quant.data$xmax/x.units - quant.data$xmin/x.units
  range.order = 10^round(log(range)/log(10))
  range.order5 = ifelse(abs(log((range.order/2)/range)) < abs(log((range.order*5)/range)), range.order/2, range.order*5)
  range.prec = ifelse(abs(log(range.order/range)) < abs(log(range.order5/range)), range.order, range.order5)
  range.prec = range.prec/100

  xmin = quant.data$xmin/x.units
  xmin = round(xmin/range.prec)*range.prec
  xmin.prec = round(log(abs(xmin/range.prec))/log(10)+0.51)
  xmin.fmt = sprintf("%%-12.%dg ", xmin.prec)

  xmax = quant.data$xmax/x.units
  xmax = round(xmax/range.prec)*range.prec
  xmax.prec = round(log(abs(xmax/range.prec))/log(10)+0.51)
  xmax.fmt = sprintf("%%-12.%dg ", xmax.prec)
  
  ## cat("* ", sprintf("%12.2g ", c(quant.data$xmin/x.units, quant.data$xmax/x.units)))
  cat("* ",
      sprintf(paste(xmin.fmt, xmax.fmt, sep=""), xmin, xmax),
      as.character(quant.data$precision), " ", label.root(quant), units.label, "\n", file=fh, sep="")
  cat("# next lines are average, error, CL (or -ScaleFactor) for HFAG Averages\n", file=fh)

  ##-- print HFAG averages
  comb = plot.data[[quant]]$hfag$avg
  if (!is.null(comb)) {
    label.extra = ""
    ##
    ## special treatment for hhh results, add hhh global HFAG average by extracting
    ## the relevant line in a previously produced plot.input file
    ## moreover, add proper information after "HFAG average"
    ##
    if (any(quant == c("PimPimPipNu", "PimKmPipNu", "PimKmKpNu", "KmKmKpNu", "HmHmHpNu"))) {
      if (length(quant.names) == 1) {
        ##-- single average of a hhh mode, get also the global hhh average result
        fname.temp = paste("../TauTo3Prongs/plot-", quant, ".input", sep="")
        cat("info: get HFAG tau -> hhh nu average from ", fname.temp, "\n", sep="")
        fhtemp = pipe(paste("grep \"^&\" ", fname.temp, " | head -1", sep=""))
        average = readLines(fhtemp)
        close(fhtemp)
        if (length(average) >0) {
          average = gsub("(\\S+\\s+\\S+\\s+\\S+\\s+\\S+)\\s*.*$", "\\1", average)
          cat(format(average, width=12*3+4),
              "HFAG Average [B(#tau^{-} #rightarrow h^{-}h^{-}h^{+}#nu) modes combined]\n", file=fh)
        }
        label.extra = paste(" [", label.root(quant), " mode only]", sep="")
      } else {
        ##-- combined average of all hhh modes
        label.extra = " [B(#tau^{-} #rightarrow h^{-}h^{-}h^{+}#nu) modes combined]"
      }
    } ## any(quant == c("PimPimPipNu", "PimKmPipNu", "PimKmKpNu", "KmKmKpNu"))

    if (any(quant == c("HmNu", "HmPizNu", "KmNu", "KmPizNu", "MmNumbarNu", "PimNu", "PimPizNu"))) {
      if (length(quant.names) == 1) {
        ##-- single average of a 1-prong mode, get also the global 1-prong average result
        fname.temp = paste("../TauTo1Prong/plot-", quant, ".input", sep="")
        cat("info: get HFAG tau -> 1-prong average from ", fname.temp, "\n", sep="")
        fhtemp = pipe(paste("grep \"^&\" ", fname.temp, " | head -1", sep=""))
        average = readLines(fhtemp)
        close(fhtemp)
        if (length(average) >0) {
          average = gsub("(\\S+\\s+\\S+\\s+\\S+\\s+\\S+)\\s*.*$", "\\1", average)
          cat(format(average, width=12*3+4),
              "HFAG Average [B(#tau^{-} #rightarrow 1-prong) modes combined]\n", file=fh)
        }
        label.extra = paste(" [", label.root(quant), " mode only]", sep="")
      } else {
        ##-- combined average of all 1-prong modes
        label.extra = " [B(#tau^{-} #rightarrow 1-prong) modes combined]"
      }
    } ## any(quant == c("HmNu", "HmPizNu", "KmNu", "KmPizNu", "MmNumbarNu", "PimNu", "PimPizNu"))
    
    cat("& ", sprintf("%-12.6g ", c(comb/x.units, plot.data[[quant]]$hfag$CL.plot)),
        "HFAG Average", label.extra, "\n", file=fh, sep="")
    cat("# next lines are average, error, Scale Factor for PDG Averages; Scale==0 means none quoted\n", file=fh)
  }

  ##-- PDG average
  comb = plot.data[[quant]]$pdg$avg
  if (!is.null(comb)) {
    cat("% ", sprintf("%-12.6g ", c(comb/x.units, plot.data[[quant]]$pdg$scale)), "PDG'09 Average\n", file=fh, sep="")
  } else {
    x.outside.plot = 2*quant.data$xmax - quant.data$xmin
    cat("% ", sprintf("%-12.6g ", c(x.outside.plot, 0, 0)/x.units), ">>> NOT FOUND <<< PDG'09 Average\n", file=fh, sep="")
  }
  
  ##-- measurements
  cat("# next lines are measurement, stat pos-error, neg-error[with negative sign],",
      "syst pos-error, neg-error[with negative sign], experiment-name\n", file=fh)
  for (exp in plot.data[[quant]]$expts) {
    bibitem = exp$bibitem
    ##-- capitalize 1st word
    bibitem = sub("^(\\S+)", "\\1", bibitem, perl=TRUE)
    ##
    ## get the year from the "where" field in the Combos cards
    ## please end the "where" field either with (YYYY) or with YYYY
    ##
    bibitem = sub("^(\\S+)\\s+.*(\\(|)(\\d\\d\\d\\d)(\\)|)$", "\\1 \\3", bibitem, perl=TRUE)
    ##--- change PDG-like years with possible letter in 4-digit years
    bibitem = sub("^(\\S+)\\s+.*\\s+([0-2]\\d)([A-Z]+|)\\s*$", "\\1 20\\2", bibitem, perl=TRUE)
    bibitem = sub("^(\\S+)\\s+.*\\s+([3-9]\\d)([A-Z]+|)\\s*$", "\\1 19\\2", bibitem, perl=TRUE)
    if (regexpr("\\d\\d\\d\\d$", bibitem, perl=TRUE) == -1) {
      cat("warning: could not find year in \"where\" Combos field for", quant, "\n")
      cat("  (", exp$bibitem, ")\n", sep="")
    }
    cat("  ", sprintf("%-12.6g ", exp$value/x.units), bibitem, "\n", file=fh, sep="")
  }
  close(fh)
  cat("file", fname, "created\n")
}

##++} ## aluplot()
