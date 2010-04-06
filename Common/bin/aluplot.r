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
  if (str == str.orig) {
    str = gsub("(^|[^A-Z])([A-Z][a-y]*)([z+-]*)b", "\\1#bar{\\2}\\3", str, perl=TRUE)
    str = gsub("F1", "f_{1}",str)
    str = gsub("Pi", "#pi",str)
    str = gsub("Nu", "#nu",str)
    str = gsub("H", "h",str)
    str = gsub("m($|[#}A-Zh])", "^{-}\\1", str, perl=TRUE)
    str = gsub("p($|[#}A-Zh])", "^{+}\\1", str, perl=TRUE)
    str = gsub("z($|[#}A-Zh])", "^{0}\\1", str, perl=TRUE)
    str = paste("B(#tau^{-} #rightarrow ", str, ")", sep="")
  }
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
    "awk 2>/dev/null '/CHI2_SYM: CHI2, NMEFF, CHI2\\/NDOF/ {flag++}; flag>1'",
    file, sep=" "))
  average = readLines(fh)
  close(fh)
  average = paste(average, collapse=";")
  patt = paste(".*CHI2_SYM\\s+:\\s+([^=\\s]+)\\s*=\\s*([.0-9eE+-]+)",
    "\\s*\\+-\\s*([.0-9eE+-]+)\\s*CL\\s*=\\s*([.0-9eE+-]+).*", sep="")
  if (regexpr(patt,average,perl=TRUE) == -1) {
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
  rc.chi2sym = get.combos.chi2sym(file)
  rc.chi2nsym = get.combos.chi2nsym(file)
  if (length(rc.chi2nsym)) {
    if (length(rc.chi2sym)) {
      warning("both chi2_sym and chi2_nsym results in ", file)
    }
    return(rc.chi2nsym)
  }
  return(rc.chi2sym)
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
    prev.options = options()
    fields = unlist(strsplit(lines[li], "\\s+", perl=TRUE))
    options(warn=-1)
    values = as.numeric(fields[-1])
    options(prev.options)
    if (length(values) == 0 || any(is.na(values))) {
      break
    }
    data = rbind(data, values)
    rows = c(rows, fields[1])
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
    warning("cannot find '", pattern, "' in ", file)
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

  lines = suppressWarnings(try(get.file.lines(log.file("average_alucomb.log")), silent=TRUE))
  if (inherits(lines, "try-error")) {
    warning("Cannot open / read file ", file)
    return(ol)
  }
  
  rc = get.alucomb.section(lines, "^#+\\s+S-factors accounting", file, 1)
  if (is.null(rc)) {
    return(ol)
  }
  li = 1 + rc$lines
  ol$chisq = rc$data["updated", "chisq"]
  ol$dof = rc$data["updated", "dof"]
  
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
  ol$err = data["upd.error",]
  ol$sfact = data["S-factor",]
  ol$chisq.all = data["chisq",]
  ol$dof.all = data["dof",]

  return(ol)
}

## ////////////////////////////////////////////////////////////////////////////
## code

file = "average.input"

##++ aluplot = function(file = "") {

rc = alucomb.read(file)
measurements = rc$measurements
combination = rc$combination

log.dir = file.path("../../../Data", sub("^.*/([^/]*/[^/]*/[^/]*)$", "\\1", getwd()))
log.file = function(fname) {file.path(log.dir, fname)}

##-- get quantities measured by each experiment
meas.names = names(measurements)

##++ gets names in Combos MEASUREMENT card, sometimes incorrect
meas.quantities.raw = unlist(lapply(measurements, function(x) names(x$value)))
##++ get accurate names for quantities measured by experiments
meas.quantities = sub("[.][^.]+$", "", sub("^[^.]+[.]", "", names(measurements)))

##-- quantities we want to determine from measurements
## quant.names = combination$quantities
quant.names = unique(meas.quantities)

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

plot.data = list()
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
  for (meas in meas.names[quant == meas.quantities]) {
    exp.meas = list()
    value = measurements[[meas]]$value
    stat = measurements[[meas]]$stat
    syst = measurements[[meas]]$syst
    bibitem = measurements[[meas]]$bibitem

    error = sqrt(sum(c(stat), syst)^2)
    x.mins = c(x.mins, value - error)
    x.maxs = c(x.maxs, value + error)
    x.values = c(x.values, value)

    exp.meas$value = c(value, +stat, -stat, +syst, -syst)
    exp.meas$bibitem = paste(c(bibitem[1], bibitem[-(1:3)]), collapse=" ")

    if (!is.null(meas.asymm[[meas]])) {
      ##-- override measurement with info on asymmetric errors
      cat("ovverride", meas, "\n")
      cat("  from", sprintf("%10.4g ", exp.meas$value), "\n")
      cat("    to", sprintf("%10.4g ", meas.asymm[[meas]]), "\n")
      exp.meas$value = meas.asymm[[meas]]
    }

    all.exp.meas[[meas]] = exp.meas
  }

  ##
  ## get average info either from alucomb or Combos
  ##
  value = plot.data$hfag$val[quant]
  error = plot.data$hfag$err[quant]
  if (is.null(value)) {
    warning("could not find HFAG average for ", quant)
  } else {
    chisq = plot.data$hfag$chisq.all[quant]
    dof = plot.data$hfag$dof.all[quant]
    scale.factor = plot.data$hfag$sfact[quant]
    if (!is.null(chisq)) {
      conf.lev = pchisq(chisq, df=dof, lower.tail=FALSE)
      conf.lev.plot = conf.lev
      ##-- if using alucomb, set S-factor if larger than one
      if (scale.factor > 1) {
        conf.lev.plot = -scale.factor
      }
    } else {
      ##-- alucomb data not present, use combos common data
      conf.lev = plot.data$hfag$CL
      conf.lev.plot = conf.lev
      scale.factor = plot.data$hfag$scale
      ##++ combos data rule for setting scale, ask Swagato
      if (conf.lev < 1e-3) {
        ##++ remove following "if" when scale.factor for chi_sym fixed
        if (conf.lev != -scale.factor) {
          error = error*scale.factor
          conf.lev.plot = -scale.factor
        }
      }
    }
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
  x.padding.left = 0.05+2
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
  cat("* ", sprintf("%10.4g ", c(quant.data$xmin/x.units, quant.data$xmax/x.units)),
      as.character(quant.data$precision), " ", label.root(quant), units.label, "\n", file=fh, sep="")
  cat("# next lines are average, error, CL (or -ScaleFactor) for HFAG Averages\n", file=fh)

  ##-- print HFAG averages
  comb = plot.data[[quant]]$hfag$avg
  if (!is.null(comb)) {
    cat("& ", sprintf("%-10.4g ", c(comb/x.units, plot.data[[quant]]$hfag$CL.plot)), "HFAG Average\n", file=fh, sep="")
    cat("# next lines are average, error, Scale Factor for PDG Averages; Scale==0 means none quoted\n", file=fh)
  }

  ##-- PDG average
  comb = plot.data[[quant]]$pdg$avg
  if (!is.null(comb)) {
    cat("% ", sprintf("%-10.4g ", c(comb/x.units, plot.data[[quant]]$pdg$scale)), "PDG'09 Average\n", file=fh, sep="")
  } else {
    cat("% ", sprintf("%-10.4g ", c((2*quant.data$xmax - quant.data$xmin),0,0)/x.units), ">>> NOT FOUND <<< PDG'09 Average\n", file=fh, sep="")
  }
  
  ##-- measurements
  cat("# next lines are measurement, stat pos-error, neg-error[with negative sign], syst pos-error, neg-error[with negative sign], experiment-name\n", file=fh)
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
    cat("  ", sprintf("%-10.4g ", exp$value/x.units), bibitem, "\n", file=fh, sep="")
  }
  close(fh)
  cat("file", fname, "created\n")
}

##++} ## aluplot()
