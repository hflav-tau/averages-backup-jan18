#!/usr/bin/env Rscript

##
## aluplot
##

library(methods)

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

##++alucomb = function(file = "") {

file = "average.input"

if (match.nocase("^\\s*$", file)) {
  stop("alucomb: please provide as argument the card input file\n")
}

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
  if (regexpr("^\\s*$", line) != -1 ||
      regexpr("^[*#;]", line) != -1) {
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

if (FALSE) {
##-- include measurements that correspond to combination of quantities
##++ should probably check also that _all_ quantities are in the combination
for (meas in names(combination@meas.lin.combs)) {
  if (sum(combination@quantities %in% names(combination@meas.lin.combs[[meas]])) != 0) {
    cat("meas", meas, "included, as linear combination\n")
    meas.list[meas] = TRUE
  }
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
  } else if (syst.contribs > (1+1e-4)*measurements[[meas]]@syst) {
    cat("warning: sum of syst. terms slightly larger than syst. error\n  ",
        syst.contribs, " vs. ", measurements[[meas]]@syst, "\n", sep="")
  }
}

##-- get list of stat and syst errors
meas.stat.true = unlist(lapply(measurements, function(x) x@stat))
names(meas.stat.true) = meas.names
meas.syst.true = unlist(lapply(measurements, function(x) x@syst))
names(meas.syst.true) = meas.names

##
## quantities we want to determine from measurements
##
quant.names = combination@quantities

##++} ##-- end function alucomb

## args <- commandArgs(TRUE)
## if (length(args) > 0) alucomb(file = args[1]) 

for (quant in quant.names) {
  ## show(quant)
  ## show(meas.names[quant == meas.quantities])
}

label.root = function(str) {
  str.orig = str
  # str = sub("PimKzbNu", "B(#tau^{-} #rightarrow h^{-} #nu)", str)
  if (str == str.orig) {
    str = gsub("(^|[^A-Z])([A-Z][a-z])*b","\\1#bar{\\2}", str, perl=TRUE)
    str = gsub("Pi","#pi",str)
    str = gsub("Nu","#nu",str)
    str = gsub("m($|[#}A-Z])","^{-}\\1", str, perl=TRUE)
    str = gsub("p($|[#}A-Z])","^{+}\\1", str, perl=TRUE)
    str = gsub("z($|[#}A-Z])","^{0}\\1", str, perl=TRUE)
    str = paste("B(#tau^{-} #rightarrow ", str, ")", sep="")
  }
  return(str)
}

##
## enter here the PDG 2008 averages
## - value, error, scale factor
##
pdg.averages.text =
"PimPimPipNu  8.85e-2    0.13e-2       1
PimKmPipNu   0.280e-2    0.019e-2      2.1
PimKmKpNu    1.37e-3     0.06e-3       1.8
KmKmKpNu     1.58e-5     0.1769181e-5  1
PimKzbNu     0.852e-2    0.032e-2      1"

pdg.averages = list()
for (fields in strsplit(unlist(strsplit(pdg.averages.text, "\n")), "\\s+", perl=TRUE)) {
  pdg.averages[[fields[1]]] = as.numeric(fields[-1])
}

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
  if (regexpr(patt,average,perl=TRUE) == -1) return(list())
  fields = unlist(strsplit(sub(patt, "\\1;\\2;\\3", average, perl=TRUE), ";"))
  rc$scale = fields[1]
  rc$value = as.numeric(fields[2:3])
  return(rc)
}

get.results.chi2.sym = function(file) {
  rc = list()
  fh = pipe(paste(
    "awk 2>/dev/null '/CHI2_SYM: CHI2, NMEFF, CHI2\\/NDOF/ {flag++}; flag>1'",
    file, sep=" "))
  average = readLines(fh)
  close(fh)
  average = paste(average, collapse=";")
  patt = paste(".*CHI2_SYM\\s+:\\s+([^=\\s]+)\\s*=\\s*([.0-9eE+-]+)",
    "\\s*\\+-\\s*([.0-9eE+-]+)\\s*CL\\s*=\\s*([.0-9eE+-]+).*", sep="")
  if (regexpr(patt,average,perl=TRUE) == -1) return(list())
  fields = unlist(strsplit(sub(patt, "\\1;\\2;\\3;\\4", average, perl=TRUE), ";"))
  rc$CL = as.numeric(fields[4])
  ##++ need to recover also sqrt(chisq/dof)
  rc$scale = -rc$CL
  combs = list()
  label = sub("^M_", "", fields[1])
  combs[[label]] = as.numeric(fields[2:3])
  rc$combs = combs
  return(rc)
}

get.results.chi2.nsym = function(file) {
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
  patt = paste(".*CHI2_N_SYM\\s+:\\s+([^=\\s]+)\\s*=\\s*([.0-9eE+-]+)",
    "\\s*\\+-\\s*([.0-9eE+-]+)\\s*CL\\s*=\\s*([.0-9eE+-]+).*", sep="")
  patt2 = paste(".*\\s+([^=\\s]+)\\s*=\\s*([.0-9eE+-]+)",
    "\\s*\\+-\\s*([.0-9eE+-]+)\\s*.*", sep="")
  flag = FALSE  
  combs = list()
  for (line in average) {
    if (regexpr(patt, line) != -1) {
      fields = unlist(strsplit(sub(patt, "\\1;\\2;\\3;\\4", line, perl=TRUE), ";"))
      flag = TRUE
    } else if (flag) {
      if (regexpr(patt2, line) != -1) {
        fields = unlist(strsplit(sub(patt2, "\\1;\\2;\\3", line, perl=TRUE), ";"))
      }
    } else {
      flag = FALSE
    }
    if (flag) {
      quant = fields[1]
      quant = sub("^M_", "", quant)
      combs[[quant]] = as.numeric(fields[2:3])
    }
  }
  rc$combs = combs
  return(rc)
}

plot.data = list()
for (quant in quant.names) {
  plot.data[[quant]] = list()
}

rc = get.pdg.average("log/pdg_average.log")
if (length(rc) > 0) {
   plot.data$pdg = rc
}

rc = get.results.chi2.sym("log/average.log")
if (length(rc) > 0) {
  if (!is.null(plot.data[[quant.names[1]]]$hfag)) {
    cat("warning, double hfag def of", quant.names[1], "\n")
  }
  plot.data$hfag = rc
}

rc = get.results.chi2.nsym("log/average.log")
if (length(rc) > 0) {
  if (!is.null(plot.data[[quant.names[1]]]$hfag)) {
    cat("warning, double hfag def of", quant.names[1], "\n")
  }
  plot.data$hfag = rc
}

for (quant in quant.names) {
  x.mins = numeric(0)
  x.maxs = numeric(0)
  all.exp.meas = list()
  for (meas in meas.names[quant == meas.quantities]) {
    exp.meas = list()
    value = measurements[[meas]]@value
    stat = measurements[[meas]]@stat
    syst = measurements[[meas]]@syst
    bibitem = measurements[[meas]]@bibitem

    error = sqrt(sum(c(stat), syst)^2)
    x.mins = c(x.mins, value - error)
    x.maxs = c(x.maxs, value + error)

    exp.meas$value = c(value, -stat, +stat, -syst, +syst)
    exp.meas$bibitem = paste(c(bibitem[1], bibitem[-(1:3)]), collapse=" ")
    all.exp.meas[[meas]] = exp.meas
  }
  x.min = min(x.mins)
  x.max = max(x.maxs)
  x.mean = (x.min+x.max)/2
  x.units = exp(log(10)*round(log(x.mean)/log(10)))
  if (round(1/x.units) == 1000 || round(1/x.units) == 10) {
    x.units = 1/100
  }
  x.padding.left = 0.05+2
  x.padding.right = 0.05
  x.min = x.min - x.padding.left*(x.max-x.min)
  x.max = x.max + x.padding.right*(x.max-x.min)
  plot.data[[quant]]$xmin = x.min
  plot.data[[quant]]$xmax = x.max
  plot.data[[quant]]$expts = all.exp.meas
  plot.data[[quant]]$units = x.units
}

for (quant in quant.names) {
  fname = "plot.input"
  if (length(quant.names) > 1) {
    fname = paste("plot-", quant, ".input", sep="")
  }
  fh = file(fname, "w")
  quant.data = plot.data[[quant]]
  x.units = quant.data$units
  cat("# first line is xmin, xmax, units, title\n", file=fh)
  if (round(1/x.units) == 100) {
    units.label = "%"
  } else {
    units.label = sprintf("#times%.0e", x.units)
  }
  cat("* ", sprintf("%10.4g ", c(quant.data$xmin/x.units, quant.data$xmax/x.units)),
      sprintf("%10.0e ", x.units), label.root(quant), " [", units.label, "]\n", file=fh, sep="")

  ##-- print HFAG averages
  comb = plot.data$hfag$combs[[toupper(quant)]]
  if (is.null(comb)) {
    cat("warning: could not find HFAG average(s) for", quant, "\n")
  } else {
    conf.lev = plot.data$hfag$CL
    scale.factor = plot.data$hfag$scale
    if (conf.lev < 1e-3) {
      ##++ remove "if" when scale.factor for chi_sym fixed
      if (conf.lev != -scale.factor) {
        comb[2] = comb[2]*scale.factor
        conf.lev = -scale.factor
      }
    }
    cat("& ", sprintf("%10.4g ", c(comb/x.units, conf.lev)), "HFAG Average\n", file=fh, sep="")
    cat("# next lines are average, error, Scale Factor for HFAG Averages; Scale==0 means none quoted\n", file=fh)
  }

  ##-- PDG average
  pdgav = pdg.averages[[quant]]
  scale.factor = pdgav[3]
  if (scale.factor == 1) {
    scale.factor = 0
  }
  if (!is.null(pdgav)) {
    cat("% ", sprintf("%10.4g ", c(pdgav[1:2]/x.units, scale.factor)), "PDG'08 Average\n", file=fh, sep="")
  } else {
    cat("% ", sprintf("%10.4g ", c(0,0,0)/x.units), ">>> NOT FOUND <<< PDG'08 Average\n", file=fh, sep="")
  }
  
  ##-- measurements
  cat("# next lines are measurement, stat-error, syst-error, experiment-name\n", file=fh)
  for (exp in plot.data[[quant]]$expts) {
    bibitem = exp$bibitem
    ##-- capitalize 1st word
    bibitem = sub("^(\\S+)", "\\U\\1", bibitem, perl=TRUE)
    ##
    ## get the year from the "where" field in the Combos cards
    ## please end the "where" field either with (YYYY) or with YYYY
    ##
    bibitem = sub("^(\\S+)\\s+.*(\\(|)(\\d\\d\\d\\d)(\\)|)$", "\\U\\1 \\3", bibitem, perl=TRUE)
    ##--- change PDG-like years with possible letter in 4-digit years
    bibitem = sub("^(\\S+)\\s+.*\\s+([0-2]\\d)([A-Z]|)\\s*$", "\\U\\1 20\\2", bibitem, perl=TRUE)
    bibitem = sub("^(\\S+)\\s+.*\\s+([3-9]\\d)([A-Z]|)\\s*$", "\\U\\1 19\\2", bibitem, perl=TRUE)
    if (regexpr("\\d\\d\\d\\d$", bibitem, perl=TRUE) == -1) {
      cat("warning: could not find year in \"where\" Combos field for", quant, "\n")
      cat("  (", exp$bibitem, ")\n", sep="")
    }
    cat("  ", sprintf("%10.4g ", exp$value/x.units), bibitem, "\n", file=fh, sep="")
  }
  close(fh)
  cat("file", fname, "created\n")
}

