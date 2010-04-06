#!/usr/bin/env Rscript

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
## representation of Combos input card files data
##

##
## a "measurement" is stored in a list with fields
## - tag = "character",
## - value = "numeric",
## - stat = "numeric",
## - syst = "numeric",
## - bibitem = "character",
## - params = "list",
## - syst.terms = "numeric",
## - corr.terms = "list",
## - corr.terms.tot = "list"
## 

##
## a "combination" is stored as a list with fields:
##
## - tag = "character",
## - bibitem = "character",
## - quantities = "character",
## - params = "list",
## - meas.lin.combs = "list"
##

## ////////////////////////////////////////////////////////////////////////////
## code

##
## alucomb.read
##

alucomb.read = function(file = "") {

if (match.nocase("^\\s*$", file)) {
  stop("alucomb: please provide as argument the card input file\n")
}

flag.replace.corr = FALSE
flag.build.delta = FALSE

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
combination = list()

flag.in.meas = FALSE
flag.in.data = FALSE
flag.in.params = FALSE
flag.in.combine = FALSE
flag.in.sumofmeas = FALSE

for (line in lines) {
  ## cat(line,"\n")
  if (regexpr("^\\s*$", line, perl=TRUE) != -1 ||
      regexpr("^[*#;]", line, perl=TRUE) != -1) {
    next
  }
  line = gsub("\\s*!.*", "", line, perl=TRUE)
  fields = unlist(strsplit(line, "\\s+", perl=TRUE))
  if (match.nocase("^BEGIN$", fields[1])) {
    if (flag.in.meas) {
      cat("error, BEGIN keyword inside definition of", meas$tag, "\n")
    } else {
      ##--- init measurement list
      meas = list()
      meas$tag = paste(as.character(fields[2:4]), collapse=".")
      meas$bibitem = as.character(fields[-1])
      meas$value = numeric()
      meas$stat = numeric()
      meas$syst = numeric()
      meas$params = list()
      meas$syst.terms = numeric()
      meas$corr.terms = list()
      meas$corr.terms.tot = list()

      meas.labels = character()
      data.labels = character()
      data.values = numeric()
      sumofmeas.values = numeric()
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
      sumofmeas.values = numeric()
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
      ##-- get measurement values from following DATA values
      meas.values = numeric(length(meas.labels))
      names(meas.values) = meas.labels
      for (meas.label in meas.labels) {
        if (is.na(data.values[meas.label])) {
          stop("missing MEASUREMENT data for ", meas.label, "\n")
        }
        meas.values[meas.label] = data.values[meas.label]
      }

      ##-- edit measurement value, stat, syst labels
      labels = meas.labels
      labels = gsub("^m_", "", labels, ignore.case=TRUE)
      labels = gsub("statistical", "stat", labels, ignore.case=TRUE)
      labels = gsub("systematic", "syst", labels, ignore.case=TRUE)
      names(meas.values) = labels

      if (meas$bibitem[2] != meas.labels[1]) {
        ##-- when combining multiple quantities, "method" should same as the measured quantity name
        cat("warning: measurement method '", meas$bibitem[2], "' does not match value '", data.labels[1], "'\n", sep="")
      }
      ##-- deal with stat & syst errors expressed as percentage of value
      patt.perc = "([[:alnum:]]+[^[:alnum:]]*)(%)([^[:alnum:]]*)$"
      for (i in 2:length(meas.values)) {
        if (regexpr(patt.perc, meas.labels[i]) != -1) {
          meas.labels[i] = gsub(patt.perc, "\\1\\3", meas.labels[i])
          meas.values[i] = meas.value[1] * meas.values[i] /100
        }
      }
      meas$value = meas.values[1]
      meas$stat = meas.values[2]
      meas$syst = meas.values[3]
      meas$syst.terms = data.values[!data.labels %in% meas.labels]
      measurements[[meas$tag]] = meas
    } else {
      ##
      ## combination cards
      ##
      combination$tag = meas$tag
      combination$bibitem = meas$bibitem
      combination$quantities = meas.labels
      combination$params = meas$params
      combination$meas.lin.combs = measlincombs.list
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
    meas$corr.terms = c(meas$corr.terms,corr)
    next
  }
  if (match.nocase("^ERROR_CORR_WITH$", fields[1])) {
    corr = as.numeric(fields[5])
    names(corr) =  paste(as.character(fields[2:4]), collapse=".")
    meas$corr.terms.tot = c(meas$corr.terms.tot, corr)
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
    meas$params = c(meas$params, list(unlist(params.data[2:4])))
    names(meas$params)[length(meas$params)] = params.data[[1]]
    next
  }
  if (match.nocase("^SUMOFMEAS$", fields[1]) ||
      (flag.in.sumofmeas && match.nocase("^\\s*$", fields[1]))) {
    flag.in.sumofmeas = TRUE
    sumofmeas.values = c(sumofmeas.values, as.character(fields[-1]))
  }
}

list(combination=combination, measurements=measurements)
} ##-- end of alubomb.read
