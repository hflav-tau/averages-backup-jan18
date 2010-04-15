#!/usr/bin/env Rscript

## ////////////////////////////////////////////////////////////////////////////
##
## aluelab.r
##
## ////////////////////////////////////////////////////////////////////////////

library(methods)

source("../../../Common/bin/alu-utils.r")

## ////////////////////////////////////////////////////////////////////////////
## definitions

log.dir = file.path("../../../Data", sub("^.*/([^/]*/[^/]*/[^/]*)$", "\\1", getwd()))
log.dir.base = dirname(log.dir)
cur.dir = basename(log.dir)
log.file = function(fname) {file.path(log.dir, fname)}

conf.lev.min = pnorm(3, mean=0, sd=1, lower.tail=FALSE)

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
## get numeric data from alucomb log file
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
## get a section of numeric data from alucomb log
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
  rc$lines = rc$lines + li -1
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
  ol$chisq = rc$data["original", "chisq"]
  ol$dof = rc$data["original", "dof"]
  
  rc.av = get.alucomb.section(lines[-(1:(li-1))], "^Averaged quantities", file)
  if (is.null(rc.av)) {
    return(ol)
  }
  li = li + rc.av$lines
  data = rc.av$data
  
  rc.corr = get.alucomb.section(lines[-(1:(li-1))], "^correlation", file)
  if (!is.null(rc.corr)) {
    li = li + rc.corr$lines
    ol$corr = rc.corr$data
  }
  
  rc.corr2 = get.alucomb.section(lines[-(1:(li-1))], "^correlation, S-factor inflated", file)
  if (!is.null(rc.corr2)) {
    li = li + rc.corr2$lines
    ol$corr2 = rc.corr2$data
  }
  
  rc.extra = get.alucomb.section(lines[-(1:(li-1))], "^Non-averaged measurement types", file)
  if (!is.null(rc.extra)) {
    li = li + rc.extra$lines
    data = cbind(data, rc.extra$data)
  }
  
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

meas.val = numeric(0)
meas.err = numeric(0)
meas.corr = matrix(0,0,0)

meas.add = function(label, val, err) {
  names(val) = label
  names(err) = label
  meas.add.alucomb(val, err)
}

##
## add a vector of measurements and their correnation to the list
##
meas.add.alucomb = function(add.val, add.err, add.corr=NULL) {
  ##
  ## if a correlation matrix is found, use it to determine which quantities where averaged
  ## and which quantities are non-averaged combined quantities
  ## retain only averaged quantities and their correlation
  ##
  if (!is.null(add.corr)) {
    quant.names.averaged = rownames(add.corr)
  } else {
    quant.names.averaged = names(add.val)
    add.corr = as.matrix(diag.m(rep(1, length(add.val))))
    rownames(add.corr) = names(add.val)
    colnames(add.corr) = names(add.val)
  }

  ##-- assemble correlation matrix
  corr.right = matrix(0, dim(meas.corr)[1], dim(add.corr)[2])
  colnames(corr.right) = colnames(add.corr)
  corr.top = cbind(meas.corr, corr.right)
  corr.left = matrix(0, dim(add.corr)[1], dim(meas.corr)[2])
  rownames(corr.left) = rownames(add.corr)
  corr.bottom = cbind(corr.left, add.corr)
  meas.corr <<- rbind(corr.top, corr.bottom)

  ##-- assemble averaged quantities
  quant.averaged.sel = names(add.val) %in% quant.names.averaged
  meas.val <<- c(meas.val, add.val[quant.averaged.sel])
  meas.err <<- c(meas.err, add.err[quant.averaged.sel])
}

for (dir in c("TauTo1Prong",
              "TauTo3Prongs",
              ## "TauToF1Nu",
              "TauToHmHmHpNu",
              "TauToKmEtaNu",
              "TauToKmKstarzNu",
              "TauToKmPhiNu",
              "TauToKmPizEtaNu",
              ## "TauToKmPizKstarzNu",
              "TauToKstarmEtaNu",
              "TauToPimF1Nu",
              "TauToPimKzbEtaNu",
              "TauToPimKzbNu",
              ## "TauToPimPhiNu",
              "TauToPimPimPipEtaNu",
              "TauToPimPizEtaNu",
              "TauToPimPizKzbNu"
              )) {
  ##-- get alucomb results
  rc = get.alucomb(file.path(log.dir.base, dir, "average_alucomb.log"))
  if (length(rc) == 0) {
    cat("error: cannot read alucomb log file\n  ", file.path(log.dir.base, dir, "average_alucomb.log"), "\n", sep="")
    ## cat("make -C", dir, " update_alucomb\n")
    ## warning("Cannot read alucomb log file\n  ", file.path(log.dir.base, dir, "average_alucomb.log"))
    next
  }

  ##-- default correlation matrix (NULL if not read)
  corr.current = rc$corr
  conf.lev = pchisq(rc$chisq, df=rc$dof, lower.tail=FALSE)
  if (conf.lev < conf.lev.min) {
    ##-- inflate errors if CL < CL_min
    rc$err = rc$err * rc$sfact
    ##-- correlation after S-factor error inflation
    corr.current = rc$corr2
  }

  meas.add.alucomb(rc$val, rc$err, corr.current)
}

## "TauToF1Nu",
## "TauToKmPizKstarzNu",
## "TauToPimPhiNu",
meas.add("TauToPimPhiNu", 3.42e-5, 0.604152e-5)

##-- names of all measurements
meas.names = names(meas.val)

##-- variance matrix
meas.cov = meas.corr * (meas.err %o% meas.err)

## show(rbind(value=meas.val, error=meas.err))
## show(meas.corr)
## show(meas.cov)

##
## compute linear combinations
##
linear.comb = function(lc) {
  meas.lc = meas.names %in% names(lc)
  lc.val = meas.val[meas.lc] %*% lc
  lc.err = sqrt(drop(lc %*% meas.cov[meas.lc, meas.lc] %*% lc))
  return(c(val=lc.val, err=lc.err))
}

##-- compute B(tau -> hhh nu) as linear combination of averaged quantities
meas.hhh = linear.comb(c(PimPimPipNu=1, PimKmPipNu=1, PimKmKpNu=1, KmKmKpNu=1))
show(rbind(HmHmHpNu=meas.hhh))

##++} ## aluelab()
