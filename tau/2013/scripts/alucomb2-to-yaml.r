#!/usr/bin/env Rscript

'usage: alocomb2-to-yaml.r <R data file> [<YAML file>]

options:

' -> doc

## /////////////////////////////////////////////////////////////////////////////
##
##	alucomb2-to-yaml.r
##
## - read R data file produced by alocomb2.r and convert it to YAML
##
## /////////////////////////////////////////////////////////////////////////////

suppressPackageStartupMessages(require(docopt))
require(yaml, quietly=TRUE)
require(stringr, quietly=TRUE)

## /////////////////////////////////////////////////////////////////////////////
##
##	functions
##
## /////////////////////////////////////////////////////////////////////////////

##--- convert list to data.frame
list.to.df = function(lst, keyname=NA, stringsAsFactors=FALSE) {
  rc = do.call(rbind, lapply(lst, function(x) as.data.frame(x, stringsAsFactors=stringsAsFactors)))
  if (!is.na(keyname)) {
    rc[[keyname]] = as.character(names(lst))
  }
  rc
}

##--- convert data.frame to list
df.to.list = function(df, keyfield=NA, keyname="key") {
  if (!is.na(keyfield)) {
    rc = split(df[, -keyfield], seq(length.out=nrow(df)))
    names(rc) = df[, keyfield]
  } else {
    rc = split(df, seq(length.out=nrow(df)))
  }
  rc
}

load.in.list <- function(.file.name) { load(.file.name); as.list(environment()) }

## /////////////////////////////////////////////////////////////////////////////
##
##	code
##
## /////////////////////////////////////////////////////////////////////////////

opts <- docopt(doc)

ifname = opts[["R data file"]]
ofname = opts[["YAML file"]]
if (is.null(ofname)) {
  ofname = paste0(sub("[.][^.]*$", "", ifname, perl=TRUE), ".yaml")
}

data.list = load.in.list(ifname)

data.list$combination$constr.all.expr = NULL
data.list$combination$constr.all.expr.input = NULL

meas.val = data.list$meas.val
meas.err = data.list$meas.err
meas.corr = data.list$meas.corr

quant.val = data.list$quant.val
quant.err = data.list$quant.err
quant.corr = data.list$quant.corr

meas.names = names(meas.val)
meas.quant = sapply(data.list$measurements[meas.names], function(meas) meas$quant)

quant.names = names(quant.val)
quant.node = sapply(data.list$combination$quantities, function(quant) {ifelse(is.null(quant$node), "", quant$node)})
quant.node = quant.node[quant.names]

cat(file=ofname,
    as.yaml(list(
      meas.names = meas.names,
      meas.quant = meas.quant,
      meas.val = meas.val,
      meas.err = meas.err,
      meas.corr = meas.corr,
      quant.names = quant.names,
      quant.node = quant.node,
      quant.val = quant.val,
      quant.err = quant.err,
      quant.corr = quant.corr
    ),
    column.major = FALSE))

cat("file '", ofname, "' created\n", sep="")
