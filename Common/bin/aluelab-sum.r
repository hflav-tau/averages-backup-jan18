#!/usr/bin/env Rscript

## ////////////////////////////////////////////////////////////////////////////
##
## aluelab-sum.r
##
## - compute sum of averages in current directory
##
## ////////////////////////////////////////////////////////////////////////////

source("../../../Common/bin/aluelab.r")

args <- commandArgs(TRUE)
if (length(args) > 0) {
  items = args
} else {
  items = basename(getwd())
}

##-- collect data for current directory
aeb.collect.data(items)

##-- compute B(tau -> hhh nu) as linear combination of averaged quantities
quant.sum.lc = quant.val*0 + rep(1, length(quant.val))
quant.sum = aeb.linear.comb.glob(quant.sum.lc)

cat("summing following measurements:\n")
show(rbind(val=quant.val, err=quant.err))
cat("result:\n")
show(rbind(sum=quant.sum))
