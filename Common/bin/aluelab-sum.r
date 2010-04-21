#!/usr/bin/env Rscript

## ////////////////////////////////////////////////////////////////////////////
##
## aluelab-sum.r
##
## - compute sum of averages in current directory
##
## ////////////////////////////////////////////////////////////////////////////

source("../../../Common/bin/aluelab.r")

##-- collect data for current directory
aeb.collect.data(basename(getwd()))

##-- compute B(tau -> hhh nu) as linear combination of averaged quantities
meas.sum.lc = meas.val*0 + rep(1, length(meas.val))
meas.sum = aeb.linear.comb.glob(meas.sum.lc)

cat("summing following measurements:\n")
show(rbind(val=meas.val, err=meas.err))
cat("result:\n")
show(rbind(sum=meas.sum))
