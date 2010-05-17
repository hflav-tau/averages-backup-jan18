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
## meas.sum.lc = meas.val*0 + rep(1, length(meas.val))
meas.sum = aeb.linear.comb.glob(c(
  Gamma3=1,   Gamma5=1,   Gamma9=1,   Gamma10=1,  Gamma14=1,  Gamma16=1, 
  Gamma20=1,  Gamma23=1,  Gamma27=1,  Gamma28=1,  Gamma30=1,  Gamma35=1, 
  Gamma37=1,  Gamma40=1,  Gamma42=1,  Gamma47=1,  Gamma48=1,  Gamma62=1, 
  Gamma70=1,  Gamma77=1,  Gamma78=1,  Gamma85=1,  Gamma89=1,  Gamma93=1, 
  Gamma94=1,  Gamma103=1, Gamma104=1, Gamma126=1, Gamma128=1, Gamma130=1,
  Gamma132=1, Gamma150=1, Gamma152=1)) 

cat("aluelab-taufit.r: sum all basic modes:\n")
show(rbind(val=meas.val, err=meas.err))

cat("aluelab-taufit.r: result of sum of all basic modes:\n")
show(rbind(sum=meas.sum))
