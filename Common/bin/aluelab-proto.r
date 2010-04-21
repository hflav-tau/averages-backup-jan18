#!/usr/bin/env Rscript

## ////////////////////////////////////////////////////////////////////////////
##
## aluelab-proto.r
##
## - example on how to compute a linear combination of averages
##
## ////////////////////////////////////////////////////////////////////////////

source("../../../Common/bin/aluelab.r")

##-- collect data in global variables
aeb.collect.data(c("TauTo1Prong",
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
                   ))
## "TauToF1Nu",
## "TauToKmPizKstarzNu",
## "TauToPimPhiNu",
aeb.meas.add.single("TauToPimPhiNu", 3.42e-5, 0.604152e-5)

## show(rbind(value=meas.val, error=meas.err))
## show(meas.corr)
## show(meas.cov)

##-- compute B(tau -> hhh nu) as linear combination of averaged quantities
meas.hhh = aeb.linear.comb.glob(c(PimPimPipNu=1, PimKmPipNu=1, PimKmKpNu=1, KmKmKpNu=1))
show(rbind(HmHmHpNu=meas.hhh))
