#!/usr/bin/env Rscript

##
## aluelab-dump.r [flags] [<.rdata file>]
## dump some saved info from the fit
##

## ////////////////////////////////////////////////////////////////////////////
## functions

## ////////////////////////////////////////////////////////////////////////////
## code

aluelab.dump = function(args) {
  if (length(args) > 0) {
    file.name = args[1]
  } else {
    file.name = "average.rdata"
  }

  ##--- get alucomb results and data
  load(file.name)

  ## print(quant.corr)

  corr.flag = abs(quant.corr) >= 1.005
  gamma.flag = apply(corr.flag, 2, function(x) {sum(x)}) != 0
  gamma.flag = gamma.flag[gamma.flag]
  gamma.flag.names = names(gamma.flag)

  print(quant.corr[gamma.flag.names, gamma.flag.names])
}

args = commandArgs(TRUE)
aluelab.dump(args)
