#!/usr/bin/env Rscript

## ////////////////////////////////////////////////////////////////////////////
##
## aluelab-show.r
##
## - print some alucomb data as requested by Swagato
##
## ////////////////////////////////////////////////////////////////////////////

library(methods)

args <- commandArgs(TRUE)

if (any(args == "-h")) {
  ##--- use S-factors
  cat("aluelab_show.r\n")
  args = args[args != "-h"]
}

if (length(args) > 0) {
  file.name = args[1]
} else {
  file.name = "average_alucomb.rdata"
}

##--- get alucomb results and data
load(file.name)

options.save = options()
options(width=10)

meas.names = names(meas.val)
meas.sel.babar.belle = regexpr("babar|belle", meas.names, perl=TRUE, ignore.case=TRUE) != -1

cat("\nbabar belle measurements\n\n")
show(rbind(val=meas.val[meas.sel.babar.belle],
          err=meas.err[meas.sel.babar.belle]))
cat("\nbabar belle correlation\n\n")
show(meas.corr[meas.sel.babar.belle, meas.sel.babar.belle])

options(options.save)
