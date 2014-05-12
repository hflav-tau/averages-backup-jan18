#!/usr/bin/env Rscript

##
## alucomb2.r
##
## uses alucomb2-fit.code to average measurements
##

source("../../../Common/bin/alucomb2-fit.r")

## ////////////////////////////////////////////////////////////////////////////
## code

method="solnp"
method="alabama"
method="alucomb"
method="alucomb2"

args <- commandArgs(TRUE)
if (length(args) == 1 && exists("alucomb")) {
  alucomb(filename = args[1], method = method)
} else {
  cat("Usage: alucomb2.r <input file>\n")
}
