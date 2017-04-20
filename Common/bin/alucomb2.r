#!/usr/bin/env Rscript

##
## alucomb2.r
##
## uses alucomb2-fit.r code to average measurements
##

getScriptPath <- function() {
  cmd.args <- commandArgs()
  m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
  script.dir <- dirname(regmatches(cmd.args, m))
  if (length(script.dir) == 0) {
    script.dir = try({dirname(sys.frame(1)$ofile)}, silent=TRUE)
    if (inherits(script.dir, "try-error")) {
      stop("cannot determine script dir")
    }
  }
  return(script.dir[1])
}

source(file.path(getScriptPath(), "alucomb2-utils.r"))
source(file.path(getScriptPath(), "alucomb2-fit.r"))

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
