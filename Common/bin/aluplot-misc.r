#!/usr/bin/env Rscript

## ////////////////////////////////////////////////////////////////////////////
##
## aluplot-misc.r
##
## - demo code, get data from plot.input file
##   - print PDG average contained in plot.input of the current dir
##
## ////////////////////////////////////////////////////////////////////////////

library(methods)

source("../../../Common/bin/alu-utils.r")

## ////////////////////////////////////////////////////////////////////////////
## definitions

log.dir = file.path("../../../Data", sub("^.*/([^/]*/[^/]*/[^/]*)$", "\\1", getwd()))
log.dir.base = dirname(log.dir)
cur.dir = basename(log.dir)
quant.name = gsub("TauTo", "", cur.dir)

lines = suppressWarnings(try(readLines("plot.input"), silent=TRUE))
if (!inherits(lines, "try-error")) {
  for (line in lines) {
    fields = unlist(strsplit(line, "\\s+", perl=TRUE))
    if (fields[1] == "%") {
      value = as.numeric(fields[2])
      error = as.numeric(fields[3])
      s.factor =  as.numeric(fields[4])
      if (s.factor == 0) s.factor = 1
      cat(sprintf("%-15s %-16.6g%-16.6g%d\n", quant.name, value, error, s.factor))
    }
  }
}
