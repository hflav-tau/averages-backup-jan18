#!/usr/bin/env Rscript

'usage: alucomb2_cmp.r <file1> <file2>

options:

' -> doc

suppressPackageStartupMessages(require(docopt))

##
## functions
##

load.in.list <- function(.file.name) { load(.file.name); as.list(environment()) }

##
## main code
##

opts <- docopt(doc)

f1 = load.in.list(opts$file1)
f2 = load.in.list(opts$file2)

quant.name.f1 = names(f1$quant.val)
quant.name.f2 = names(f2$quant.val)
quant.name.common = intersect(quant.name.f1, quant.name.f2)
in.f1.not.f2 = setdiff(quant.name.f1, quant.name.f2)
if (length(in.f1.not.f2) > 0) {
  cat("--- fitted values in file1 but not in file2\n")
  cat(paste(in.f1.not.f2, collapse=" "))
  cat("\n--- end\n")
}
in.f2.not.f1 = setdiff(quant.name.f2, quant.name.f1)
if (length(in.f2.not.f1) > 0) {
  cat("--- fitted values in file2 but not in file1\n")
  cat(paste(in.f2.not.f1, collapse=" "))
  cat("\n--- end\n")
}

f1$quant.val = f1$quant.val[quant.name.common]
f2$quant.val = f2$quant.val[quant.name.common]
f1$quant.err = f1$quant.err[quant.name.common]
f2$quant.err = f2$quant.err[quant.name.common]

val.delta.sigma = sqrt((f1$quant.err^2 + f2$quant.err^2)/2)

cmp.df = data.frame(
  val.delta.perc = 100 * (f2$quant.val - f1$quant.val) / val.delta.sigma,
  err.delta.perc = 100 * (f2$quant.err - f1$quant.err) / val.delta.sigma,
  val.delta.sigma = val.delta.sigma,
  v1 = f1$quant.val,
  v2 = f2$quant.val,
  s.v1 = f1$quant.err,
  s.v2 = f2$quant.err,
  descr = sapply(quant.name.common, function(qn) {f1$combination$quantities[[qn]]$descr}),
  gamma = quant.name.common)

rm(val.delta.sigma)

val.delta.order = order(abs(cmp.df$val.delta.perc), decreasing=TRUE)
cmp.df = cmp.df[val.delta.order, ]

cat(paste0("compare v1=", opts$file1, " with v2=", opts$file2, "\n"))
cat("\n")
cat("(v2-v1)\n")
cat("--------        v1         v2   s(v1-v2)      [ s(v2-v1) = sqrt[(s_v1^2 + s_v2^2)/2] ]\n")
cat("s(v2-v1)\n")
cat("  (%)\n")
cat("\n")
rc = lapply(split(cmp.df, 1:nrow(cmp.df)),
  function(cmp) {
    cat(sprintf(
      "%7.2f %#10.4g %#10.4g %#10.4g  %-16s %-s\n",
      cmp$val.delta.perc, cmp$v1, cmp$v2, cmp$val.delta.sigma,
      cmp$gamma, cmp$descr))
  })

err.delta.order = order(abs(cmp.df$err.delta.perc), decreasing=TRUE)
cmp.df = cmp.df[err.delta.order, ]

cat("\n")
cat("(s_v2-s_v1)\n")
cat("-----------   s_v1       s_v2   s(v1-v2)      [ s(v2-v1) = sqrt[(s_v1^2 + s_v2^2)/2] ]\n")
cat(" s(v2-v1)\n")
cat("  (%)\n")
cat("\n")
rc = lapply(split(cmp.df, 1:nrow(cmp.df)),
  function(cmp) {
    cat(sprintf(
      "%7.2f %#10.4g %#10.4g %#10.4g  %-16s %-s\n",
      cmp$err.delta.perc, cmp$s.v1, cmp$s.v2, cmp$val.delta.sigma,
      cmp$gamma, cmp$descr))
  })
