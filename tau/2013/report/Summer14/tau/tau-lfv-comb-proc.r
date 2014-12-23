#!/usr/bin/env Rscript

require(stringr, quietly=TRUE)

## /////////////////////////////////////////////////////////////////////////////
##
##	functions
##
## /////////////////////////////////////////////////////////////////////////////

##--- texify
texify = function(text) {
  text = gsub("([aeiou])(['`])", "\\\\\\2\\1", rlines, perl=TRUE, ignore.case=TRUE)
  text = gsub("([_^])", "\\\\verb|\\1|", rlines, perl=TRUE, ignore.case=TRUE)
  return(text)
}

df.to.list <- function( df ) {
  split(df, rownames(df))
}

##--- return latex command def with specified multi-line body
alurep.tex.cmd.nv = function(cmd, body) {
  paste("\\htdef{", cmd, "}{%\n", body, "}%", sep="")
}
alurep.tex.cmd = Vectorize(alurep.tex.cmd.nv)

##--- return latex command def with specified one-line body
alurep.tex.cmd.short.nv = function(cmd, body) {
  paste("\\htdef{", cmd, "}{", body, "}%", sep="")
}
alurep.tex.cmd.short = Vectorize(alurep.tex.cmd.short.nv)

## /////////////////////////////////////////////////////////////////////////////
##
##	code
##
## /////////////////////////////////////////////////////////////////////////////

file = "tau-lfv-comb-data.txt"
lfv.df = read.delim("tau-lfv-comb-data.txt", sep="")

rc = lapply(df.to.list(lfv.df), function(x) {
  gamma = x$Gamma
  rc = character(0)
  exp = tolower(as.vector(str_match(x$result, "([^_]*)_"))[2])
  label = paste("g", gamma, ".", exp, sep="")
  rc = c(rc, alurep.tex.cmd.short(paste(label, ".lumi", sep=""), x$lumi))
  rc = c(rc, alurep.tex.cmd.short(paste(label, ".xsec", sep=""), x$cross.section))
  rc = c(rc, alurep.tex.cmd.short(paste(label, ".eff", sep=""),
    sprintf("$%.2f \\pm %.2f$", x$efficiency*100, x$efficiency.error*100)))
  rc = c(rc, alurep.tex.cmd.short(paste(label, ".bkg", sep=""),
    sprintf("$%.2f \\pm %.2f$", x$bkg, x$bkg.error)))
  rc = c(rc, alurep.tex.cmd.short(paste(label, ".evs", sep=""), x$observed.events))

  row = paste("\\htuse{", label, ".", c("lumi", "xsec", "eff", "bkg", "evs"), "}", sep="")
  rc = c(rc, alurep.tex.cmd.short(paste(label, ".row", sep=""), paste(row, collapse=" & ")))

  rc = paste(rc, collapse="\n")
  rc
})

fname = "tau-lfv-comb-data.tex"
cat(paste(rc, collapse="\n"), "\n", file=fname)
cat("file '", fname, "' created\n", sep="")
