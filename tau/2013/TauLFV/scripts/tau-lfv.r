#!/usr/bin/env Rscript

## /////////////////////////////////////////////////////////////////////////////
##
##	tau-lfv.r
##
## /////////////////////////////////////////////////////////////////////////////

require(yaml, quietly=TRUE)
require(stringr, quietly=TRUE)

## /////////////////////////////////////////////////////////////////////////////
##
##	functions
##
## /////////////////////////////////////////////////////////////////////////////

##--- convert data.frame to list
df.to.list = function(df, keyfield=NA, keyname="key") {
  if (!is.na(keyfield)) {
    rc = split(df[, -keyfield], seq(length.out=nrow(df)))
    names(rc) = df[, keyfield]
  } else {
    rc = split(df, seq(length.out=nrow(df)))
  }
  rc
}

##--- convert list to data.frame
list.to.df = function(lst, keyname=NA, stringsAsFactor=FALSE) {
  rc = do.call(rbind, lapply(in.list, function(x) as.data.frame(x, stringsAsFactors=stringsAsFactors)))
  if (!is.na(keyname)) {
    rc[[keyname]] = as.character(names(lst))
  }
  rc
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

##
## memo on how the data was prepared
##
## labels.df = read.csv("a.txt", stringsAsFactors=FALSE)
## lfv.df = read.delim(ifname, sep="", stringsAsFactors=FALSE)
## 
## data.lst = list(gamma.labels=as.list(gamma.labels), lfv.info.df=lfv.df)
## 
## ao = as.yaml(data.lst, column.major = FALSE)
## cat(ao, file="a.yaml")
##

fname = "tau-lfv-data"
ifname = paste(fname, "yaml", sep=".")
ofname = paste(fname, "tex", sep=".")

lfv.data = yaml.load_file(ifname)

out.txt = character(0)

rc = mapply(function(gamma.num, label) {
  label = paste0("g", gamma.num)
  alurep.tex.cmd.short(paste0(label, ".texlabel"), label)
}, names(lfv.data$gamma.labels), lfv.data$gamma.labels)
out.txt = c(out.txt, rc)

rc = lapply(lfv.data$lfv.info, function(x) {
  rc = character(0)
  gamma.num = x$gamma
  exp = tolower(as.vector(str_match(x$result, "([^_]*)_"))[2])
  label = paste("g", gamma.num, ".", exp, sep="")
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
out.txt = c(out.txt, rc)

cat(paste(out.txt, collapse="\n"), "\n", file=ofname)
cat("file '", ofname, "' created\n", sep="")
