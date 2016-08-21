#!/usr/bin/env Rscript

## /////////////////////////////////////////////////////////////////////////////
##
##	tau-lfv.r
##
## - read tau LFV data file from .yaml file
## - output TeX file with definitions useful to report tau LFV data
##
## /////////////////////////////////////////////////////////////////////////////

require(yaml, quietly=TRUE)
require(stringr, quietly=TRUE)

## /////////////////////////////////////////////////////////////////////////////
##
##	functions
##
## /////////////////////////////////////////////////////////////////////////////

##--- convert list to data.frame
list.to.df = function(lst, keyname=NA, stringsAsFactors=FALSE) {
  rc = do.call(rbind, lapply(lst, function(x) as.data.frame(x, stringsAsFactors=stringsAsFactors)))
  if (!is.na(keyname)) {
    rc[[keyname]] = as.character(names(lst))
  }
  rc
}

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

tau.lfv.data = yaml.load_file("tau-lfv-data.yaml")
ofname = "tau-lfv-data.tex"

gamma.df = list.to.df(tau.lfv.data$gamma)
names(tau.lfv.data$gamma) = gamma.df$gamma

out.txt = character(0)

rc = mapply(function(gamma, descr) {
  key = paste0("g", gamma)
  alurep.tex.cmd.short(paste0(key, ".descr"), paste0("\\ensuremath{", descr, "}"))
}, gamma.df$gamma, gamma.df$descr)
out.txt = c(out.txt, rc)

if (FALSE) {
rc = lapply(tau.lfv.data$combs.extra, function(x) {
  rc = character(0)
  label = paste0("g", x$gamma, ".", x$ref)
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
}

tau.lfv.data$combs.extra = tau.lfv.data$combs.extra[
  order(sapply(tau.lfv.data$combs.extra, function(el) el$gamma))]

rc = lapply(tau.lfv.data$combs.extra, function(br) {
  gamma.info = tau.lfv.data$gamma[[as.character(br$gamma)]]
  rc = paste0(
    "\\htCombExtraLine%\n",
    "  {\\ensuremath{\\Gamma_{", br$gamma, "} = ", gamma.info$descr, "}}%\n",
    "  {", br$exp, "}%\n",
    "  {", br$ref, "}%\n",
    ## "  {", br$lumi, "}%\n",
    ## "  {", br$cross.section, "}%\n",
    "  {", br$num.tau, "}%\n",    
    "  {", sprintf("\\ensuremath{%.2f \\pm %.2f}", br$efficiency*100, br$efficiency.error*100), "}%\n",
    "  {", sprintf("\\ensuremath{%.2f \\pm %.2f}", br$bkg, br$bkg.error), "}%\n",
    "  {", br$observed.events, "}"
  )
})

rc = alurep.tex.cmd(
  "CombExtraLines",
  paste0(rc, collapse="%\n")
)
out.txt = c(out.txt, rc)

tau.lfv.data$combs = tau.lfv.data$combs[
  order(sapply(tau.lfv.data$combs, function(el) el$gamma))]

descr.last = ""
type.last = ""
rc = lapply(tau.lfv.data$combs, function(br) {
  gamma.info = tau.lfv.data$gamma[[as.character(br$gamma)]]
  rc = character(0)
  descr = ifelse(gamma.info$descr != descr.last, gamma.info$descr, "")
  descr.last <<- gamma.info$descr
  type = ifelse(gamma.info$type != type.last, gamma.info$type, "")
  type.last <<- gamma.info$type

  if (type != "") {
    rc = c(rc, "\\hline")
  }
  
  rc = c(
    rc,
    paste0(
      "\\htCombLimitLine%\n",
      "  {\\ensuremath{\\Gamma_{", br$gamma, "} = ", gamma.info$descr, "}}%\n",
      "  {\\ensuremath{", type, "}}%\n",
      "  {\\ensuremath{", sub("^(.*)e([+-]*)0*(\\d*)$", "\\1 \\\\cdot 10^{\\2\\3}", sprintf("%.1e", br$limit)), "}}"
    ))
})

rc = alurep.tex.cmd(
  "CombLines",
  paste0(unlist(rc), collapse="%\n")
)
out.txt = c(out.txt, rc)

cat(paste(out.txt, collapse="\n"), "\n", file=ofname, sep="")
cat("file '", ofname, "' created\n", sep="")
