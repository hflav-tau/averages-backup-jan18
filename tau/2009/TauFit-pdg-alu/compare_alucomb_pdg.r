#!/usr/bin/env Rscript

## ////////////////////////////////////////////////////////////////////////////
##
## compare_alucomb_pdg.r
##
## - get averages from alucomb.r .rdata file
## - get PDG averages from pdginput/pdgread.cc
##
## ////////////////////////////////////////////////////////////////////////////

library(methods)

source("../../../Common/bin/alu-utils.r")

## ////////////////////////////////////////////////////////////////////////////
## definitions

## ////////////////////////////////////////////////////////////////////////////
## code

args <- commandArgs(TRUE)
if (length(args) > 0) {
  file.name = args[1]
} else {
  file.name = "average_alucomb.rdata"
}
file.name.pdg.swb = "pdginput/readpdg.cc"
file.name.pdg = "pdginput/s035-fit-no-babar-belle.data"

##--- get alucomb results and data
load(file.name)

##--- get pdg data
ll = readLines(file.name.pdg)
labels = tolower(unlist(strsplit(ll[2], "\\s+", perl=TRUE)))
labels[1] = "num"
labels[7] = paste(labels[4], labels[7], sep=".")
labels[6] = paste(labels[4], labels[6], sep=".")
labels[4] = paste(labels[4], labels[5], sep=".")
labels[10] = paste(labels[8], labels[10], sep=".")
labels[8] = paste(labels[8], labels[9], sep=".")
labels = labels[c(-5, -9)]
labels = c(labels, "descr")
ll.sel = grep("^\\*", ll, perl=TRUE)
ll = ll[ll.sel]
ll =
  gsub(paste("^.\\s*(\\d+)", paste(rep("\\s+(\\S+)", 7), sep="", collapse=""), "\\s+(.*\\S+)\\s*$", sep="", collapse=""),
       paste("\\", 1:9, sep="", collapse=";"), ll, perl=TRUE)
ll =
  gsub(paste("(", paste(rep("[^;\\s]*", 9), collapse=";"), ")",
             paste(rep("\\s+(\\S+)", 5), collapse=""), "\\s+(.*\\S+)\\s*$", sep="", collapse=""),
       paste("\\", 1:7, sep="", collapse=";"), ll, perl=TRUE)
con = textConnection(ll)
pdg = read.table(con, sep=";", col.names=labels)
close(con)
pdg$quant.name = paste("Gamma", pdg$gamma, sep="")

##--- get data from Swagato .cc file
ll = readLines(file.name.pdg.swb)
ll.sel = grep("push_back.*basefitvalue", ll, perl=TRUE, useBytes=TRUE)
ll = ll[ll.sel]

ll = gsub(paste("^.*basegamma.push_back.\\s*(\\S+)\\s*.;",
  ".*baseseed.ibase.\\s*=\\s*(\\S+)\\s*;",
  ".*basefitvalue.ibase.\\s*=\\s*(\\S+)\\s*;",
  ".*basefiterror.ibase.\\s*=\\s*(\\S+)\\s*;",
  ".*baserescalederror.ibase.\\s*=\\s*(\\S+)\\s*;",
  ".*basescalefactor.ibase.\\s*=\\s*(\\S+)\\s*;",
  ".*$",
  sep=""),
  "\\1;\\2;\\3;\\4;\\5;\\6", ll, perl=TRUE, useBytes=TRUE)

con = textConnection(ll)
pdg.swb = read.table(con, sep=";",
  col.names = c("quant.id", "quant.seed", "quant", "quant.err", "quant3.err", "sfact3.types"))
close(con)
pdg.swb$quant.name = paste("Gamma", pdg.swb$quant.id, sep="")

quant.names = names(quant)

pdg.select = quant.names %in% pdg.swb$quant.name
quant.names.pdg = quant.names[pdg.select]
pdg.swb.sel = match(quant.names.pdg, pdg.swb$quant.name)

## pdg.rows = sapply(unique(pdg$quant.name), function(x) which(x == pdg$quant.name)[1])
pdg.rows = match(quant.names.pdg, pdg$quant.name)
pdg.rows = pdg.rows[!is.na(pdg.rows)]

cat("             quant     err       err3      sf3  pdg       pdg.err   pdg.err3  pdg.sf3 Dquant%    Derr%       Derr2%    Dsf% \n")
rc = mapply(function(
  name,
  quant,
  quant.err,
  quant3.err,
  sfact3.types,
  pdg.quant,
  pdg.quant.err,
  pdg.quant3.err,
  pdg.sfact3.types
  ) {
  cat(sprintf("%12s %9.7f %9.7f %9.7f %4.2f %9.7f %9.7f %9.7f %4.2f %9.4f%% %9.4f%% %9.4f%% %9.4f%%\n",
              name,
              quant,
              quant.err,
              quant3.err,
              sfact3.types,
              pdg.quant,
              pdg.quant.err,
              pdg.quant3.err,
              pdg.sfact3.types,
              (quant/pdg.quant-1)*100,
              (quant.err/pdg.quant.err-1)*100,
              (quant3.err/pdg.quant3.err-1)*100,
              (sfact3.types/pdg.sfact3.types-1)*100
              ))
},
  quant.names.pdg,
  quant[quant.names.pdg],
  quant.err[quant.names.pdg],
  quant3.err[quant.names.pdg],
  sfact3.types[quant.names.pdg],
  ## pdg.quant$fit.value[pdg.rows],
  ## pdg.quant$fit.error,
  ## pdg.quant$fit.scalederr,
  ## pdg.quant$scale
  pdg.swb$quant[pdg.swb.sel],
  pdg.swb$quant.err[pdg.swb.sel],
  pdg.swb$quant3.err[pdg.swb.sel],
  pdg.swb$sfact3.types[pdg.swb.sel]
  )

if (FALSE) {
show(rbind(fit       = quant[quant.names.pdg],
           pdg       = pdg.quant[quant.names.pdg],
           fit.err   = quant.err[quant.names.pdg],
           pdg.err   = pdg.quant.err[quant.names.pdg],
           fit.err2  = quant2.err[quant.names.pdg],
           pdg.err2  = pdg.quant2.err[quant.names.pdg],
           fit.sf    = sfact3.types[quant.names.pdg],
           pdg.sf    = pdg.sfact3.types[quant.names.pdg],
           "Dquant%" = (quant[quant.names.pdg]/pdg.quant[quant.names.pdg]-1)*100,
           "Derr%"   = (quant.err[quant.names.pdg]/pdg.quant.err[quant.names.pdg]-1)*100,
           "Derr2%"  = (quant2.err[quant.names.pdg]/pdg.quant2.err[quant.names.pdg]-1)*100,
           "Dsf2%"   = (sfact3.types[quant.names.pdg]/pdg.sfact3.types[quant.names.pdg]-1)*100
           ))

show(rbind(fit       = quant[quant.names.pdg],
           pdg       = pdg.quant[quant.names.pdg],
           fit.err   = quant.err[quant.names.pdg],
           pdg.err   = pdg.quant.err[quant.names.pdg],
           fit.err2  = quant2.err[quant.names.pdg],
           pdg.err2  = pdg.quant2.err[quant.names.pdg],
           fit.sf    = sfact3.types[quant.names.pdg],
           pdg.sf    = pdg.sfact3.types[quant.names.pdg],
           "Dquant%" = sprintf("%7.5f", (quant[quant.names.pdg]/pdg.quant[quant.names.pdg]-1)*100),
           "Derr%"   = sprintf("%7.5f", (quant.err[quant.names.pdg]/pdg.quant.err[quant.names.pdg]-1)*100),
           "Derr2%"  = sprintf("%7.5f", (quant2.err[quant.names.pdg]/pdg.quant2.err[quant.names.pdg]-1)*100),
           "Dsf2%"   = sprintf("%7.5f", (sfact3.types[quant.names.pdg]/pdg.sfact3.types[quant.names.pdg]-1)*100)
           ))
}
