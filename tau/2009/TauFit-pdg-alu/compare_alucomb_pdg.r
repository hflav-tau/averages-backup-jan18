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
file.name.pdg = "pdginput/readpdg.cc"

load(file.name)
lines = readLines(file.name.pdg)
lines.results = suppressWarnings(grep("push_back.*basefitvalue", lines, perl=TRUE))
lines = lines[lines.results]

pdg = gsub(paste("^.*basegamma.push_back.\\s*(\\S+)\\s*.;",
           ".*baseseed.ibase.\\s*=\\s*(\\S+)\\s*;",
           ".*basefitvalue.ibase.\\s*=\\s*(\\S+)\\s*;",
           ".*basefiterror.ibase.\\s*=\\s*(\\S+)\\s*;",
           ".*baserescalederror.ibase.\\s*=\\s*(\\S+)\\s*;",
           ".*basescalefactor.ibase.\\s*=\\s*(\\S+)\\s*;",
           ".*$",
           sep=""),
     "\\1;\\2;\\3;\\4;\\5;\\6", lines, perl=TRUE)

pdg.num = t(sapply(strsplit(pdg, ";"), as.numeric))
pdg.quant.names = sprintf("Gamma%d", pdg.num[,1])

pdg.quant.seed = pdg.num[,2]
pdg.quant = pdg.num[,3]
pdg.quant.err = pdg.num[,4]
pdg.quant2.err = pdg.num[,5]
pdg.sfact.row = pdg.num[,6]

names(pdg.quant.seed) = pdg.quant.names
names(pdg.quant) = pdg.quant.names
names(pdg.quant.err) = pdg.quant.names
names(pdg.quant2.err) = pdg.quant.names
names(pdg.sfact.row) = pdg.quant.names

quant.names = names(quant)
pdg.select = quant.names %in% pdg.quant.names
quant.names.pdg = quant.names[pdg.select]

if (FALSE) {
show(rbind(fit       = quant[quant.names.pdg],
           pdg       = pdg.quant[quant.names.pdg],
           fit.err   = quant.err[quant.names.pdg],
           pdg.err   = pdg.quant.err[quant.names.pdg],
           fit.err2  = quant2.err[quant.names.pdg],
           pdg.err2  = pdg.quant2.err[quant.names.pdg],
           fit.sf    = sfact.row[quant.names.pdg],
           pdg.sf    = pdg.sfact.row[quant.names.pdg],
           "Dquant%" = (quant[quant.names.pdg]/pdg.quant[quant.names.pdg]-1)*100,
           "Derr%"   = (quant.err[quant.names.pdg]/pdg.quant.err[quant.names.pdg]-1)*100,
           "Derr2%"  = (quant2.err[quant.names.pdg]/pdg.quant2.err[quant.names.pdg]-1)*100,
           "Dsf2%"   = (sfact.row[quant.names.pdg]/pdg.sfact.row[quant.names.pdg]-1)*100
           ))

show(rbind(fit       = quant[quant.names.pdg],
           pdg       = pdg.quant[quant.names.pdg],
           fit.err   = quant.err[quant.names.pdg],
           pdg.err   = pdg.quant.err[quant.names.pdg],
           fit.err2  = quant2.err[quant.names.pdg],
           pdg.err2  = pdg.quant2.err[quant.names.pdg],
           fit.sf    = sfact.row[quant.names.pdg],
           pdg.sf    = pdg.sfact.row[quant.names.pdg],
           "Dquant%" = sprintf("%7.5f", (quant[quant.names.pdg]/pdg.quant[quant.names.pdg]-1)*100),
           "Derr%"   = sprintf("%7.5f", (quant.err[quant.names.pdg]/pdg.quant.err[quant.names.pdg]-1)*100),
           "Derr2%"  = sprintf("%7.5f", (quant2.err[quant.names.pdg]/pdg.quant2.err[quant.names.pdg]-1)*100),
           "Dsf2%"   = sprintf("%7.5f", (sfact.row[quant.names.pdg]/pdg.sfact.row[quant.names.pdg]-1)*100)
           ))
}

cat("             quant     err       err2      sf   pdg       pdg.err   pdg.err3  pdg.sf Dquant%    Derr%       Derr2%    Dsf% \n")
rc = mapply(function(
  name,
  quant,
  quant.err,
  quant2.err,
  sfact.row,
  pdg.quant,
  pdg.quant.err,
  pdg.quant2.err,
  pdg.sfact.row
  ) {
  cat(sprintf("%12s %9.7f %9.7f %9.7f %4.2f %9.7f %9.7f %9.7f %4.2f %9.4f%% %9.4f%% %9.4f%% %9.4f%%\n",
              name,
              quant,
              quant.err,
              quant2.err,
              sfact.row,
              pdg.quant,
              pdg.quant.err,
              pdg.quant2.err,
              pdg.sfact.row,
              (quant/pdg.quant-1)*100,
              (quant.err/pdg.quant.err-1)*100,
              (quant2.err/pdg.quant2.err-1)*100,
              (sfact.row/pdg.sfact.row-1)*100
              ))
},
  quant.names.pdg,
  quant[quant.names.pdg],
  quant.err[quant.names.pdg],
  quant2.err[quant.names.pdg],
  sfact.row[quant.names.pdg],
  pdg.quant[quant.names.pdg],
  pdg.quant.err[quant.names.pdg],
  pdg.quant2.err[quant.names.pdg],
  pdg.sfact.row[quant.names.pdg]
  )
