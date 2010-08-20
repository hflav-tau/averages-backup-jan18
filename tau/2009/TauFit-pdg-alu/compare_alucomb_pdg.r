#!/usr/bin/env Rscript

## ////////////////////////////////////////////////////////////////////////////
##
## compare_alucomb_pdg.r
##
## - get averages from alucomb.r .rdata file
## - get PDG averages from pdginput/pdgread.cc and pdginput/s035-fit-no-babar-belle.datx
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
pdg$quant.names = paste("Gamma", pdg$gamma, sep="")

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
  col.names = c("quant.id", "quant.seed", "quant.val", "quant.err", "quant.sf.err", "quant.sf.sfact"))
close(con)
pdg.swb$quant.names = paste("Gamma", pdg.swb$quant.id, sep="")

quant.names = names(quant.val)

pdg.select = quant.names %in% pdg.swb$quant.names
quant.names.pdg = quant.names[pdg.select]
pdg.swb.sel = match(quant.names.pdg, pdg.swb$quant.names)

## pdg.rows = sapply(unique(pdg$quant.names), function(x) which(x == pdg$quant.names)[1])
pdg.rows = match(quant.names.pdg, pdg$quant.names)
pdg.rows = pdg.rows[!is.na(pdg.rows)]

##
## compare results pdg vs. alucomb
##
compare.1 = function() {

cat("           quant     err       err.sf    sf   pdg       pdg.err   pdg.err.sf pdg.sf  Dquant%    Derr%      Derr.sf%   Dsf%\n")
rc = mapply(function(
  name,
  quant.val,
  quant.err,
  quant.sf.err,
  quant.sf.sfact,
  pdg.quant.val,
  pdg.quant.err,
  pdg.quant.sf.err,
  pdg.quant.sf.sfact
  ) {
  cat(sprintf("%10s %9.7f %9.7f %9.7f %4.2f %9.7f %9.7f %9.7f %5.2f %9.4f%% %9.4f%% %9.4f%% %9.4f%%\n",
              name,
              quant.val,
              quant.err,
              quant.sf.err,
              quant.sf.sfact,
              pdg.quant.val,
              pdg.quant.err,
              pdg.quant.sf.err,
              pdg.quant.sf.sfact,
              (quant.val/pdg.quant.val-1)*100,
              (quant.err/pdg.quant.err-1)*100,
              (quant.sf.err/pdg.quant.sf.err-1)*100,
              (quant.sf.sfact/pdg.quant.sf.sfact-1)*100
              ))
},
  quant.names.pdg,
  quant.val[quant.names.pdg],
  quant.err[quant.names.pdg],
  quant.sf.err[quant.names.pdg],
  quant.sf.sfact[quant.names.pdg],
  pdg.swb$quant.val[pdg.swb.sel],
  pdg.swb$quant.err[pdg.swb.sel],
  pdg.swb$quant.sf.err[pdg.swb.sel],
  pdg.swb$quant.sf.sfact[pdg.swb.sel]
  )
}

##
## compare results pdg vs. alucomb
##
compare.2 = function() {

cat("            pdg   or.sc or.fq al.sc al.fq or.fu al.fu\n")
rc = mapply(function(
  name,
  pdg.quant.sf.sfact,
  orin.sc.quant.sfact,
  orin.fq.quant.sfact,
  alu.sc.quant.sfact,
  alu.fq.quant.sfact,
  orin.full.quant.sfact,
  alu.full.quant.sfact
  ) {
  cat(sprintf("%10s %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n",
              name,
              pdg.quant.sf.sfact,
              orin.sc.quant.sfact,
              orin.fq.quant.sfact,
              alu.sc.quant.sfact,
              alu.fq.quant.sfact,
              orin.full.quant.sfact,
              alu.full.quant.sfact
              ## (quant.val/pdg.quant.val-1)*100,
              ## (quant.err/pdg.quant.err-1)*100,
              ## (quant.sf.err/pdg.quant.sf.err-1)*100,
              ## (quant.sf.sfact/pdg.quant.sf.sfact-1)*100
              ))
},
  quant.names.pdg,
  pdg.swb$quant.sf.sfact[pdg.swb.sel],
  orin.sc$quant.sfact[quant.names.pdg],
  orin.fq$quant.sfact[quant.names.pdg],
  alu.sc$quant.sfact[quant.names.pdg],
  alu.fq$quant.sfact[quant.names.pdg],
  orin.full$quant.sfact[quant.names.pdg],
  alu.full$quant.sfact[quant.names.pdg]
  )
}

compare.2()

##
## obsolete code follows
##

if (FALSE) {

show(rbind(fit       = quant.val[quant.names.pdg],
           pdg       = pdg.quant.val[quant.names.pdg],
           fit.err   = quant.err[quant.names.pdg],
           pdg.err   = pdg.quant.err[quant.names.pdg],
           fit.err2  = quant2.err[quant.names.pdg],
           pdg.err2  = pdg.quant2.err[quant.names.pdg],
           fit.sf    = quant.sf.sfact[quant.names.pdg],
           pdg.sf    = pdg.quant.sf.sfact[quant.names.pdg],
           "Dquant%" = (quant[quant.names.pdg]/pdg.quant.val[quant.names.pdg]-1)*100,
           "Derr%"   = (quant.err[quant.names.pdg]/pdg.quant.err[quant.names.pdg]-1)*100,
           "Derr2%"  = (quant2.err[quant.names.pdg]/pdg.quant2.err[quant.names.pdg]-1)*100,
           "Dsf2%"   = (quant.sf.sfact[quant.names.pdg]/pdg.quant.sf.sfact[quant.names.pdg]-1)*100
           ))

show(rbind(fit       = quant[quant.names.pdg],
           pdg       = pdg.quant.val[quant.names.pdg],
           fit.err   = quant.err[quant.names.pdg],
           pdg.err   = pdg.quant.err[quant.names.pdg],
           fit.err2  = quant2.err[quant.names.pdg],
           pdg.err2  = pdg.quant2.err[quant.names.pdg],
           fit.sf    = quant.sf.sfact[quant.names.pdg],
           pdg.sf    = pdg.quant.sf.sfact[quant.names.pdg],
           "Dquant%" = sprintf("%7.5f", (quant.val[quant.names.pdg]/pdg.quant.val[quant.names.pdg]-1)*100),
           "Derr%"   = sprintf("%7.5f", (quant.err[quant.names.pdg]/pdg.quant.err[quant.names.pdg]-1)*100),
           "Derr2%"  = sprintf("%7.5f", (quant2.err[quant.names.pdg]/pdg.quant2.err[quant.names.pdg]-1)*100),
           "Dsf2%"   = sprintf("%7.5f", (quant.sf.sfact[quant.names.pdg]/pdg.quant.sf.sfact[quant.names.pdg]-1)*100)
           ))

##---

print.wo.rownames = function(df) {
  matr = as.matrix(df)
  rownames(matr) = rep("", nrow(matr))
  print(matr, quote=FALSE, right=TRUE)
}

options.save = options()
options(width=132)

rc = print.wo.rownames(data.frame(
  Gamma=as.character(pdg$quant.names),
  Descr=as.character(pdg$descr)))

options(options.save)

##---

pdg.rows = sapply(unique(pdg$quant.names[order(pdg$quant.names)]), function(x) which(x == pdg$quant.names)[1])
quant.name.maxwidth = max(nchar(pdg$quant.names[pdg.rows]))
node.maxwidth = max(nchar(as.character(pdg$node[pdg.rows])))
rc = mapply(function(node, quant.name, descr) {
  cat(sprintf(paste("%-", node.maxwidth+1, "s %-", quant.name.maxwidth+1, "s %s\n", sep=""), node, quant.name, descr))
}, pdg$node[pdg.rows], pdg$quant.names[pdg.rows], pdg$descr[pdg.rows])

}

