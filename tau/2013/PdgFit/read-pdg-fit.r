#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(stringr))
suppressPackageStartupMessages(require(readr))
source("~/repo/hfagtau-averages/Common/bin/alucomb2-utils.r")

##
## functions
##

##--- read saved variables inside a list
load.in.list <- function(.file.name) { load(.file.name); as.list(environment()) }

##--- compute relative difference between two numbers
rel.diff = function(x, y) {
  (y-x)/(abs(x) + abs(y))
}

##
## deparse expression and produce single line
##
deparse.one.line = function(expr) {
  paste(gsub("^\\s+|\\s+$", "", sapply(as.expression(expr), function(x) deparse(x)), perl=TRUE), collapse="")
}

##
## return numeric id for sorting labels like "Gamma5", "Gamma3by5"
## <n>by<m> are sorted after <n> in ascending order ov <m>
##
alurep.gamma.num.id = function(gamma.name) {
  gamma.name = sub("Unitarity", "Gamma1000", gamma.name, fixed=TRUE)
  gamma.name = sub("GammaAll", "Gamma999", gamma.name, fixed=TRUE)
  rc = str_match(gamma.name, "\\D+(\\d+)by(\\d+)|\\D+(\\d+)")
  if (length(rc) == 0) return(numeric(0))
  num = ifelse(!is.na(rc[,4]), as.numeric(rc[,4])*1000, as.numeric(rc[,2])*1000 + as.numeric(rc[,3]))
  return(num)
}

##
## main code
##

##
## get a section of the PDG fit log
##
get.section = function(input.txt, match.beg, match.end, skip.beg=2, skip.end=1) {
  for (match.str in match.beg) {
    l1 = grep(match.str, input.txt)[1]
    if (is.na(l1)) {
      stop("Cannot find sequence matches '", paste(match.neg, collapse="/"), "'")
    }
    input.txt = input.txt[-(1:l1)]
  }
  input.txt = input.txt[-(1:skip.beg)]
  l2 = grep(match.end, input.txt)[1]
  if (is.na(l2)) {
    stop("Cannot find end match '", match.end, "'")
  }
  input.txt[1:(l2-skip.end)]
}

##
## convert string to double or zero if blank
##
parse_double_or_blank = function(str) {
  ifelse(grepl("^\\s*$", str), 0., parse_double(str))
}

##
## main code
##

##--- get PDG fit log
fit.txt = readLines("s035-fit-ns-all-2015-01-28.fit")

##
## get HFAG data
##

##--- HFAG fit reproducing PDG
hfag.data.fname = "pdgfit-pdg-meas.rdata"
ht = load.in.list(hfag.data.fname)
cat(paste0("HFAG reproducing PDG fit file '", hfag.data.fname, "' read\n"))

##--- HFAG reference fit
hfag.data.ref.fname = "../TauFit/average2-aleph-hcorr_13-ref.rdata"
htref = load.in.list(hfag.data.ref.fname)
cat(paste0("HFAG fit reference file '", hfag.data.ref.fname, "' read\n"))

##--- HFAG fit using PDG measurements
hfag.data.ref2.fname = "pdgfit-hfag-meas2.rdata"
htref2 = load.in.list(hfag.data.ref2.fname)
cat(paste0("HFAG fit reference file '", hfag.data.ref2.fname, "' read\n"))

##--- HFAG fit, no hcorr, no Belle incl K0S, with unitarity
hfag.data.ref3.fname = "pdgfit-step-0.rdata"
htref3 = load.in.list(hfag.data.ref3.fname)
cat(paste0("HFAG fit reference file '", hfag.data.ref3.fname, "' read\n"))

if (FALSE) {
meas.val = sapply(ht$measurements, function(x) {x$value})
meas.stat = sapply(ht$measurements, function(x) {x$stat})
meas.stat.p = sapply(ht$measurements, function(x) {x$stat.p})
meas.stat.n = sapply(ht$measurements, function(x) {x$stat.n})
meas.syst = sapply(ht$measurements, function(x) {x$syst})
meas.syst.p = sapply(ht$measurements, function(x) {x$syst.p})
meas.syst.n = sapply(ht$measurements, function(x) {x$syst.n})
meas.err = sqrt(meas.stat^2 + meas.syst^2)
}

meas.gamma = sapply(ht$measurements, function(x) {x$tags[2]})
meas.val = sapply(ht$measurements, function(x) {ifelse(is.null(x$value.orig), x$value, x$value.orig)})
meas.stat = sapply(ht$measurements, function(x) {x$stat})
meas.stat.p = sapply(ht$measurements,
  function(x) {ifelse(attr(x$stat.p, "input") != "", parse_double(attr(x$stat.p, "input")), x$stat)})
meas.stat.n = sapply(ht$measurements,
  function(x) {ifelse(attr(x$stat.n, "input") != "", parse_double(attr(x$stat.n, "input")), -x$stat)})
meas.syst = sapply(ht$measurements, function(x) {ifelse(is.null(x$syst.orig), x$syst, x$syst.orig)})
meas.syst.p = sapply(ht$measurements,
  function(x) {ifelse(attr(x$syst.p, "input") != "",
                      parse_double(attr(x$syst.p, "input")),
                      ifelse(is.null(x$syst.orig), x$syst, x$syst.orig))})
meas.syst.n = sapply(ht$measurements,
  function(x) {ifelse(attr(x$syst.n, "input") != "",
                      parse_double(attr(x$syst.n, "input")),
                      ifelse(is.null(x$syst.orig), -x$syst, -x$syst.orig))})
meas.err = sqrt(meas.stat^2 + meas.syst^2)

##
## get input parameters
##
input.params = get.section(fit.txt, "^ INPUT PARAMETERS", "^\\s*$")

df.params = data.frame(
  count =   parse_number(str_sub(input.params,  1,  7)),
  parcode = str_trim(str_sub(input.params,  8, 16)),
  param = parse_number(str_sub(input.params, 17, 25)),
  seed =  parse_double(str_sub(input.params, 26, 39)),
  descr = str_trim(str_sub(input.params, 40)),
  stringsAsFactors = FALSE
  )

##--- set input parameters parcode to preceding row parcode if missing
df.params$parcode = mapply(
  function(i) {
    parcode = df.params$parcode[i]
    previous.parcodes = df.params$parcode[i:1]
    ifelse(parcode=="",
           previous.parcodes[which("" != previous.parcodes)[1]],
           parcode)
  }, 1:nrow(df.params))

##--- add label summarizing parcode and param
df.params$plab = paste0(df.params$parcode, "_p", df.params$param)

##--- change description to canonical text description form
df.params$descr = gsub("tau- --> (.*)", "G(\\1)/G(total)", df.params$descr)

##
## get input nodes
##
input.nodes = get.section(fit.txt, "^ INPUT NODES", "^\\s*$")

df.nodes = data.frame(
  node.count =   parse_number(str_sub(input.nodes,  1,  7)),
  node = str_trim(str_sub(input.nodes,  8, 15)),
  eq = str_trim(str_sub(input.nodes, 16, 17)),
  adj =  parse_double(str_sub(input.nodes, 18, 21)),
  un = str_trim(str_sub(input.nodes, 22, 27)),
  count = parse_number(str_sub(input.nodes, 28, 32)),
  parcode = str_trim(str_sub(input.nodes, 33, 42)),
  param = parse_number(str_sub(input.nodes, 43, 47)),
  sum = parse_number(str_sub(input.nodes, 51, 53)),
  coeff = parse_double(str_sub(input.nodes, 55, 65)),
  coeff.extra = parse_double(str_sub(input.nodes, 66, 76)),
  descr = str_trim(str_sub(input.nodes, 80)),
  stringsAsFactors = FALSE
  )

##--- add label summarizing parcode and param
df.nodes$plab = paste0(df.nodes$parcode, "_p", df.nodes$param)

##
## get input measurements
##
input.meas = get.section(fit.txt, "^ INPUT MEASUREMENTS", "^\\s*$")

##--- initially, put stat and syst errors in positive errors
df.meas = data.frame(
  count =   parse_number(str_sub(input.meas,  1,  7)),
  node = str_trim(str_sub(input.meas,  8, 15)),
  author = str_trim(str_sub(input.meas, 16, 27)),
  year =  str_trim(str_sub(input.meas, 28, 32)),
  val = parse_double(str_sub(input.meas, 35, 45)),
  stat.sign = str_trim((str_sub(input.meas, 47, 48))),
  statp = parse_double_or_blank(str_sub(input.meas, 49, 56)),
  syst.sign = str_trim((str_sub(input.meas, 58, 59))),
  systp = parse_double_or_blank(str_sub(input.meas, 60, 67)),
  descr = str_trim(str_sub(input.meas, 89)),
  stringsAsFactors = FALSE
  )

##--- negative errors set equal to positive errors
df.meas$statn = df.meas$statp
df.meas$systn = df.meas$systp

##--- identify where are asymmetric errors
stat.asymm.p = which(df.meas$stat.sign == "+")
stat.asymm.n = which(df.meas$stat.sign == "-")
if (any(stat.asymm.n != stat.asymm.p+1)) {
  stop("asymmetric errors do not respect having first + then -")
}

df.meas[stat.asymm.p, "statn"] = df.meas[stat.asymm.n, "statn"]
df.meas = df.meas[-stat.asymm.n, ]

syst.asymm.p = which(df.meas$syst.sign == "+")
syst.asymm.n = which(df.meas$syst.sign == "-")
if (any(syst.asymm.n != syst.asymm.p+1)) {
  stop("asymmetric errors do not respect having first + then -")
}

##--- set negative error on line where the positive error was set
df.meas[syst.asymm.p, "systn"] = df.meas[syst.asymm.n, "systn"]
##--- remove lines where the negative error was set
df.meas = df.meas[-syst.asymm.n, ]

##--- compute squared root averages of positive and negative errors
df.meas$stat = sqrt( (df.meas$statp^2 + df.meas$statn^2)/2 )
df.meas$syst = sqrt( (df.meas$systp^2 + df.meas$systn^2)/2 )
##--- total error
df.meas$err = sqrt( df.meas$stat^2 + df.meas$syst^2 )

##--- set node on measurements that miss it using the last defined one
df.meas$node = unname(
  mapply(
    function(i) {
      node = df.meas$node[i]
      previous.nodes = df.meas$node[i:1]
      ifelse(node=="",
             previous.nodes[which("" != previous.nodes)[1]],
             node)
    }, 1:nrow(df.meas)))

##
## get input measurements correlations
##

input.meas.corr.txt = get.section(fit.txt, "^ INPUT CORRELATIONS BETWEEN MEASUREMENTS", "^ FIRST FIT", skip.beg=1, skip.end=2)
ii.beg = grep("^ INDEX", input.meas.corr.txt) + 1
ii.end = grep("^\\s*$", input.meas.corr.txt) - 1
input.meas.corr.list = mapply(
  function(section.beg, section.end) {
    section = input.meas.corr.txt[section.beg:section.end]
    ## cat(section, collapse="\n", sep="\n")
    ii.cols = parse_integer(str_match_all(section, "M\\s*(\\d+)")[[1]][,2])
    tmp = str_match(section[-(1:2)], "M\\s*(\\d+)[*]\\s+(.*)" )
    ii.rows = parse_integer(tmp[,2])
    corr = lapply(strsplit(tmp[,3], "\\s+"), parse_double)
    rc = lapply(1:length(ii.cols),
      function(ii) {
        rc = lapply(1:ii,
          function(jj) {
            c(row=ii.rows[ii], col=ii.cols[jj], corr=corr[[ii]][jj])
          })
      })
    do.call(c, rc)
  }, ii.beg, ii.end, SIMPLIFY=FALSE)
input.meas.corr.list = do.call(c, input.meas.corr.list)
df.meas$stat.corr = lapply(1:length(df.meas$node),
  function(ii.meas) {    
    which.rows = which(ii.meas == sapply(input.meas.corr.list, function(el) el["row"]))
    meas.corr.rows = lapply(which.rows, function(ii.el) { input.meas.corr.list[[ii.el]][c("col", "corr")]})
    which.cols = which(ii.meas == sapply(input.meas.corr.list, function(el) el["col"]))
    meas.corr.cols = lapply(which.cols, function(ii.el) { input.meas.corr.list[[ii.el]][c("row", "corr")]})
    meas.corr = c(meas.corr.rows, meas.corr.cols)
  })

##
## get fit result summary
##

input.fit.summary = get.section(fit.txt, "^ FINAL RESULTS", "^\\s*$")
val = str_match(input.fit.summary[2], "FINAL CHISQUARE[=\\s]+(\\d+[.]*\\d*) FOR (\\d+)")[2:3]
pdg.chisq = parse_double(val[1])
pdg.dof = parse_integer(val[2])
pdg.chisq.prob = parse_double(str_match(input.fit.summary[3], "FINAL CONFIDENCE\\D+(\\d+[.]*\\d*)")[2])

##
## get fit results for parameters
##
input.params.fit = get.section(fit.txt, c("^ FINAL RESULTS", "^ RESULTS FOR PARAMETERS"), "^\\s*$")

df.params.fit = data.frame(
  count =   parse_number(str_sub(input.params.fit,  1,  4)),
  parcode = str_trim(str_sub(input.params.fit,  6, 9)),
  param = parse_double(str_sub(input.params.fit, 11, 18)),
  val =  parse_double(str_sub(input.params.fit, 19, 31)),
  err =  parse_double(str_sub(input.params.fit, 34, 43)),
  errp =  parse_double(str_sub(input.params.fit, 46, 55)),
  errn =  parse_double(str_sub(input.params.fit, 58, 67)),
  units =  parse_double(str_sub(input.params.fit, 70, 80)),
  scale = parse_double(str_sub(input.params.fit, 81, 86)),
  descr = str_trim(str_sub(input.params.fit, 88)),
  stringsAsFactors = FALSE
  )

##--- change description to canonical text description form
df.params.fit$descr = gsub("tau- --> (.*)", "G(\\1)/G(total)", df.params.fit$descr)

##
## get fit results for nodes
##
input.nodes.fit = get.section(fit.txt, c("^ FINAL RESULTS", "^ RESULTS FOR NODES"), "^\\s*$")

df.nodes.fit = data.frame(
  count =   parse_number(str_sub(input.nodes.fit,  1,  4)),
  node = str_trim(str_sub(input.nodes.fit,  10, 18)),
  val =  parse_double(str_sub(input.nodes.fit, 19, 31)),
  err =  parse_double(str_sub(input.nodes.fit, 33, 43)),
  errp =  parse_double(str_sub(input.nodes.fit, 45, 55)),
  errn =  parse_double(str_sub(input.nodes.fit, 57, 67)),
  units =  parse_double(str_sub(input.nodes.fit, 70, 80)),
  scale = parse_double(str_sub(input.nodes.fit, 81, 86)),
  descr = str_trim(str_sub(input.nodes.fit, 88)),
  stringsAsFactors = FALSE
  )

##
## get HFAG input nodes
##

##--- HFAG nodes
quant.node = sapply(ht$combination$quantities, function(quant) {ifelse(is.null(quant$node), "", quant$node)})
quant.gamma = names(quant.node)
names(quant.gamma) = quant.node[quant.gamma]
quant.node.defined = quant.node[quant.node != ""]

##--- HFAG modes text description
quant.descr = sapply(ht$combination$quantities, function(quant) quant$descr)

##--- remove spaces to match PDG
quant.descr = gsub(" / ", "/", quant.descr, fixed=TRUE)
quant.descr = gsub("ex. ", "ex.", quant.descr, fixed=TRUE)
quant.descr = gsub(", ", ",", quant.descr, fixed=TRUE)

##--- shorter descripion in the form B()
quant.descr.br = gsub("/G(total)", "", quant.descr, fixed=TRUE)
quant.descr.br = gsub("G(", "B(", quant.descr.br, fixed=TRUE)

##--- put eta always first, to match PDG format
quant.descr = gsub("G\\((.+) eta (.*)nu\\(tau\\)", "G(eta \\1 \\2nu(tau)", quant.descr)

##--- build array to convert from descr to gamma in HFAG context
quant.gamma.descr = names(quant.descr)
names(quant.gamma.descr) = quant.descr

##
## compute distance between a PDG and a HFAG measurement
## (quadratic sum of value and error asymmetry)
##
val.cmp = outer(df.meas$val, meas.val, rel.diff)
err.cmp = outer(df.meas$err, meas.err, rel.diff)
tot.cmp = sqrt(val.cmp^2 + err.cmp^2)

##--- for each PDG measurement get closest HFAG measurement
best.matches = sapply(seq(1, nrow(val.cmp)), function(i) {order(tot.cmp[i, ])[1]})

##--- there must not be any duplication
if (any(duplicated(best.matches))) {
  stop("more than one measurement matches a PDG measurement")
}

##
## compute matched measurements, starting with largest discrepancies
##
meas.matched.names = names(meas.val[best.matches])

meas.diff = sapply(1:nrow(tot.cmp), function(i) {tot.cmp[i, best.matches[i]]})
meas.diff.order = order(meas.diff, decreasing=TRUE)
meas.diff.order = meas.diff.order[meas.diff[meas.diff.order] > 1e-4]

## ggplot(data.frame(x=meas.diff[abs(meas.diff)>1e-12]), aes(x = x)) + geom_histogram() + scale_x_log10()

df.meas.cmp = data.frame(
  pdg = paste(df.meas$author, df.meas$year),
  node = df.meas$node,
  hfag = meas.matched.names,
  gamma = meas.gamma[best.matches],

  pdg.val = df.meas$val,
  hfag.val = meas.val[best.matches],

  pdg.err = df.meas$err,
  hfag.err = meas.err[best.matches],

  pdg.stat = df.meas$stat,
  hfag.stat = meas.stat[best.matches],

  pdg.syst = df.meas$syst,
  hfag.syst = meas.syst[best.matches],

  pdg.statp = df.meas$statp,
  hfag.statp = meas.stat.p[best.matches],

  pdg.statn = -df.meas$statn,
  hfag.statn = meas.stat.n[best.matches],

  pdg.systp = df.meas$systp,
  hfag.systp = meas.syst.p[best.matches],

  pdg.systn = -df.meas$systn,
  hfag.systn = meas.syst.n[best.matches],

  stringsAsFactors=FALSE
  )
rownames(df.meas.cmp) = NULL

##
## PDG input nodes, check consistency with HFAG
##

##--- get list of defined input nodes
nodes.def = (df.nodes$descr != "")
## nodes.def = df.nodes$eq == "+" | df.nodes$eq == "/"

quant.pdg.descr = df.nodes$descr[nodes.def]
quant.pdg.node = df.nodes$node[nodes.def]
quant.pdg.plab = df.nodes$plab[nodes.def]
##--- can convert from description (useful for fit parameters)
names(quant.pdg.node) = quant.pdg.descr
names(quant.pdg.descr) = quant.pdg.node
names(quant.pdg.plab) = quant.pdg.node

##--- node definitions to be added to HFAG
nodes.to.define = which(! (quant.pdg.node %in% quant.node.defined))
if (length(nodes.to.define) > 0) {
  cat("please define:\n")
  print(cbind(quant.pdg.node[nodes.to.define]))
  stop()
}

##--- check that HFAG and PDG descriptions match
quant.hfag.descr = quant.descr[quant.gamma[quant.pdg.node]]
non.matching = gsub("\\s+", "", quant.hfag.descr) != gsub("\\s+", "", quant.pdg.descr)

##
## elaboration on PDG parameters
##

##--- recover node from description
df.params$node = quant.pdg.node[df.params$descr]
df.params.fit$node = quant.pdg.node[df.params.fit$descr]

##--- utility vectors to translate between plab and node for parameters
params.node = df.params$node
names(params.node) = df.params$plab
params.plab = df.params$plab
names(params.plab) = df.params$node

##--- PDG nodes in HFAG gamma notation
quant.gamma.nodes = quant.gamma[quant.pdg.node]
##--- PDG fit parameters in HFAG gamma notation
quant.gamma.params = quant.gamma[params.node]
##--- PDG nodes that are not fit parameters in HFAG gamma notation
quant.gamma.nodes.nonpar = setdiff(quant.gamma.nodes, quant.gamma.params)

##
## print PDG fit parameters definitions
##

cat("##\n")
cat("## PDG fit parameters definitions\n")
cat("##\n")

cat(paste(
  format(names(params.plab), width=10),
  format(params.plab, width=10),
  quant.pdg.descr[names(params.plab)], collapse="\n"
  ))
cat("\n")

##
## print PDG nodes definitions (nodes that are not fit parameters)
##

cat("##\n")
cat("## PDG definitions of nodes that are not fit parameters\n")
cat("##\n")

cat(paste(
  format(quant.node[quant.gamma.nodes.nonpar], width=10),
  format(quant.pdg.plab[quant.node[quant.gamma.nodes.nonpar]], width=10),
  quant.pdg.descr[quant.node[quant.gamma.nodes.nonpar]], collapse="\n"
  ))
cat("\n")

##
## print PDG vs. HFAG input measurements mismatches
##

cat("##\n")
cat("## PDG vs. HFAG input measurements mismatches >1e-4, in decreasing size\n")
cat("## (mismatches are quadrature sum of the relative discrepancies in value and total uncertainty)\n")
cat("##\n")
## print(df.meas.cmp[meas.diff.order,])

if (length(meas.diff.order) > 0) {
rc = lapply(split(df.meas.cmp[meas.diff.order,], 1:length(meas.diff.order)),
  function(cmp) {
    cat(cmp$node, cmp$pdg, cmp$gamma)
    cat("\n")
    print(cmp)
    if (cmp$pdg.statp != cmp$pdg.stat || -cmp$pdg.statn != cmp$pdg.stat) {
      pdg.stat = paste(sprintf("%+.8g", c(cmp$pdg.statp, cmp$pdg.statn)))
    } else {
      pdg.stat = paste(sprintf("+-%.8g", c(cmp$pdg.stat)))
    }

    if (cmp$pdg.systp != cmp$pdg.syst || -cmp$pdg.systn != cmp$pdg.syst) {
      pdg.syst = paste(sprintf("%+.8g", c(cmp$pdg.systp, cmp$pdg.systn)))
    } else {
      pdg.syst = paste(sprintf("+-%.8g", c(cmp$pdg.syst)))
    }
    
    if (cmp$hfag.statp != cmp$hfag.stat || -cmp$hfag.statn != cmp$hfag.stat) {
      hfag.stat = paste(sprintf("%+.8g", c(cmp$hfag.statp, cmp$hfag.statn)))
    } else {
      hfag.stat = paste(sprintf("+-%.8g", c(cmp$hfag.stat)))
    }

    if (cmp$hfag.systp != cmp$hfag.syst || -cmp$hfag.systn != cmp$hfag.syst) {
      hfag.syst = paste(sprintf("%+.8g", c(cmp$hfag.systp, cmp$hfag.systn)))
    } else {
      hfag.syst = paste(sprintf("+-%.8g", c(cmp$hfag.syst)))
    }
    
    cat("  pdg  ", format(c(paste(sprintf("%.8g", c(cmp$pdg.val, cmp$pdg.err))), pdg.stat, pdg.syst), width=13)) 
    cat("\n")
    cat("  hfag ", format(c(paste(sprintf("%.8g", c(cmp$hfag.val, cmp$hfag.err))), hfag.stat, hfag.syst), width=13)) 
    cat("\n")
  })
}

cat("##\n")
cat("## PDG-HFAG mismatches in node descriptions\n")
cat("##\n")
rc = mapply(
  function(descr.pdg, descr.hfag, gamma.hfag) {
    cat(gamma.hfag, "\n  pdg  ", descr.pdg, "\n  hfag ", descr.hfag, "\n", sep="")
  },
  quant.pdg.descr[non.matching],
  quant.hfag.descr[non.matching],
  names(quant.hfag.descr[non.matching]))

##
## build list of MODMEAS KEEP cards to define measurements to be used
##

pdg.meas.tags = sapply(ht$measurements[meas.matched.names], function(meas) {paste(meas$tags, collapse=" ")})

if (FALSE) {
  fname = "pdg-meas-cards.input"
  keep.cards.fh = file(fname, "w")
  cat(paste("MODMEAS KEEP", pdg.meas.tags, collapse="\n"), file=keep.cards.fh)
  cat("\n", file=keep.cards.fh)
  close(keep.cards.fh)
  cat("file", fname, "produced\n")
} else {
  cat("##\n")
  cat("## HFAG measurements matching PDG input measurements\n")
  cat("##\n")
  cat(paste("MODMEAS KEEP", pdg.meas.tags, collapse="\n"))
  cat("\n")
}

##
## list HFAG fit quantities matching PDG input nodes
##

quant.gamma.combine = c(quant.gamma.params, quant.gamma.nodes.nonpar)

quant.gamma.combine.maxlen = max(nchar(quant.gamma.combine))
gamma.per.line = trunc((80-2) / (quant.gamma.combine.maxlen+2))
range.beg = seq(1, length(quant.gamma.combine), by=gamma.per.line)
range.end = ifelse(range.beg+gamma.per.line-1 > length(quant.gamma.combine),
  length(quant.gamma.combine), range.beg+gamma.per.line-1)
cat("##\n")
cat("## HFAG quantities corresponding to PDG fit nodes\n")
cat("##\n")
cat("COMBINE\n")
rc = mapply(
  function(beg, end) {
    cat("  ", paste(format(quant.gamma.combine[beg:end], width=quant.gamma.combine.maxlen), collapse=" "), "\n", sep="")
  }, range.beg, range.end)
## cat(paste(quant.gamma[quant.pdg.node], collapse=" "))
cat("\n")

##
## PDG input nodes, find parameter expressions matching nodes definitions coefficients
##

##--- HFAG parameters for expressions substitutions
ht.params = lapply(ht$combination$params, function(x) unname(x["value"]))

##
## parameter expressions to replace numerical PDG coefficients in nodes
##
params.expr.txt = c(
  "1", "2", "(1/2)",
  "(BR_KS_2piz*BR_KS_2piz)",
  "(BR_KS_2piz/2)", "(BR_KS_pimpip/2)",
  "(2*BR_KS_pimpip*BR_KS_2piz)",
  "(1+BR_KS_2piz*BR_KS_2piz)",
  "((1+BR_KS_2piz)/2)",
  "(BR_om_pimpippiz+BR_om_pimpip)",
  "(BR_eta_pimpippiz+BR_eta_3piz)",
  "BR_eta_2gam", "BR_eta_neutral",
  "BR_eta_3piz", "BR_eta_pimpippiz", "BR_eta_charged", "BR_KS_2piz",
  "BR_KS_pimpip", "BR_om_pimpippiz", "BR_om_pimpip", "BR_om_pizgamma",
  "BR_phi_KmKp", "BR_phi_KSKL", "BR_f1_2pizpippim", "BR_f1_2pip2pim"
  )

##--- evaluate parameter expressions
params.expr.val = sapply(params.expr.txt,
  function(expr.txt) {
    eval(esub.expr(parse(text=expr.txt), ht.params))
  })

##--- find expression best approximating PDG coefficients in nodes
rc = sapply(df.nodes$coeff,
  function(coeff) {
    ##+++ patch to get 0.09 identified as omega -> pi0 gamma
    if (coeff == 0.09) {coeff = eval(substitute(BR_om_pizgamma, ht.params))}
    ##+++ patch to get 0.3431 identified as 1/2 * KS -> pi+pi-
    if (coeff == 0.3431) {coeff = eval(substitute(BR_KS_pimpip/2, ht.params))}
    mism = abs(coeff - params.expr.val)*2 / (abs(coeff) + abs(params.expr.val))
    mism = mism[order(mism)[1]]
  })

##--- store expressions for PDG coefficients in nodes
df.nodes$coeff.expr.txt = names(rc)

##
## log how coefficients are replaced
##
cat("##\n")
cat("## Parameter expressions for PDG nodes definitions coefficients\n")
cat("## - 0.09 forced to be matched with B(omega -> pi0 gamma)\n")
cat("## - 0.3431 forced to be matched with 1/2*B(KS -> pi+pi-)\n")
cat("##\n")
coeff.unique.indices = which(!duplicated(df.nodes$coeff))
rc.unique = rc[coeff.unique.indices]
rc.unique.order = order(rc.unique, decreasing=TRUE)
coeff.unique.indices = coeff.unique.indices[rc.unique.order]
cat("      coeff   closest expression                  expression value mismatch\n\n")
rc2 = lapply(coeff.unique.indices,
  function(coeff.ii) {
    cat(paste0(sprintf("%11.6g", df.nodes$coeff[coeff.ii]), " = ",
               format(names(rc[coeff.ii]), width=40),
               " ", sprintf("%11.4g", params.expr.val[names(rc[coeff.ii])]),
               " ", sprintf("%#7.4f%%", 100*rc[coeff.ii]), "\n"))
  })

##
## PDG input nodes, get constraints
##

##--- get begin of node definitions
node.def.beg = which(!is.na(df.nodes$node.count))
##--- get end of node definitions
node.def.end = c(node.def.beg[-1], nrow(df.nodes)+1) - 1

if (FALSE) {
  ##
  ## get PDG constraints with numeric coefficients
  ##
  constr.pdg = mapply(
    function(def.beg, def.end) {
      df.nodes.def = df.nodes[def.beg:def.end, ]
      sel.above = which(df.nodes.def$sum == 1)
      sel.below = which(df.nodes.def$sum == 2)
      ##--- what node is defined
      node.defined = df.nodes[def.beg, ]$node
      ##--- what parameters and coeff. above fraction line
      coeff.above = df.nodes.def[sel.above, ]$coeff
      names(coeff.above) = df.nodes.def[sel.above, ]$plab
      ##--- what parameters and coeff. below fraction line
      coeff.below = df.nodes.def[sel.below, ]$coeff
      names(coeff.below) = df.nodes.def[sel.below, ]$plab
      ##--- return node defined, coefficients and parameters in numerator and denominator
      list(node=node.defined, above=coeff.above, below=coeff.below)
    }, node.def.beg, node.def.end, SIMPLIFY=FALSE)
} else {
  ##
  ## get PDG constraints with symbolic coefficients
  ##
  constr.pdg = mapply(
    function(def.beg, def.end) {
      df.nodes.def = df.nodes[def.beg:def.end, ]
      sel.above = which(df.nodes.def$sum == 1)
      sel.below = which(df.nodes.def$sum == 2)
      ##--- what node is defined
      node.defined = df.nodes[def.beg, ]$node
      ##--- what parameters and coeff. above fraction line
      coeff.above = df.nodes.def[sel.above, ]$coeff.expr.txt
      names(coeff.above) = df.nodes.def[sel.above, ]$plab
      ##--- what parameters and coeff. below fraction line
      coeff.below = df.nodes.def[sel.below, ]$coeff.expr.txt
      names(coeff.below) = df.nodes.def[sel.below, ]$plab
      ##--- return node defined, coefficients and parameters in numerator and denominator
      list(node=node.defined, above=coeff.above, below=coeff.below)
    }, node.def.beg, node.def.end, SIMPLIFY=FALSE)
}

constr.hfag = lapply(constr.pdg,
  function(constr) {
    constr$node = unname(quant.gamma[constr$node])
    names(constr$above) = quant.gamma[params.node[names(constr$above)]]
    names(constr$below) = quant.gamma[params.node[names(constr$below)]]
    constr
  })

constr.human = lapply(constr.pdg,
  function(constr) {
    constr$node = unname(quant.descr.br[quant.gamma[constr$node]])
    names(constr$above) = quant.descr.br[quant.gamma[params.node[names(constr$above)]]]
    names(constr$below) = quant.descr.br[quant.gamma[params.node[names(constr$below)]]]
    constr
  })

constr.to.text = function(constr.list) {
  rc = sapply(constr.list,
    function(constr) {
      if (length(constr$above) == 0) {
        above.txt = "1"
      } else {
        above.txt = paste0(constr$above, "*", names(constr$above))
        ##--- simplify multiplication by 1 and -1
        above.txt = paste0(gsub("^(-|)1[*]", "\\1", above.txt), collapse=" + ")
      }
      if (length(constr$below) == 0) {
        full.txt = above.txt
      } else {
        if (length(constr$above) > 1) {
          above.txt = paste0("(", above.txt, ")")
        }
        below.txt = paste0(constr$below, "*", names(constr$below))
        ##--- simplify multiplication by 1 and -1
        below.txt = paste0(gsub("^(-|)1[*]", "\\1", below.txt), collapse=" + ")
        if (length(constr$below) > 1) {
          below.txt = paste0("(", below.txt, ")")
        }
        full.txt = paste0(above.txt, " / ", below.txt)
      }
      paste0(constr$node, " = ", full.txt)
    })}

constr.pdg.txt = constr.to.text(constr.pdg)
constr.hfag.txt = constr.to.text(constr.hfag)
constr.human.txt = constr.to.text(constr.human)

constr.hfag.trivial = (
  gsub("(\\S+)\\s*=\\s*(.*)$", "\\1", constr.hfag.txt)
  ==
  gsub("(\\S+)\\s*=\\s*(.*)$", "\\2", constr.hfag.txt)
  )
constr.hfag.txt.nontrivial = constr.hfag.txt[!constr.hfag.trivial]
constr.hfag.txt.nontrivial =
  gsub("(\\S+)\\s*=\\s*(.*)$",
       "NLCONSTRAINT \\1.c 0 \"-\\1 + (\\2)\"",
       constr.hfag.txt.nontrivial)

cat("##\n")
cat("## PDG constraints with PDG nodes and parameters\n")
cat("##\n")
cat(paste(constr.pdg.txt, collapse="\n"))
cat("\n")

cat("##\n")
cat("## PDG constraints in human form\n")
cat("##\n")
cat(paste(constr.human.txt, collapse="\n\n"))
cat("\n")

cat("##\n")
cat("## PDG constraints in HFAG gamma notation\n")
cat("##\n")
cat(paste(constr.hfag.txt, collapse="\n"))
cat("\n")

cat("##\n")
cat("## non-identity PDG constraints in HFAG cards notation\n")
cat("##\n")
cat(paste(constr.hfag.txt.nontrivial, collapse="\n"))
cat("\n")
cat("NLCONSTRAINT Unitarity.c 1 \"", paste(quant.gamma[params.node], collapse=" + "), "\"\n", sep="")

##
## compare PDG constraints with HFAG constraints
##

##--- recover NLCONTRAINT directive
htref.constr.str.expr = paste0(
  "NLCONSTRAINT ",
  names(htref$combination$constr.nl.str.expr), " ",
  as.character(htref$combination$constr.nl.str.val), " \"",
  htref$combination$constr.nl.str.expr, "\""
  )
htref.constr.str.expr = gsub("Unitarity ", "Unitarity.c ", htref.constr.str.expr)
##--- otbain equation in simplest form
htref.constr.str.expr = gsub("^NLCONSTRAINT .*[.]c (\\d+) \"(.*)\"", "\\1 = \\2", htref.constr.str.expr)
htref.constr.str.expr = gsub("0 = -(\\S+) [+] \\((.*)\\)", "\\1 = \\2", htref.constr.str.expr)
htref.constr.str.expr = gsub("0 = -(\\S+) [+] (.*)", "\\1 = \\2", htref.constr.str.expr)
##--- get list of nodes in right part of equation
htref.constr.gammas = str_extract_all(htref.constr.str.expr, "(Gamma\\d+by\\d+|Gamma\\d+)")
##--- form string with sorted right side nodes
rc = sapply(htref.constr.gammas,
  function(gammas) {
    paste(sort(gammas[-1]), collapse=" ")
  })
names(rc) = sapply(htref.constr.gammas, function(constr) {constr[1]})
htref.constr.gammas = sort(rc)
names(htref.constr.str.expr) = names(rc)

##--- otbain equation in simplest form
pdg.constr.str.expr = gsub("^NLCONSTRAINT .*[.]c (\\d+) \"(.*)\"", "\\1 = \\2", constr.hfag.txt.nontrivial)
pdg.constr.str.expr = gsub("0 = -(\\S+) [+] \\((.*)\\)", "\\1 = \\2", pdg.constr.str.expr)
##--- get list of nodes in right part of equation
pdg.constr.gammas = str_extract_all(pdg.constr.str.expr, "(Gamma\\d+by\\d+|Gamma\\d+)")
##--- form string with sorted right side nodes
rc = sapply(pdg.constr.gammas,
  function(gammas) {
    paste(sort(gammas[-1]), collapse=" ")
  })
names(rc) = sapply(pdg.constr.gammas, function(constr) {constr[1]})
pdg.constr.gammas = sort(rc)
names(pdg.constr.str.expr) = names(rc)

##
## take two lists of constraint equations in text format
## for each item of the two lists
## - get all the HFAG Gamma symbols
##   first left side, then the right side, sorted by Gamma number
## - return the array of Gamma symbols, left first
## if one list is empty, use just the first one
##
get.gammas.from.two.lists = function(expr.list.1, expr.list.2=NULL) {
  if (!is.null(expr.list.2)) {
    rc = mapply(
      function(elem.1, elem.2) {
        gammas.1 = str_extract_all(elem.1, "(Gamma\\d+by\\d+|Gamma\\d+)")[[1]]
        gammas.2 = str_extract_all(elem.2, "(Gamma\\d+by\\d+|Gamma\\d+)")[[1]]
        right = unique(c(gammas.1[-1], gammas.2[-1]))
        right = right[order(alurep.gamma.num.id(right))]
        left = unique(c(gammas.1[1], gammas.2[1]))
        left = left[order(alurep.gamma.num.id(left))]
        gammas = c(left, right)
      }, expr.list.1, expr.list.2, SIMPLIFY=FALSE)
  } else {
    rc = mapply(
      function(elem.1) {
        gammas = str_extract_all(elem.1, "(Gamma\\d+by\\d+|Gamma\\d+)")[[1]]
        gammas.right = unique(gammas[-1])
        gammas.right = gammas.right[order(alurep.gamma.num.id(gammas.right))]
        gammas = c(gammas[1], gammas.right)
      }, expr.list.1, SIMPLIFY=FALSE)
  }
  rc
}

##
## take a list where each element is a list of Gamma symbols
## for each element of the list
## - assemble text equating the Gamma symbols to their rescription
##
make.gamma.descr.comments = function(expr.list.1, expr.list.2=NULL) {
  gammas = get.gammas.from.two.lists(expr.list.1, expr.list.2)
  rc = sapply(gammas,
    function(gamma.list) {
      rc = paste(gamma.list, "=", quant.descr.br[gamma.list])
      rc = paste(rc, collapse="\n# ")
      paste0("#\n# ", rc, "\n#\n")
    })
  rc
}

##
## now we can compare for the same left side the sorted list of nodes in the right side
##

##--- constraints that have same left side
constr.left.common = intersect(names(htref.constr.gammas), names(pdg.constr.gammas))
constr.non.equal = names(which(htref.constr.gammas[constr.left.common] != pdg.constr.gammas[constr.left.common]))
constr.equal = names(which(htref.constr.gammas[constr.left.common] == pdg.constr.gammas[constr.left.common]))

##--- in PDG but not HFAG
constr.pdg.not.hfag = setdiff(names(pdg.constr.gammas), names(htref.constr.gammas))
##--- in HFAG but not PDG
constr.hfag.not.pdg = setdiff(names(htref.constr.gammas), names(pdg.constr.gammas))

cat("##\n")
cat("## Constraints that are the same on PDG and HFAG\n")
cat("## PDG first, HFAG second\n")
cat("##\n")

cat(
  paste0(
    paste0(
      make.gamma.descr.comments(htref.constr.str.expr[constr.equal]),
      pdg.constr.str.expr[constr.equal], "\n",
      htref.constr.str.expr[constr.equal],
      sep="\n"),
    collapse="\n"))
cat("\n")

##--- flag array indicating what constraints HFAG uses
hfag.constr.flag = htref$combination$constr.all.nl | htref$combination$constr.all.lin
hfag.constr.used = names(htref.constr.str.expr[hfag.constr.flag])

##--- non-equal AND used in HFAG
constr.non.equal.used = intersect(constr.non.equal, hfag.constr.used)
##--- non-equal AND NOT used in HFAG
constr.non.equal.not.used = setdiff(constr.non.equal, hfag.constr.used)

cat("##\n")
cat("## Constraints that are NOT the same on PDG and HFAG, used in HFAG\n")
cat("## PDG first, HFAG second\n")
cat("##\n")

cat(
  paste0(
    paste0(
      make.gamma.descr.comments(pdg.constr.str.expr[constr.non.equal.used], htref.constr.str.expr[constr.non.equal.used]),
      pdg.constr.str.expr[constr.non.equal.used], "\n",
      htref.constr.str.expr[constr.non.equal.used],
      sep="\n"),
    collapse="\n"))
cat("\n")

cat("##\n")
cat("## Constraints that are NOT the same on PDG and HFAG, NOT used in HFAG\n")
cat("## PDG first, HFAG second\n")
cat("##\n")

cat(
  paste0(
    paste0(
      make.gamma.descr.comments(pdg.constr.str.expr[constr.non.equal.not.used], htref.constr.str.expr[constr.non.equal.not.used]),
      pdg.constr.str.expr[constr.non.equal.not.used], "\n",
      htref.constr.str.expr[constr.non.equal.not.used],
      sep="\n"),
    collapse="\n"))
cat("\n")

cat("##\n")
cat("## Constraints in PDG but not HFAG\n")
cat("##\n")

cat(
  paste0(
    paste0(
      make.gamma.descr.comments(pdg.constr.str.expr[constr.pdg.not.hfag]),
      pdg.constr.str.expr[constr.pdg.not.hfag],
      sep="\n"),
    collapse="\n"))
cat("\n")

cat("##\n")
cat("## Constraints in HFAG but not PDG\n")
cat("##\n")

cat(
  paste0(
    paste0(
      make.gamma.descr.comments(htref.constr.str.expr[constr.hfag.not.pdg]),
      htref.constr.str.expr[constr.hfag.not.pdg],
      sep="\n"),
    collapse="\n"))
cat("\n")

##--- replace Gamma HFAG notation with description (not used)
hfag.gamma.to.descr = function(str.expr) {
  sapply(gsub("['`]", "", str.expr),
    function(expr.txt) {
      gsub("`", "",
           deparse.one.line(esub.expr(parse(text=expr.txt), lapply(gsub("['`]", "", quant.descr.br), as.symbol))),
           fixed=TRUE)
    })
}
## rc = hfag.gamma.to.descr(htref.constr.str.expr)

##
## print PDG measurements in HFAG format
##

##--- collect PDG measurements in HFAG format
pdg.measurements = lapply(1:length(df.meas$node),
  function(ii.meas) {
    corr.terms.tot = sapply(unname(df.meas$stat.corr[[ii.meas]]),
      function(x) setNames(x[2]/100, paste(ht$measurements[[best.matches[x[1]]]]$tags, collapse=".")))
    if (length(corr.terms.tot) > 0) {
      corr.terms.tot = corr.terms.tot[order(alurep.gamma.num.id(
        str_match(names(corr.terms.tot), "[^.]+[.](Gamma(\\d+)(by|)\\d*)[.]")[,2]))]
    }
    meas = list(
      value = df.meas$val[ii.meas],
      
      stat = df.meas$stat[ii.meas],
      stat.p = df.meas$statp[ii.meas],
      stat.n = -df.meas$statn[ii.meas],
      
      syst = df.meas$syst[ii.meas],
      syst.p = df.meas$systp[ii.meas],
      syst.n = -df.meas$systn[ii.meas],
      
      quant = ht$measurements[[best.matches[ii.meas]]]$quant,
      tags = ht$measurements[[best.matches[ii.meas]]]$tags,
      params = list(),
      corr.terms.tot = corr.terms.tot[corr.terms.tot != 0]
      )

    attr(meas$value, "input") = sprintf("%.8g", meas$value)

    if (meas$stat.p != -meas$stat.n) {
      attr(meas$stat, "input") = ""
      attr(meas$stat.p, "input") = sprintf("%+.8g", meas$stat.p)
      attr(meas$stat.n, "input") = sprintf("%+.8g", meas$stat.n)
    } else {
      attr(meas$stat, "input") = sprintf("%.8g", meas$stat)
      attr(meas$stat.p, "input") = ""
      attr(meas$stat.n, "input") = ""
    }

    if (meas$syst.p != -meas$syst.n) {
      attr(meas$syst, "input") = ""
      attr(meas$syst.p, "input") = sprintf("%+.8g", meas$syst.p)
      attr(meas$syst.n, "input") = sprintf("%+.8g", meas$syst.n)
    } else {
      attr(meas$syst, "input") = sprintf("%.8g", meas$syst)
      attr(meas$syst.p, "input") = ""
      attr(meas$syst.n, "input") = ""
    }

    attr(meas$corr.terms.tot, "input") = setNames(sprintf("%+.8g", meas$corr.terms.tot), names(meas$corr.terms.tot))

    meas
  })

dump.measurements = function(measurements, file.dir=".") {
  file.naming.list = list(
    babar="babar",
    belle="belle",
    cleo=c("cleo", "cleo3"),
    aleph="aleph",
    opal="opal",
    l3="l3",
    delphi="delphi"
    )
  
  meas.exp = tolower(sapply(measurements, function(meas) {meas$tags[1]}))
  meas.exp.unique = unique(meas.exp)
  file.naming.list$other = setdiff(meas.exp.unique, do.call(c, file.naming.list))

  if (!dir.exists(file.dir)) {dir.create(file.dir)}
  
  for(exp.tag in names(file.naming.list)) {
    exp.fname = paste0(file.path(file.dir, "measurements-"), exp.tag, ".input")
    meas.fh = file(exp.fname, "w")
    sink(file=meas.fh)
    meas.fname = measurements[meas.exp %in% file.naming.list[[exp.tag]]]
    meas.fname.gamma = sapply(meas.fname, function(meas) {meas$tags[2]})
    meas.fname.pub = tolower(sapply(meas.fname, function(meas) {paste0(meas$tags[4:5], collapse="")}))
    meas.fname = meas.fname[order(alurep.gamma.num.id(meas.fname.gamma), meas.fname.pub)]
    for(meas in meas.fname) {
      alucomb2.print.meas(meas, ht$combination$quantities)
      cat("\n")
    }
    sink()
    close(meas.fh)
    cat(paste0("file ", exp.fname, " created\n"))
  }
}

##--- sort stat correlations on gamma number to facilitate comparisons
ht$measurements = lapply(ht$measurements,
  function(meas) {
    if (length(meas$corr.terms.stat) > 0) {
      terms.order = order(alurep.gamma.num.id(
        str_match(names(meas$corr.terms.stat), "[^.]+[.](Gamma(\\d+)(by|)\\d*)[.]")[,2]))
      attr.tmp = attr(meas$corr.terms.stat, "input")[terms.order]
      meas$corr.terms.stat = meas$corr.terms.stat[terms.order]
      attr(meas$corr.terms.stat, "input") = attr.tmp
    }
    if (length(meas$corr.terms.tot) > 0) {
      terms.order = order(alurep.gamma.num.id(
        str_match(names(meas$corr.terms.tot), "[^.]+[.](Gamma(\\d+)(by|)\\d*)[.]")[,2]))
      attr.tmp = attr(meas$corr.terms.tot, "input")[terms.order]
      meas$corr.terms.tot = meas$corr.terms.tot[terms.order]
      attr(meas$corr.terms.tot, "input") = attr.tmp
    }
    meas
  })

dump.measurements(pdg.measurements, "pdg")
dump.measurements(ht$measurements, "hfag")

##
## compare fit results
##

df.nodes.fit = within(df.nodes.fit,
  {
    hfag.quant.val = ht$quant.val[quant.gamma[node]]
    hfag.quant.err = ht$quant.err[quant.gamma[node]]
    mism.val = (val - hfag.quant.val) / sqrt((err^2+hfag.quant.err^2)/2)
    mism.err = (err - hfag.quant.err) / sqrt((err^2+hfag.quant.err^2)/2)
    mism = sqrt(mism.val^2 + mism.err^2)
  })

cat("##\n")
cat("## PDG fit summary\n")
cat("##\n")

cat("chisq / dof =", pdg.chisq, "/", pdg.dof, " ( Prob = ", pdg.chisq.prob, ")\n")
cat("measurements:", nrow(df.meas),
    "fit parameters:", nrow(df.params),
    "fit nodes:", nrow(df.nodes.fit),
    "\n")

cat("##\n")
cat("## HFAG fit summary\n")
cat("##\n")

cat("chisq / dof =", ht$chisq, "/", ht$dof, " ( Prob = ", ht$chisq.prob, ")\n")
cat("measurements:", length(ht$meas.val),
    "fit quantities:", length(ht$quant.val),
    "constraints:", ht$constr.num,
    "\n")

cat("##\n")
cat("## Fit results PDG-HFAG mismatches\n")
cat("## - first line:\n")
cat("##   - (PDG_value - HFAG_value) / <error on value>\n")
cat("##   - (PDG_error - HFAG_error) / <error on value>\n")
cat("##   - PDG fit results asymmetric errors if existing\n")
cat("##\n")
rc = lapply(split(df.nodes.fit, 1:nrow(df.nodes.fit))[order(df.nodes.fit$mism, decreasing=TRUE)],
  function(cmp.fit) {
    line.1 = paste(
      paste(sprintf("%12s", c(cmp.fit$node, quant.gamma[cmp.fit$node])), collapse=" "),
      paste(sprintf("%13.6e", c(cmp.fit$mism.val, cmp.fit$mism.err)), collapse=" "),
      sep=" ")
    line.2 = paste(
      sprintf("%25s", "PDG"),
      paste(sprintf("%13.6e", c(cmp.fit$val, cmp.fit$err)), collapse=" "),
      paste(sprintf("%+13.6e", c(cmp.fit$errp, -cmp.fit$errn)), collapse=" "),
      sep=" ")
    line.3 =  paste(
      sprintf("%25s", "HFAG"),
      paste(sprintf("%13.6e", c(cmp.fit$hfag.quant.val, cmp.fit$hfag.quant.err)), collapse=" "),
      sep=" ")
    cat(paste(line.1, line.2, line.3, "", sep="\n"))
  })

##--- check unitarity
## sum(df.nodes.fit$val[df.params.fit.in.df.nodes.fit])
## sum(df.nodes.fit$hfag.quant.val[df.params.fit.in.df.nodes.fit])

##
## compare constraints with reference HFAG fit
##

df.params.fit.in.df.nodes.fit = sapply(df.params.fit$node, 
  function(node) {
    which(node == df.nodes.fit$node)
  })

##
## find out which measurements have systematic errors
## HFAG 2014: 45 measurements of 191 have listed systematic terms
##

htref2.meas.with.syst = sapply(htref2$measurements,
  function(meas) {
    !is.null(meas$syst.terms)
  })

pdg.meas.with.syst = htref2.meas.with.syst[meas.matched.names]

cat("##\n")
cat("## PDG measurements that have systematics in HFAG \n")
cat("##\n")

cat("  ", paste(names(pdg.meas.with.syst[pdg.meas.with.syst]), collapse="\n  "), "\n", sep="")

##
## list measurements used in HFAG fit but not PDG
## HFAG fit without Belle inclusive K0, no hcorr
##

hfag.meas.not.in.pdg = setdiff(names(htref3$meas.val), meas.matched.names)
meas.name.width.max = max(nchar(hfag.meas.not.in.pdg))

cat("##\n")
cat("## measurements used in HFAG 2014 but not in PDG fit\n")
cat("##\n")

rc = lapply(hfag.meas.not.in.pdg,
  function(meas.name) {
    cat(format(meas.name, width=meas.name.width.max),
        "  ",
        quant.descr.br[htref3$measurements[[meas.name]]$quant],
        "\n", sep="")
  })

##
## list measurements used in PDG fit but not HFAG
## HFAG fit without Belle inclusive K0, no hcorr
##

pdg.meas.not.in.hfag = setdiff(meas.matched.names, names(htref3$meas.val))
meas.name.width.max = max(nchar(pdg.meas.not.in.hfag))

cat("##\n")
cat("## measurements used in PDG fit but not in HFAG 2014\n")
cat("##\n")

rc = lapply(pdg.meas.not.in.hfag,
  function(meas.name) {
    cat(format(meas.name, width=meas.name.width.max),
        "  ",
        quant.descr.br[htref3$measurements[[meas.name]]$quant],
        "\n", sep="")
  })
