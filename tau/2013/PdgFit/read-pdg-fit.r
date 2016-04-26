#!/usr/bin/env Rscript

##
## read-pdg-fit.r
##
## - read PDG 2015 fit output with no scale factors, s035-fit-ns-all-2016-01-28.fit
## - read HFAG 2014 fit data
## - read HFAG fit with same inputs as PDG 2015 fit, pdgfit-pdg-meas.rdata
##

"Usage:
  read-pdg-fit.r [-l] [--last]

Arguments:

Options:
  -l --last     read last HFAG PDG 2016 fit candidate data
" -> doc

require(docopt, quiet=TRUE)
require(stringr, quiet=TRUE)
require(readr, quiet=TRUE)

source("../../../Common/bin/alucomb2-utils.r")

##
## functions
##

##--- read saved variables inside a list
load.in.list <- function(.file.name) { load(.file.name); as.list(environment()) }

##
## assign list return fields directly to variables
## example: list[a, b] = list(1, 2)
##
list <- structure(NA, class="result")
"[<-.result" <- function(x,...,value) {
  args <- as.list(match.call())
  args <- args[-c(1:2,length(args))]
  length(value) <- length(args)
  for(i in seq(along=args)) {
    a <- args[[i]]
    if(!missing(a)) eval.parent(substitute(a <- v,list(a=a, v=value[[i]])))
  }
  x
}

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
alucomb2.gamma.num.id = function(gamma.name) {
  gamma.name = sub("Unitarity", "Gamma1000", gamma.name, fixed=TRUE)
  gamma.name = sub("GammaAll", "Gamma999", gamma.name, fixed=TRUE)
  rc = str_match(gamma.name, "\\D+(\\d+)by(\\d+)|\\D+(\\d+)")
  if (length(rc) == 0) return(numeric(0))
  num = ifelse(!is.na(rc[,4]), as.numeric(rc[,4])*1000, as.numeric(rc[,2])*1000 + as.numeric(rc[,3]))
  return(num)
}

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
## setup information from a HFAG fit dataset
## - measurement data from input cards data (without HFAG tweaks)
## - quantity nodes, descr
##
get.hfag.input = function(hfag.data) {
  hfag.data$meas.in.gamma = sapply(hfag.data$measurements, function(x) {x$tags[2]})
  hfag.data$meas.in.val = sapply(hfag.data$measurements, function(x) {ifelse(is.null(x$value.orig), x$value, x$value.orig)})
  hfag.data$meas.in.stat = sapply(hfag.data$measurements, function(x) {x$stat})
  hfag.data$meas.in.stat.p = sapply(hfag.data$measurements,
    function(x) {ifelse(attr(x$stat.p, "input") != "", parse_double(attr(x$stat.p, "input")), x$stat)})
  hfag.data$meas.in.stat.n = sapply(hfag.data$measurements,
    function(x) {ifelse(attr(x$stat.n, "input") != "", parse_double(attr(x$stat.n, "input")), -x$stat)})
  hfag.data$meas.in.syst = sapply(hfag.data$measurements, function(x) {ifelse(is.null(x$syst.orig), x$syst, x$syst.orig)})
  hfag.data$meas.in.syst.p = sapply(hfag.data$measurements,
    function(x) {ifelse(attr(x$syst.p, "input") != "",
                        parse_double(attr(x$syst.p, "input")),
                        ifelse(is.null(x$syst.orig), x$syst, x$syst.orig))})
  hfag.data$meas.in.syst.n = sapply(hfag.data$measurements,
    function(x) {ifelse(attr(x$syst.n, "input") != "",
                        parse_double(attr(x$syst.n, "input")),
                        ifelse(is.null(x$syst.orig), -x$syst, -x$syst.orig))})
  hfag.data$meas.in.err = sqrt(hfag.data$meas.in.stat^2 + hfag.data$meas.in.syst^2)

  ##--- HFAG nodes
  quant.node = sapply(hfag.data$combination$quantities, function(quant) {ifelse(is.null(quant$node), "", quant$node)})
  ##--- fix some HFAG-specific nodes
  quant.node["Gamma998"] = "B(unitarity residual)"
  quant.node["GammaAll"] = "B(total)"
  if (!is.na(quant.node["Gamma805"]) && quant.node["Gamma805"] == "")
    quant.node["Gamma805"] = "B(a1(pi- gamma) nu(tau))"
  if (!is.na(quant.node["Gamma910"]) && quant.node["Gamma910"] == "")
    quant.node["Gamma910"] = "B(2pi- pi+ eta(3pi0) nu(tau) (ex. K0))"
  if (!is.na(quant.node["Gamma911"]) && quant.node["Gamma911"] == "")
    quant.node["Gamma911"] = "B(pi- 2pi0 eta(pi+ pi- pi0) nu(tau))"
  if (!is.na(quant.node["Gamma930"]) && quant.node["Gamma930"] == "")
    quant.node["Gamma930"] = "B(2pi- pi+ eta(pi+ pi- pi0) nu(tau) (ex. K0))"
  if (!is.na(quant.node["Gamma944"]) && quant.node["Gamma944"] == "")
    quant.node["Gamma944"] = "B(2pi- pi+ eta(gamma gamma) nu(tau) (ex. K0))"
  if (!is.na(quant.node["Gamma997"]) && quant.node["Gamma997"] == "")
    quant.node["Gamma997"] = "B(a1- -> pi- gamma)"

  quant.gamma = names(quant.node)
  names(quant.gamma) = quant.node[quant.gamma]

  ##--- HFAG modes text description
  quant.descr = sapply(hfag.data$combination$quantities, function(quant) quant$descr)

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

  hfag.data$quant.node = quant.node
  hfag.data$quant.gamma = quant.gamma
  hfag.data$quant.descr = quant.descr
  hfag.data$quant.descr.br = quant.descr.br
  hfag.data$quant.gamma.descr = quant.gamma.descr

  return(hfag.data)
}

##
## get fit result summary
##
get.pdg.summary = function(pdg.data) {
  input.fit.summary = get.section(pdg.data$fit.txt, "^ FINAL RESULTS", "^\\s*$")
  val = str_match(input.fit.summary[2], "FINAL CHISQUARE[=\\s]+(\\d+[.]*\\d*) FOR (\\d+)")[2:3]
  pdg.data$chisq = parse_double(val[1])
  pdg.data$dof = parse_integer(val[2])
  pdg.data$chisq.prob = parse_double(str_match(input.fit.summary[3], "FINAL CONFIDENCE\\D+(\\d+[.]*\\d*)")[2])
  return(pdg.data)
}

##
## get PDG input parameters
##
get.pdg.input.params = function(pdg.data) {
  input.params = get.section(pdg.data$fit.txt, "^ INPUT PARAMETERS", "^\\s*$")

  pdg.data$params.in.count =   parse_number(str_sub(input.params,  1,  7))
  pdg.data$params.in.parcode = str_trim(str_sub(input.params,  8, 16))
  pdg.data$params.in.param = parse_number(str_sub(input.params, 17, 25))
  pdg.data$params.in.seed =  parse_double(str_sub(input.params, 26, 39))
  pdg.data$params.in.descr = str_trim(str_sub(input.params, 40))

  pdg.data$params.in = data.frame(
    count =   parse_number(str_sub(input.params,  1,  7)),
    parcode = str_trim(str_sub(input.params,  8, 16)),
    param = parse_number(str_sub(input.params, 17, 25)),
    seed =  parse_double(str_sub(input.params, 26, 39)),
    descr = str_trim(str_sub(input.params, 40)),
    stringsAsFactors = FALSE
    )

  ##--- set input parameters parcode to preceding row parcode if missing
  pdg.data$params.in.parcode = mapply(
    function(i) {
      parcode = pdg.data$params.in.parcode[i]
      previous.parcodes = pdg.data$params.in.parcode[i:1]
      ifelse(parcode=="",
             previous.parcodes[which("" != previous.parcodes)[1]],
             parcode)
    }, 1:nrow(pdg.data$params.in))

  ##--- add label summarizing parcode and param
  pdg.data$params.in.plab = paste0(pdg.data$params.in.parcode, "_p", pdg.data$params.in.param)

  ##--- change description to canonical text description form
  pdg.data$params.in.descr = gsub("tau- --> (.*)", "G(\\1)/G(total)", pdg.data$params.in.descr)

  return(pdg.data)
}

##
## get PDG fit results for parameters
##
get.pdg.fitted.params = function(pdg.data) {
  input.params.fit = get.section(pdg.data$fit.txt, c("^ FINAL RESULTS", "^ RESULTS FOR PARAMETERS"), "^\\s*$")

  pdg.data$params.fit = data.frame(
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
  pdg.data$params.fit$descr = gsub("tau- --> (.*)", "G(\\1)/G(total)", pdg.data$params.fit$descr)

  return(pdg.data)
}

##
## get PDG input nodes
##
get.pdg.input.nodes = function(pdg.data) {
  input.nodes = get.section(pdg.data$fit.txt, "^ INPUT NODES", "^\\s*$")

  pdg.data$nodes.in.node.count = parse_number(str_sub(input.nodes,  1,  7))
  pdg.data$nodes.in.node = str_trim(str_sub(input.nodes,  8, 15))
  pdg.data$nodes.in.eq = str_trim(str_sub(input.nodes, 16, 17))
  pdg.data$nodes.in.adj =  parse_double(str_sub(input.nodes, 18, 21))
  pdg.data$nodes.in.un = str_trim(str_sub(input.nodes, 22, 27))
  pdg.data$nodes.in.count = parse_number(str_sub(input.nodes, 28, 32))
  pdg.data$nodes.in.parcode = str_trim(str_sub(input.nodes, 33, 42))
  pdg.data$nodes.in.param = parse_number(str_sub(input.nodes, 43, 47))
  pdg.data$nodes.in.sum = parse_number(str_sub(input.nodes, 51, 53))
  pdg.data$nodes.in.coeff = parse_double(str_sub(input.nodes, 55, 65))
  pdg.data$nodes.in.coeff.extra = parse_double(str_sub(input.nodes, 66, 76))
  pdg.data$nodes.in.descr = str_trim(str_sub(input.nodes, 80))

  pdg.data$nodes.in.plab = paste0(pdg.data$nodes.in.parcode, "_p", pdg.data$nodes.in.param)

  ##--- get begin of node definitions
  pdg.data$node.in.def.beg = which(!is.na(pdg.data$nodes.in.node.count))
  ##--- get end of node definitions as beginning of next node minus one
  pdg.data$node.in.def.end = c(pdg.data$node.in.def.beg[-1], length(pdg.data$nodes.in.node.count)+1) - 1

  ##--- for each node def, store description and node
  pdg.data$quant.descr = pdg.data$nodes.in.descr[pdg.data$node.in.def.beg]
  pdg.data$quant.node = pdg.data$nodes.in.node[pdg.data$node.in.def.beg]

  ##--- for the nodes that correspond to a fit parameter, store parameter label else empty string
  pdg.data$quant.plab = ifelse(pdg.data$node.in.def.beg != pdg.data$node.in.def.end, "",
    paste0(pdg.data$nodes.in.parcode[pdg.data$node.in.def.beg], "_p", pdg.data$nodes.in.param))

  ##--- set labels for conversion from descr to node
  names(pdg.data$quant.node) = pdg.data$quant.descr
  ##--- set labels for conversion from node to descr
  names(pdg.data$quant.descr) = pdg.data$quant.node
  ##--- set labels for conversion from node to parameter label
  names(pdg.data$quant.plab) = pdg.data$quant.node

  ##
  ## get PDG constraints with numeric coefficients
  ##
  pdg.data$constr.numeric = mapply(
    function(def.beg, def.end) {
      sel.above = which(pdg.data$nodes.in.sum[def.beg:def.end] == 1)
      sel.below = which(pdg.data$nodes.in.sum[def.beg:def.end] == 2)
      ##--- what node is defined
      node.defined = pdg.data$nodes.in.node[def.beg]
      ##--- what parameters and coeff. above fraction line
      coeff.above = pdg.data$nodes.in.coeff[def.beg:def.end][sel.above]
      names(coeff.above) = pdg.data$nodes.in.plab[def.beg:def.end][sel.above]
      ##--- what parameters and coeff. below fraction line
      coeff.below = pdg.data$nodes.in.coeff[def.beg:def.end][sel.below]
      names(coeff.below) = pdg.data$nodes.in.plab[def.beg:def.end][sel.below]
      ##--- return node defined, coefficients and parameters in numerator and denominator
      list(node=node.defined, above=coeff.above, below=coeff.below)
    }, pdg.data$node.in.def.beg, pdg.data$node.in.def.end, SIMPLIFY=FALSE)

  return(pdg.data)
}

##
## get PDG fit results for nodes
##
get.pdg.fitted.nodes = function(pdg.data) {

  input.nodes.fit = get.section(pdg.data$fit.txt, c("^ FINAL RESULTS", "^ RESULTS FOR NODES"), "^\\s*$")

  pdg.data$nodes.fit = data.frame(
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

  return(pdg.data)
}

##
## get input measurements
##
get.pdg.input.meas = function(pdg.data) {
  ##
  ## input measurements
  ##
  input.meas = get.section(pdg.data$fit.txt, "^ INPUT MEASUREMENTS", "^\\s*$")

  ##--- initially, put stat and syst errors in positive errors
  pdg.data$meas.in = data.frame(
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
  pdg.data$meas.in$statn = pdg.data$meas.in$statp
  pdg.data$meas.in$systn = pdg.data$meas.in$systp

  ##--- identify where are asymmetric errors
  stat.asymm.p = which(pdg.data$meas.in$stat.sign == "+")
  stat.asymm.n = which(pdg.data$meas.in$stat.sign == "-")
  if (any(stat.asymm.n != stat.asymm.p+1)) {
    stop("asymmetric errors do not respect having first + then -")
  }

  pdg.data$meas.in[stat.asymm.p, "statn"] = pdg.data$meas.in[stat.asymm.n, "statn"]
  pdg.data$meas.in = pdg.data$meas.in[-stat.asymm.n, ]

  syst.asymm.p = which(pdg.data$meas.in$syst.sign == "+")
  syst.asymm.n = which(pdg.data$meas.in$syst.sign == "-")
  if (any(syst.asymm.n != syst.asymm.p+1)) {
    stop("asymmetric errors do not respect having first + then -")
  }

  ##--- set negative error on line where the positive error was set
  pdg.data$meas.in[syst.asymm.p, "systn"] = pdg.data$meas.in[syst.asymm.n, "systn"]
  ##--- remove lines where the negative error was set
  pdg.data$meas.in = pdg.data$meas.in[-syst.asymm.n, ]

  ##--- compute squared root averages of positive and negative errors
  pdg.data$meas.in$stat = sqrt( (pdg.data$meas.in$statp^2 + pdg.data$meas.in$statn^2)/2 )
  pdg.data$meas.in$syst = sqrt( (pdg.data$meas.in$systp^2 + pdg.data$meas.in$systn^2)/2 )
  ##--- total error
  pdg.data$meas.in$err = sqrt( pdg.data$meas.in$stat^2 + pdg.data$meas.in$syst^2 )

  ##--- set node on measurements that miss it using the last defined one
  pdg.data$meas.in$node = unname(
    mapply(
      function(i) {
        node = pdg.data$meas.in$node[i]
        previous.nodes = pdg.data$meas.in$node[i:1]
        ifelse(node=="",
               previous.nodes[which("" != previous.nodes)[1]],
               node)
      }, 1:nrow(pdg.data$meas.in)))

  ##
  ## input measurements correlations
  ##

  input.meas.corr.txt =
    get.section(pdg.data$fit.txt, "^ INPUT CORRELATIONS BETWEEN MEASUREMENTS", "^ FIRST FIT", skip.beg=1, skip.end=2)

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

  pdg.data$meas.in$stat.corr = lapply(1:length(pdg.data$meas.in$node),
    function(ii.meas) {
      which.rows = which(ii.meas == sapply(input.meas.corr.list, function(el) el["row"]))
      meas.corr.rows = lapply(which.rows, function(ii.el) { input.meas.corr.list[[ii.el]][c("col", "corr")]})
      which.cols = which(ii.meas == sapply(input.meas.corr.list, function(el) el["col"]))
      meas.corr.cols = lapply(which.cols, function(ii.el) { input.meas.corr.list[[ii.el]][c("row", "corr")]})
      meas.corr = c(meas.corr.rows, meas.corr.cols)
    })

  return(pdg.data)
}

##
## compare input measurements in PDG and HFAG
##
cmp.meas.pdg.hfag = function(pdg, hfag) {
  cmp = list()

  ##
  ## compute distance between measurements in pdg vs. hfag
  ## (quadratic sum of value and error asymmetry)
  ##
  val.cmp = outer(pdg$meas.in$val, hfag$meas.in.val, rel.diff)
  err.cmp = outer(pdg$meas.in$err, hfag$meas.in.err, rel.diff)
  tot.cmp = sqrt(val.cmp^2 + err.cmp^2)

  ##--- for each PDG measurement get closest HFAG measurement
  cmp$pdg.match = sapply(seq(1, nrow(val.cmp)), function(i) {order(tot.cmp[i, ])[1]})

  ##--- there must not be any duplication
  if (any(duplicated(cmp$pdg.match))) {
    stop("more than one HFAG measurement matches a single PDG measurement")
  }

  ##--- compute best matched measurements, starting with largest discrepancies
  cmp$pdg.meas.name = names(hfag$meas.in.val[cmp$pdg.match])

  cmp$cmp = data.frame(
    pdg = paste(pdg$meas.in$author, pdg$meas.in$year),
    node = pdg$meas.in$node,
    hfag = cmp$pdg.meas.name,
    gamma = hfag$meas.in.gamma[cmp$pdg.match],

    diff = sapply(1:nrow(tot.cmp), function(i) {tot.cmp[i, cmp$pdg.match[i]]}),

    pdg.val = pdg$meas.in$val,
    hfag.val = hfag$meas.in.val[cmp$pdg.match],

    pdg.err = pdg$meas.in$err,
    hfag.err = hfag$meas.in.err[cmp$pdg.match],

    pdg.stat = pdg$meas.in$stat,
    hfag.stat = hfag$meas.in.stat[cmp$pdg.match],

    pdg.syst = pdg$meas.in$syst,
    hfag.syst = hfag$meas.in.syst[cmp$pdg.match],

    pdg.statp = pdg$meas.in$statp,
    hfag.statp = hfag$meas.in.stat.p[cmp$pdg.match],

    pdg.statn = -pdg$meas.in$statn,
    hfag.statn = hfag$meas.in.stat.n[cmp$pdg.match],

    pdg.systp = pdg$meas.in$systp,
    hfag.systp = hfag$meas.in.syst.p[cmp$pdg.match],

    pdg.systn = -pdg$meas.in$systn,
    hfag.systn = hfag$meas.in.syst.n[cmp$pdg.match],

    stringsAsFactors=FALSE
    )
  rownames(cmp) = NULL

  ##--- sort with largest discrepancies first
  cmp$cmp = cmp$cmp[order(cmp$cmp$diff, decreasing=TRUE), ]

  return(cmp)
}

##
## print input measurements mismatches between PDG and HFAG
##
pr.mism.meas.pdg.hfag = function(cmp, threshold=0) {
  above.threshold = which(cmp$cmp$diff >= threshold)
  if (length(above.threshold) > 0) {
    rc = lapply(split(cmp$cmp[above.threshold,], 1:length(above.threshold)),
      function(cmp) {
        cat(cmp$node, cmp$pdg, cmp$gamma)
        cat("\n")
        ## print(cmp)
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
}

##
## print mismatches in node description
##
pr.mism.descr.pdg.hfag = function(pdg, hfag) {
  quant.hfag.descr = hfag$quant.descr[hfag$quant.gamma[pdg$quant.node]]
  non.matching = gsub("\\s+", "", quant.hfag.descr) != gsub("\\s+", "", pdg$quant.descr)

  rc = mapply(
    function(descr.pdg, descr.hfag, gamma.hfag) {
      cat(gamma.hfag, "\n  pdg  ", descr.pdg, "\n  hfag ", descr.hfag, "\n", sep="")
    },
    pdg$quant.descr[non.matching],
    quant.hfag.descr[non.matching],
    names(quant.hfag.descr[non.matching]))
}

##
## return list of HFAG gammas corresponding to PDG nodes
##
hfag.gamma.pdg.nodes = function(hfag, pdg) {
  hfag$quant.gamma[pdg$quant.node]
}

##
## return list of HFAG gammas corresponding to PDG parameters
##
hfag.gamma.pdg.params = function(hfag, pdg) {
  hfag$quant.gamma[pdg$params.in.node]
}

##
## return list of HFAG gammas corresponding to PDG nodes that are derived from parameters
##
hfag.gamma.pdg.derived = function(hfag, pdg) {
  setdiff(hfag$quant.gamma[pdg$quant.node], hfag$quant.gamma[pdg$params.in.node])
}

##
## find out symbolic expressions for the coefficients of the PDG constraints
## the definitions of the parameters are taken from an HFAG data file
##
get.pdg.constr.params.symbolic = function(pdg, hfag) {
  hfag.params = lapply(hfag$combination$params, function(x) unname(x["value"]))

  ##
  ## parameter expressions checked for replacing of PDG coefficients in nodes
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
  pdg$params.expr.val = eval.all(esub.str.parse(params.expr.txt, hfag.params))

  ##--- find expression best approximating PDG coefficients in nodes
  mism = sapply(pdg$nodes.in.coeff,
    function(coeff) {
      ##+++ patch to get 0.09 identified as omega -> pi0 gamma
      if (coeff == 0.09) {coeff = eval(substitute(BR_om_pizgamma, hfag.params))}
      ##+++ patch to get 0.3431 identified as 1/2 * KS -> pi+pi-
      if (coeff == 0.3431) {coeff = eval(substitute(BR_KS_pimpip/2, hfag.params))}
      mism = abs(coeff - pdg$params.expr.val)*2 / (abs(coeff) + abs(pdg$params.expr.val))
      mism = mism[order(mism)[1]]
    })

  ##--- store expressions for PDG coefficients in nodes
  pdg$nodes.in.coeff.expr.txt = names(mism)
  ##--- store mismatch of found symbolic expression
  pdg$nodes.in.coeff.expr.mism = mism

  ##
  ## get PDG constraints with symbolic coefficients
  ##
  pdg$constr = mapply(
    function(def.beg, def.end) {
      sel.above = which(pdg$nodes.in.sum[def.beg:def.end] == 1)
      sel.below = which(pdg$nodes.in.sum[def.beg:def.end] == 2)
      ##--- what node is defined
      node.defined = pdg$nodes.in.node[def.beg]
      ##--- what parameters and coeff. above fraction line
      coeff.above = pdg$nodes.in.coeff.expr.txt[def.beg:def.end][sel.above]
      names(coeff.above) = pdg$nodes.in.plab[def.beg:def.end][sel.above]
      ##--- what parameters and coeff. below fraction line
      coeff.below = pdg$nodes.in.coeff.expr.txt[def.beg:def.end][sel.below]
      names(coeff.below) = pdg$nodes.in.plab[def.beg:def.end][sel.below]
      ##--- return node defined, coefficients and parameters in numerator and denominator
      list(node=node.defined, above=coeff.above, below=coeff.below)
    }, pdg$node.in.def.beg, pdg$node.in.def.end, SIMPLIFY=FALSE)

  return(pdg)
}

##
## log how PDG constraint coefficients are interpreted as symbolic expressions
##
pr.pdg.constr.params.symbolic = function(pdg) {
  cat("##\n")
  cat("## Parameter expressions for PDG nodes definitions coefficients\n")
  cat("## - 0.09 forced to be matched with B(omega -> pi0 gamma)\n")
  cat("## - 0.3431 forced to be matched with 1/2*B(KS -> pi+pi-)\n")
  cat("##\n")
  coeff.unique.indices = which(!duplicated(pdg$nodes.in.coeff))
  rc.unique = pdg$nodes.in.coeff.expr.mism[coeff.unique.indices]
  rc.unique.order = order(rc.unique, decreasing=TRUE)
  coeff.unique.indices = coeff.unique.indices[rc.unique.order]
  cat("      coeff   closest expression                  expression value mismatch\n\n")
  rc2 = lapply(coeff.unique.indices,
    function(coeff.ii) {
      cat(paste0(sprintf("%11.6g", pdg$nodes.in.coeff[coeff.ii]), " = ",
                 format(pdg$nodes.in.coeff.expr.txt[coeff.ii], width=40),
                 " ", sprintf("%11.4g", pdg$params.expr.val[pdg$nodes.in.coeff.expr.txt[coeff.ii]]),
                 " ", sprintf("%#7.4f%%", 100*pdg$nodes.in.coeff.expr.mism[coeff.ii]), "\n"))
    })
}

##
## list HFAG fit quantities matching PDG input nodes
##
pr.hfag.gamma.pdg.nodes = function(pdg, hfag) {
  quant.gamma.combine = c(hfag.gamma.pdg.params(hfag, pdg), hfag.gamma.pdg.derived(hfag, pdg))
  quant.gamma.combine.maxlen = max(nchar(quant.gamma.combine))
  gamma.per.line = trunc((80-2) / (quant.gamma.combine.maxlen+2))
  range.beg = seq(1, length(quant.gamma.combine), by=gamma.per.line)
  range.end = ifelse(range.beg+gamma.per.line-1 > length(quant.gamma.combine),
    length(quant.gamma.combine), range.beg+gamma.per.line-1)
  cat("COMBINE\n")
  rc = mapply(
    function(beg, end) {
      cat("  ", paste(format(quant.gamma.combine[beg:end], width=quant.gamma.combine.maxlen), collapse=" "), "\n", sep="")
    }, range.beg, range.end)
  cat("\n")
}

##
## save HFAG measurement structure into a file
##
save.hfag.measurements = function(measurements, hfag, file.dir=".") {
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
    meas.fname = meas.fname[order(alucomb2.gamma.num.id(meas.fname.gamma), meas.fname.pub)]
    for(meas in meas.fname) {
      alucomb2.print.meas(meas, hfag$combination$quantities)
      cat("\n")
    }
    sink()
    close(meas.fh)
    cat(paste0("file ", exp.fname, " created\n"))
  }
}

##
## save PDG measurements in HFAG format
##
save.pdg.measurements = function(pdg, hfag, pdg.match, file.dir=".") {
  ##--- collect PDG measurements in HFAG format
  pdg.measurements = lapply(1:length(pdg$meas.in$node),
    function(ii.meas) {
      corr.terms.tot = sapply(unname(pdg$meas.in$stat.corr[[ii.meas]]),
        function(x) setNames(x[2]/100, paste(hfag$measurements[[pdg.match[x[1]]]]$tags, collapse=".")))
      if (length(corr.terms.tot) > 0) {
        corr.terms.tot = corr.terms.tot[order(alucomb2.gamma.num.id(
          str_match(names(corr.terms.tot), "[^.]+[.](Gamma(\\d+)(by|)\\d*)[.]")[,2]))]
      }
      meas = list(
        value = pdg$meas.in$val[ii.meas],

        stat = pdg$meas.in$stat[ii.meas],
        stat.p = pdg$meas.in$statp[ii.meas],
        stat.n = -pdg$meas.in$statn[ii.meas],

        syst = pdg$meas.in$syst[ii.meas],
        syst.p = pdg$meas.in$systp[ii.meas],
        syst.n = -pdg$meas.in$systn[ii.meas],

        quant = hfag$measurements[[pdg.match[ii.meas]]]$quant,
        tags = hfag$measurements[[pdg.match[ii.meas]]]$tags,
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

  save.hfag.measurements(pdg.measurements, hfag, file.dir=file.dir)
}

##
## convert PDG costraints to text
##
pdg.constr.to.text = function(constr.list) {
  rc = sapply(constr.list,
    function(constr) {
      if (length(constr$above) == 0) {
        above.txt = "1"
      } else {
        above.txt = paste0(names(constr$above), "*", sub("([-].*)", "(\\1)", constr$above))
        ##--- simplify multiplication by 1 and -1
        above.txt = gsub("*1", "", above.txt, fixed=TRUE)
        above.txt = gsub("*(-1)", "", above.txt, fixed=TRUE)
        above.txt = paste0(above.txt, collapse=" + ")
      }
      if (length(constr$below) == 0) {
        full.txt = above.txt
      } else {
        if (length(constr$above) > 1) {
          above.txt = paste0("(", above.txt, ")")
        }
        below.txt = paste0(names(constr$below), "*", sub("([-].*)", "(\\1)", constr$below))
        ##--- simplify multiplication by 1 and -1
        below.txt = gsub("*1", "", below.txt, fixed=TRUE)
        below.txt = gsub("*(-1)", "", below.txt, fixed=TRUE)
        below.txt = paste0(below.txt, collapse=" + ")
        if (length(constr$below) > 1) {
          below.txt = paste0("(", below.txt, ")")
        }
        full.txt = paste0(above.txt, " / ", below.txt)
      }
      paste0(constr$node, " = ", full.txt)
    })}

##
## convert PDG constraints to text using HFAG format
##
pdg.constr.to.text.hfag = function(pdg, hfag) {
  constr.hfag = lapply(pdg$constr,
    function(constr) {
      constr$node = unname(hfag$quant.gamma[constr$node])
      names(constr$above) = hfag$quant.gamma[pdg$params.in.node[names(constr$above)]]
      names(constr$below) = hfag$quant.gamma[pdg$params.in.node[names(constr$below)]]
      constr
    })
  pdg.constr.to.text(constr.hfag)
}

##
## convert list of constraints to non-trivial list
## - remove trivia constrains of the form GammaXXX = GammaXXX
## - format as NLCONSTRAINT HFAG cards
##
pdg.constr.convert.nontrivial = function(constr.hfag.txt) {
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

  return(constr.hfag.txt.nontrivial)
}

##
## convert PDG constraints to text using human format
##
pdg.constr.to.text.human = function(pdg, hfag) {
  constr.human = lapply(pdg$constr,
    function(constr) {
      constr$node = unname(hfag$quant.descr.br[hfag$quant.gamma[constr$node]])
      names(constr$above) = hfag$quant.descr.br[hfag$quant.gamma[pdg$params.in.node[names(constr$above)]]]
      names(constr$below) = hfag$quant.descr.br[hfag$quant.gamma[pdg$params.in.node[names(constr$below)]]]
      constr
    })
  pdg.constr.to.text(constr.human)
}

##
## get constraints text espressions from HFAG
##
hfag.constr.gamma.expr = function(hfag) {
  hfag.constr.str.expr = paste0(
    "NLCONSTRAINT ",
    names(hfag$combination$constr.nl.str.expr), " ",
    as.character(hfag$combination$constr.nl.str.val), " \"",
    hfag$combination$constr.nl.str.expr, "\""
    )
  hfag.constr.str.expr = gsub("Unitarity ", "Unitarity.c ", hfag.constr.str.expr)
  ##--- otbain equation in simplest form
  hfag.constr.str.expr = gsub("^NLCONSTRAINT .*[.]c (\\d+) \"(.*)\"", "\\1 = \\2", hfag.constr.str.expr)
  hfag.constr.str.expr = gsub("0 = -(\\S+) [+] \\((.*)\\)", "\\1 = \\2", hfag.constr.str.expr)
  hfag.constr.str.expr = gsub("0 = -(\\S+) [+] (.*)", "\\1 = \\2", hfag.constr.str.expr)
  ##--- get list of nodes in right part of equation
  hfag.constr.gammas = str_extract_all(hfag.constr.str.expr, "(Gamma\\d+by\\d+|Gamma\\d+)")
  ##--- form string with sorted right side nodes
  rc = sapply(hfag.constr.gammas,
    function(gammas) {
      paste(sort(gammas[-1]), collapse=" ")
    })
  names(rc) = sapply(hfag.constr.gammas, function(constr) {constr[1]})
  hfag.constr.gammas = sort(rc)
  names(hfag.constr.str.expr) = names(rc)

  return(list(gammas=hfag.constr.gammas, expr=hfag.constr.str.expr))
}

##
## convert PDG constraint in HFAG format to simplest form
##
pdg.constr.str.expr = function(constr.pdg.hfag) {
  ##--- otbain equation in simplest form
  pdg.constr.str.expr = gsub("^NLCONSTRAINT .*[.]c (\\d+) \"(.*)\"", "\\1 = \\2", constr.pdg.hfag)
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

  return(list(gammas=pdg.constr.gammas, expr=pdg.constr.str.expr))
}

##
## take two lists of constraint equations in text format
## for each item of the two lists
## - get all the HFAG Gamma symbols
##   first left side, then the right side, sorted by Gamma number
## - return the array of Gamma symbols, left first
## if one list is empty, use just the first one
##
get.gamma.list.constr = function(expr.list.1, expr.list.2=NULL) {
  if (!is.null(expr.list.2)) {
    rc = mapply(
      function(elem.1, elem.2) {
        gammas.1 = str_extract_all(elem.1, "(Gamma\\d+by\\d+|Gamma\\d+)")[[1]]
        gammas.2 = str_extract_all(elem.2, "(Gamma\\d+by\\d+|Gamma\\d+)")[[1]]
        right = unique(c(gammas.1[-1], gammas.2[-1]))
        right = right[order(alucomb2.gamma.num.id(right))]
        left = unique(c(gammas.1[1], gammas.2[1]))
        left = left[order(alucomb2.gamma.num.id(left))]
        gammas = c(left, right)
      }, expr.list.1, expr.list.2, SIMPLIFY=FALSE)
  } else {
    rc = mapply(
      function(elem.1) {
        gammas = str_extract_all(elem.1, "(Gamma\\d+by\\d+|Gamma\\d+)")[[1]]
        gammas.right = unique(gammas[-1])
        gammas.right = gammas.right[order(alucomb2.gamma.num.id(gammas.right))]
        gammas = c(gammas[1], gammas.right)
      }, expr.list.1, SIMPLIFY=FALSE)
  }
  rc
}

##
## take a list where each element is a list of Gamma symbols
## for each element of the list
## - assemble text equating the Gamma symbols to their description
##
get.hfag.constr.comments = function(hfag, expr.list.1, expr.list.2=NULL) {
  gammas = get.gamma.list.constr(expr.list.1, expr.list.2)
  rc = sapply(gammas,
    function(gamma.list) {
      rc = paste(gamma.list, "=", hfag$quant.descr.br[gamma.list])
      rc = paste(rc, collapse="\n# ")
      paste0("#\n# ", rc, "\n#\n")
    })
  rc
}

##
## replace Gamma HFAG notation with description (not used)
##
hfag.gamma.to.node = function(str.expr, quant.node) {
  quant.node = ifelse(quant.node == "", "NA", quant.node)
  gsub("", "", esub.str(str.expr, lapply(quant.node, as.symbol)))
}

##
## substitute parameters with numeric values
##
hfag.subst.param = function(str.expr, hfag) {
  params.val = lapply(hfag$combination$params, function(x) unname(x["value"]))
  esub.str(str.expr, params.val)
}

##
## compare PDG constraints with HFAG constraints
##
cmp.pgd.hfag.constr = function(pdg, hfag, pdg.label="PDG", hfag.label="HFAG") {
  ##--- get HFAG constraint info
  list[hfag.constr.gammas, hfag.constr.str.expr] = hfag.constr.gamma.expr(hfag)
  ##--- get PDG constraint info
  list[pdg.constr.gammas, pdg.constr.str.expr] =
    pdg.constr.str.expr(
      pdg.constr.convert.nontrivial(
        pdg.constr.to.text.hfag(pdg, hfag)))
  ##
  ## compare for the same left side the sorted list of nodes in the right side
  ##
  ##--- constraints that have same left side
  constr.left.common = intersect(names(hfag.constr.gammas), names(pdg.constr.gammas))
  constr.non.equal = names(which(hfag.constr.gammas[constr.left.common] != pdg.constr.gammas[constr.left.common]))
  constr.equal = names(which(hfag.constr.gammas[constr.left.common] == pdg.constr.gammas[constr.left.common]))

  ##--- in PDG but not HFAG
  constr.pdg.not.hfag = setdiff(names(pdg.constr.gammas), names(hfag.constr.gammas))
  ##--- in HFAG but not PDG
  constr.hfag.not.pdg = setdiff(names(hfag.constr.gammas), names(pdg.constr.gammas))

  cat("##\n")
  cat("## Constraints that are the same on", pdg.label, "and", hfag.label)
  cat("\n")
  cat("## PDG first, HFAG second\n")
  cat("##\n")

  cat(
    paste0(
      paste0(
        get.hfag.constr.comments(hfag, hfag.constr.str.expr[constr.equal]),
        pdg.constr.str.expr[constr.equal], "\n",
        hfag.constr.str.expr[constr.equal], "\n",
        hfag.gamma.to.node(pdg.constr.str.expr[constr.equal], hfag$quant.node), "\n",
        hfag.gamma.to.node(hfag.constr.str.expr[constr.equal], hfag$quant.node),
        sep="\n"),
      collapse="\n"))
  cat("\n")

  ##--- flag array indicating what constraints HFAG uses
  hfag.constr.flag = hfag$combination$constr.all.nl | hfag$combination$constr.all.lin
  hfag.constr.used = names(hfag.constr.str.expr[hfag.constr.flag])

  ##--- non-equal AND used in HFAG
  constr.non.equal.used = intersect(constr.non.equal, hfag.constr.used)
  ##--- non-equal AND NOT used in HFAG
  constr.non.equal.not.used = setdiff(constr.non.equal, hfag.constr.used)

  cat("##\n")
  cat("## Constraints that are NOT the same on", pdg.label, "and", hfag.label)
  cat(", used in HFAG\n")
  cat("## PDG first, HFAG second\n")
  cat("##\n")

  cat(
    paste0(
      paste0(
        get.hfag.constr.comments(hfag, pdg.constr.str.expr[constr.non.equal.used], hfag.constr.str.expr[constr.non.equal.used]),
        pdg.constr.str.expr[constr.non.equal.used], "\n",
        hfag.constr.str.expr[constr.non.equal.used], "\n",
        hfag.gamma.to.node(pdg.constr.str.expr[constr.non.equal.used], hfag$quant.node), "\n",
        hfag.gamma.to.node(hfag.constr.str.expr[constr.non.equal.used], hfag$quant.node), "\n",
        hfag.subst.param(hfag.gamma.to.node(pdg.constr.str.expr[constr.non.equal.used], hfag$quant.node), hfag), "\n",
        hfag.subst.param(hfag.gamma.to.node(hfag.constr.str.expr[constr.non.equal.used], hfag$quant.node), hfag),
        sep="\n"),
      collapse="\n"))
  cat("\n")

  cat("##\n")
  cat("## Constraints that are NOT the same on", pdg.label, "and", hfag.label)
  cat(", NOT used in HFAG\n")
  cat("## PDG first, HFAG second\n")
  cat("##\n")

  cat(
    paste0(
      paste0(
        get.hfag.constr.comments(hfag, pdg.constr.str.expr[constr.non.equal.not.used], hfag.constr.str.expr[constr.non.equal.not.used]),
        pdg.constr.str.expr[constr.non.equal.not.used], "\n",
        hfag.constr.str.expr[constr.non.equal.not.used], "\n",
        hfag.gamma.to.node(pdg.constr.str.expr[constr.non.equal.not.used], hfag$quant.node), "\n",
        hfag.gamma.to.node(hfag.constr.str.expr[constr.non.equal.not.used], hfag$quant.node), "\n",
        hfag.subst.param(hfag.gamma.to.node(pdg.constr.str.expr[constr.non.equal.not.used], hfag$quant.node), hfag), "\n",
        hfag.subst.param(hfag.gamma.to.node(hfag.constr.str.expr[constr.non.equal.not.used], hfag$quant.node), hfag),
        sep="\n"),
      collapse="\n"))
  cat("\n")

  cat("##\n")
  cat("## Constraints in", pdg.label, "but not in", hfag.label)
  cat("\n")
  cat("##\n")

  cat(
    paste0(
      paste0(
        get.hfag.constr.comments(hfag, pdg.constr.str.expr[constr.pdg.not.hfag]),
        pdg.constr.str.expr[constr.pdg.not.hfag], "\n",
        hfag.gamma.to.node(pdg.constr.str.expr[constr.pdg.not.hfag], hfag$quant.node), "\n",
        hfag.subst.param(hfag.gamma.to.node(pdg.constr.str.expr[constr.pdg.not.hfag], hfag$quant.node), hfag),
        sep="\n"),
      collapse="\n"))
  cat("\n")

  cat("##\n")
  cat("## Constraints in", hfag.label, "but not in", pdg.label)
  cat("\n")
  cat("##\n")

  cat(
    paste0(
      paste0(
        get.hfag.constr.comments(hfag, hfag.constr.str.expr[constr.hfag.not.pdg]),
        hfag.constr.str.expr[constr.hfag.not.pdg], "\n",
        hfag.gamma.to.node(hfag.constr.str.expr[constr.hfag.not.pdg], hfag$quant.node), "\n",
        hfag.subst.param(hfag.gamma.to.node(hfag.constr.str.expr[constr.hfag.not.pdg], hfag$quant.node), hfag),
        sep="\n"),
      collapse="\n"))
  cat("\n")

  return(invisible())
}

##
## compare PDG vs. HFAG fit results
##
cmp.pdg.hfag.fit.results = function(pdg, hfag, pdg.label="PDG", hfag.label="HFAG") {
  ##--- get HFAG results in PDG structure with fit results
  pdg$nodes.fit = within(pdg$nodes.fit,
    {
      hfag.quant.val = hfag$quant.val[hfag$quant.gamma[node]]
      hfag.quant.err = hfag$quant.err[hfag$quant.gamma[node]]
      mism.val = (val - hfag.quant.val) / sqrt((err^2+hfag.quant.err^2)/2)
      mism.err = (err - hfag.quant.err) / sqrt((err^2+hfag.quant.err^2)/2)
      mism = sqrt(mism.val^2 + mism.err^2)
    })

  cat("##\n")
  cat("##", pdg.label, "fit summary\n")
  cat("##\n")

  cat("chisq / dof =", pdg$chisq, "/", pdg$dof, " ( Prob = ", pdg$chisq.prob, ")\n")
  cat("measurements:", nrow(pdg$meas.in),
      "fit parameters:", nrow(pdg$params.in),
      "fit nodes:", nrow(pdg$nodes.fit),
      "\n")

  cat("##\n")
  cat("##", hfag.label, "fit summary\n")
  cat("##\n")

  cat("chisq / dof =", hfag$chisq, "/", hfag$dof, " ( Prob = ", hfag$chisq.prob, ")\n")
  cat("measurements:", length(hfag$meas.val),
      "fit quantities:", length(hfag$quant.val),
      "constraints:", hfag$constr.num,
      "\n")

  cat("##\n")
  cat("## Fit results", pdg.label, "vs.", hfag.label, "mismatches\n")
  cat("## - first line:\n")
  cat("##   - (PDG_value - HFAG_value) / <error on value> in percent\n")
  cat("##   - (PDG_error - HFAG_error) / <error on value> in percent\n")
  cat("##   - quantity description\n")
  cat("## - second line: PDG  fit value, error, asymmetric errors if existing\n")
  cat("## - third line:  HFAG fit value, error\n")
  cat("##\n")
  rc = lapply(split(pdg$nodes.fit, 1:nrow(pdg$nodes.fit))[order(pdg$nodes.fit$mism, decreasing=TRUE)],
    function(cmp.fit) {
      line.1 = paste(
        paste(sprintf("%12s", c(cmp.fit$node, hfag$quant.gamma[cmp.fit$node])), collapse=" "),
        paste(sprintf("%+12.4f%%", 100*c(cmp.fit$mism.val, cmp.fit$mism.err)), collapse=" "),
        paste0(" ", cmp.fit$descr),
        sep=" ")
      line.2 = paste(
        sprintf("%25s", "PDG"),
        paste(sprintf("%13.6e", c(cmp.fit$val, cmp.fit$err)), collapse=" "),
        ifelse(is.na(cmp.fit$errp), "", paste(sprintf("%+14.6e", c(cmp.fit$errp, -cmp.fit$errn)), collapse=" ")),
        sep=" ")
      line.3 =  paste(
        sprintf("%25s", "HFAG"),
        paste(sprintf("%13.6e", c(cmp.fit$hfag.quant.val, cmp.fit$hfag.quant.err)), collapse=" "),
        sep=" ")
      cat(paste(line.1, line.2, line.3, "", sep="\n"))
    })
}

##
## list measurements usage mismatches PDG vs. HFAG
##
pr.meas.use.mism.pdg.hfag = function(pdg.meas.names, hfag, pdg.label="PDG", hfag.label="HFAG") {
  ##
  ## list measurements used in HFAG fit but not PDG
  ## HFAG fit without Belle inclusive K0, no hcorr
  ##
  hfag.meas.not.in.pdg = setdiff(names(hfag$meas.val), pdg.meas.names)
  meas.name.width.max = max(nchar(hfag.meas.not.in.pdg))
  cat("##\n")
  cat("## measurements used in", hfag.label, "but not in", pdg.label, "fit\n")
  cat("##\n")

  rc = lapply(hfag.meas.not.in.pdg,
    function(meas.name) {
      cat(format(meas.name, width=meas.name.width.max),
          "  ",
          hfag$quant.descr.br[hfag$measurements[[meas.name]]$quant],
          "\n", sep="")
    })

  ##
  ## list measurements used in PDG fit but not HFAG
  ## HFAG fit without Belle inclusive K0, no hcorr
  ##
  pdg.meas.not.in.hfag = setdiff(pdg.meas.names, names(hfag$meas.val))
  meas.name.width.max = max(nchar(pdg.meas.not.in.hfag))

  cat("##\n")
  cat("## measurements used in", pdg.label, "fit but not in", hfag.label, "fit\n")
  cat("##\n")

  rc = lapply(pdg.meas.not.in.hfag,
    function(meas.name) {
      cat(format(meas.name, width=meas.name.width.max),
          "  ",
          hfag$quant.descr.br[hfag$measurements[[meas.name]]$quant],
          "\n", sep="")
    })
}

##
## main code
##

opts <- docopt(doc)
## opts$last

##
## get info to translate from nodes to desig (when possible)
##
pdginfo = load.in.list("pdginfo.rdata")
pdginfo$nodes.code = pdginfo$node.df$desig
names(pdginfo$nodes.code) = pdginfo$node.df$node

##
## get output log of PDG 2015 fit with no error scaling
##
pdg15 = list()
pdg15$fit.txt = readLines("s035-fit-ns-all-2016-01-28.fit")
pdg15 = get.pdg.summary(pdg15)
pdg15 = get.pdg.input.params(pdg15)
pdg15 = get.pdg.input.nodes(pdg15)
##--- set nodes of parameters once the node - parameter relation is known
pdg15$params.in.node = pdg15$quant.node[pdg15$params.in.descr]
names(pdg15$params.in.node) = pdg15$params.in.plab
pdg15 = get.pdg.fitted.params(pdg15)
pdg15 = get.pdg.fitted.nodes(pdg15)
pdg15 = get.pdg.input.meas(pdg15)

##
## get data of HFAG fit with same inputs as PDG 2015
##
hfag.data.fname = "pdgfit-pdg-meas.rdata"
ht = load.in.list(hfag.data.fname)
cat(paste0("HFAG reproducing PDG fit file '", hfag.data.fname, "' read\n"))
##--- get un-modified HFAG input cards data
ht = get.hfag.input(ht)

##
## get data of HFAG 2014 reference fit
##
## some measurements tags are old:
##
## > grep("2012", names(htref$measurements), value=TRUE)
## [1] "BaBar.Gamma47.pub.LEES.2012Y"  "BaBar.Gamma50.pub.LEES.2012Y"
## [3] "BaBar.Gamma807.pub.LEES.2012Y" "BaBar.Gamma808.pub.LEES.2012Y"
## [5] "BaBar.Gamma811.pub.LEES.2012X" "BaBar.Gamma812.pub.LEES.2012X"
## [7] "BaBar.Gamma821.pub.LEES.2012X" "BaBar.Gamma822.pub.LEES.2012X"
## [9] "BaBar.Gamma831.pub.LEES.2012X" "BaBar.Gamma832.pub.LEES.2012X"
##[11] "BaBar.Gamma833.pub.LEES.2012X" "BaBar.Gamma910.pub.LEES.2012X"
##[13] "BaBar.Gamma911.pub.LEES.2012X" "BaBar.Gamma920.pub.LEES.2012X"
##[15] "BaBar.Gamma930.pub.LEES.2012X" "BaBar.Gamma941.pub.LEES.2012X"
##[17] "BaBar.Gamma942.pub.LEES.2012X" "BaBar.Gamma943.pub.LEES.2012X"
##[19] "BaBar.Gamma944.pub.LEES.2012X" "BaBar.Gamma950.pub.LEES.2012X"
##[21] "BaBar.Gamma951.pub.LEES.2012X" "BaBar.Gamma952.pub.LEES.2012X"
##> grep("2014vpc", names(htref$measurements), value=TRUE)
##[1] "Belle.Gamma33.pub.Ryu:2014vpc" "Belle.Gamma35.pub.Ryu:2014vpc"
##[3] "Belle.Gamma37.pub.Ryu:2014vpc" "Belle.Gamma40.pub.Ryu:2014vpc"
##[5] "Belle.Gamma42.pub.Ryu:2014vpc" "Belle.Gamma47.pub.Ryu:2014vpc"
##[7] "Belle.Gamma50.pub.Ryu:2014vpc"
##
hfag.data.ref.fname = "../TauFit/average2-aleph-hcorr_13-ref.rdata"
htref = load.in.list(hfag.data.ref.fname)
cat(paste0("HFAG fit reference file '", hfag.data.ref.fname, "' read\n"))
##--- get un-modified HFAG input cards data
htref = get.hfag.input(htref)

##
## get data of HFAG fit with HFAG inputs and PDG constraints with symbolic parameters, updated to PDG 2015
##
hfag.data.ref2.fname = "pdgfit-hfag-meas2.rdata"
htref2 = load.in.list(hfag.data.ref2.fname)
cat(paste0("HFAG fit reference file '", hfag.data.ref2.fname, "' read\n"))
##--- get un-modified HFAG input cards data
htref2 = get.hfag.input(htref2)

##
## get data of HFAG fit, no aleph-hcorr, no Belle incl K0S, with unitarity
##
hfag.data.ref3.fname = "pdgfit-step-0.rdata"
htref3 = load.in.list(hfag.data.ref3.fname)
cat(paste0("HFAG fit reference file '", hfag.data.ref3.fname, "' read\n"))
##--- get un-modified HFAG input cards data
htref3 = get.hfag.input(htref3)

##
## get data of HFAG PDG 2016 release candidate fit
##
hfag.data.rc.fname = "pdgfit-step-6.rdata"
htrc = load.in.list(hfag.data.rc.fname)
cat(paste0("HFAG-PDG 2016 file '", hfag.data.rc.fname, "' read\n"))
##--- get un-modified HFAG input cards data
htrc = get.hfag.input(htrc)

##
## if PDG nodes are not in HFAG print them and stop
##
nodes.to.define = which(! (pdg15$quant.node %in% ht$quant.node[ht$quant.node != ""]))
if (length(nodes.to.define) > 0) {
  cat("please define:\n")
  print(cbind(pdg15$quant.node[nodes.to.define]))
  stop()
}

##
## print PDG fit parameters definitions
##
cat("##\n")
cat("## PDG 2015 fit parameters definitions\n")
cat("##\n")

cat(paste(
  format(pdg15$params.in.node, width=10),
  format(pdg15$params.in.plab, width=10),
  pdg15$quant.descr[pdg15$params.in.node], collapse="\n"
  ))
cat("\n")

##
## print PDG nodes definitions (nodes that are not fit parameters)
##
cat("##\n")
cat("## PDG 2015 definitions of nodes that are not fit parameters\n")
cat("##\n")

cat(paste(
  format(ht$quant.node[hfag.gamma.pdg.derived(ht, pdg15)], width=10),
  format(pdg15$quant.plab[ht$quant.node[hfag.gamma.pdg.derived(ht, pdg15)]], width=10),
  pdg15$quant.descr[ht$quant.node[hfag.gamma.pdg.derived(ht, pdg15)]], collapse="\n"
  ))
cat("\n")

##
## print PDG 2015 vs. HFAG-like-PDG-2015 input measurements mismatches
##
cmp.meas.pdg15.ht = cmp.meas.pdg.hfag(pdg15, ht)
minimum.mism = 1e-6
cat("##\n")
cat("## PDG 2015 vs. HFAG 2014 input measurements mismatches >"); cat(minimum.mism); cat(", in decreasing size\n")
cat("## (mismatches are quadrature sum of the relative discrepancies in value and total uncertainty)\n")
cat("##\n")
pr.mism.meas.pdg.hfag(cmp.meas.pdg15.ht, minimum.mism)

##
## print PDG 2015 vs. HFAG-PDG 2016 input measurements mismatches
##
cmp.meas.pdg15.htrc = cmp.meas.pdg.hfag(pdg15, htrc)
minimum.mism = 1e-6
cat("##\n")
cat("## PDG 2015 vs. HFAG-PDG 2016 input measurements mismatches >"); cat(minimum.mism); cat(", in decreasing size\n")
cat("## (mismatches are quadrature sum of the relative discrepancies in value and total uncertainty)\n")
cat("##\n")
pr.mism.meas.pdg.hfag(cmp.meas.pdg15.htrc, minimum.mism)

##
## print PDG vs. HFAG node descr mismatches
##
cat("##\n")
cat("## PDG 2015 vs. HFAG 2014 mismatches in node descriptions\n")
cat("##\n")
pr.mism.descr.pdg.hfag(pdg15, ht)

##
## build list of MODMEAS KEEP cards to define measurements to be used
##
pdg.meas.tags = sapply(ht$measurements[cmp.meas.pdg15.ht$pdg.meas.name],
  function(meas) {paste(meas$tags, collapse=" ")})
cat("##\n")
cat("## HFAG 2014 measurements matching PDG 2015 input measurements\n")
cat("##\n")
cat(paste("MODMEAS KEEP", pdg.meas.tags, collapse="\n"))
cat("\n")
rm(pdg.meas.tags)

##
## COMBINE cards to replicate PDG 2015 with HFAG
##
cat("##\n")
cat("## HFAG 2014 quantities corresponding to PDG 2015 fit nodes\n")
cat("##\n")
pr.hfag.gamma.pdg.nodes(pdg15, ht)

##
## get symbolic expressions for PDG 2015 constraint coefficients
##
pdg15 = get.pdg.constr.params.symbolic(pdg15, ht)
pr.pdg.constr.params.symbolic(pdg15)

cat("##\n")
cat("## PDG constraints with PDG nodes and parameters\n")
cat("##\n")
cat(paste(pdg.constr.to.text(pdg15$constr), collapse="\n"))
cat("\n")

cat("##\n")
cat("## PDG constraints in human form\n")
cat("##\n")
cat(paste(pdg.constr.to.text.human(pdg15, ht), collapse="\n\n"))
cat("\n")

cat("##\n")
cat("## PDG constraints in HFAG gamma notation\n")
cat("##\n")
cat(paste(pdg.constr.to.text.hfag(pdg15, ht), collapse="\n"))
cat("\n")

cat("##\n")
cat("## non-identity PDG constraints in HFAG cards notation\n")
cat("##\n")
rc = cat(paste(pdg.constr.convert.nontrivial(pdg.constr.to.text.hfag(pdg15, ht)), collapse="\n"))
cat("\n")
cat("NLCONSTRAINT Unitarity.c 1 \"", paste(ht$quant.gamma[pdg15$params.in.node], collapse=" + "), "\"\n", sep="")

##
## find out which measurements have systematic errors
## HFAG 2014: 45 measurements of 191 have listed systematic terms
##
if (FALSE) {
  meas.with.syst = sapply(
    htref2$measurements[names(htref2$measurements) %in% cmp.meas.pdg15.ht$pdg.meas.name],
    function(meas) {!is.null(meas$syst.terms)})

  ## names(meas.with.syst) = gsub("Ryu:2014vpc", "RYU.14vpc", names(meas.with.syst), fixed=TRUE)
  ## names(meas.with.syst) = gsub("LEES.2012", "LEES.12", names(meas.with.syst), fixed=TRUE)

  cat("##\n")
  cat("## PDG measurements that have systematics in HFAG\n")
  cat("##\n")

  cat("  ", paste(names(meas.with.syst), collapse="\n  "), "\n", sep="")
}

##
## comparisons between PDG 2015 and HFAG-PDG 2016
##

##
## list measurement usage differences between PDG 2015 and HFAG-PDG 2016
##
cmp.meas.pdg15.htrc = cmp.meas.pdg.hfag(pdg15, htrc)
pr.meas.use.mism.pdg.hfag(cmp.meas.pdg15.htrc$pdg.meas.name, htrc, "PDG 2015", "HFAG-PDG 2016")

##
## compare constraints for PDG 2015 vs. HFAG-PDG 2016
##
cmp.pgd.hfag.constr(pdg15, htrc, "PDG 2015", "HFAG-PDG 2016")

##
## quantities fitted in PDG 2015 but not fitted in HFAG-PDG 2016
##
cat("##\n")
cat("## quantities fitted in PDG 2015 but not fitted in HFAG-PDG 2016\n")
cat("##\n")

in.pdg.not.in.hfag = !(pdg15$nodes.fit$node %in% htrc$quant.node[names(htrc$quant.val)])
rc = lapply(split(pdg15$nodes.fit[in.pdg.not.in.hfag, ], 1:nrow(pdg15$nodes.fit[in.pdg.not.in.hfag, ])),
  function(cmp.fit) {
    cat(sprintf("%-12s %-12s %s\n", cmp.fit$node, htrc$quant.gamma[cmp.fit$node], cmp.fit$descr))
  })

##
## quantities fitted in HFAG-PDG 2016 but not fitted in PDG 2015
##
cat("##\n")
cat("## quantities fitted in HFAG-PDG 2016 but not fitted in PDG 2015\n")
cat("##\n")

in.hfag.not.in.pdg = setdiff(names(htrc$quant.val), htrc$quant.gamma[pdg15$nodes.fit$node])
rc = lapply(in.hfag.not.in.pdg,
  function(gamma) {
    cat(sprintf("%-12.12s %-12s %s\n", htrc$quant.node[gamma], gamma, htrc$quant.descr[gamma]))
  })

##
## compare fit results for PDG 2015 vs. HFAG-PDG 2016
##
cmp.pdg.hfag.fit.results(pdg15, ht, "PDG 2015", "HFAG-PDG 2016")

##
## base quantities
##
quant.base.names = alucomb2.base.quant(htrc$combination)

cat("##\n")
cat("## base quantities of HFAG-PDG 2016 fit\n")
cat("## number of base quantities:  =", length(quant.base.names)); cat("\n")
cat("## measurements - fit n.d.o.f. =", htrc$meas.num - htrc$dof, "\n")
cat("## measurements                =", htrc$meas.num, "\n")
cat("## n.d.o.f.                    =", htrc$dof, "\n")
cat("##\n")
## alucomb2.clean.print(quant.base.names)
rc = lapply(quant.base.names,
  function(gamma) {
    rc = cat(sprintf("%-12.12s %-12.12s %3d %s\n",
      gamma,
      htrc$quant.node[gamma],
      pdginfo$nodes.code[htrc$quant.node[gamma]],
      htrc$quant.descr[gamma]))
  })

##
## unitarity constraint quantities
##
quant.uni.names = alucomb2.unitarity.quant(htrc$combination)

cat("##\n")
cat("## quantities in unitarity constraint (", length(quant.uni.names), ")\n")
cat("##\n")
rc = lapply(quant.uni.names,
  function(gamma) {
    rc = cat(sprintf("%-12.12s %-12.12s %3d %#7.4f +-%#7.4f %s\n",
      gamma,
      htrc$quant.node[gamma],
      pdginfo$nodes.code[htrc$quant.node[gamma]],
      100*htrc$quant.val[gamma], 100*htrc$quant.err[gamma],
      htrc$quant.descr[gamma]))
  })

##
## unitarity constraint quantities, non base quantities
##
quant.uni.names = alucomb2.unitarity.quant(htrc$combination)

cat("##\n")
cat("## non basis quantities in unitarity constraint\n")
cat("##\n")
rc = lapply(setdiff(quant.uni.names, quant.base.names),
  function(gamma) {
    rc = cat(sprintf("%-12.12s %-12.12s %3d %#7.4f +-%#7.4f %s\n",
      gamma,
      htrc$quant.node[gamma],
      pdginfo$nodes.code[htrc$quant.node[gamma]],
      100*htrc$quant.val[gamma], 100*htrc$quant.err[gamma],
      htrc$quant.descr[gamma]))
    index = grep(paste0(gamma,".c"), names(htrc$combination$constr.all.val))
    ##print(gamma)
    ##print(index)
    ##print(str_match(htrc$combination$constr.all.str.expr[[index]], "^-\\S+ [+] \\((.*)\\)$"))
    gamma.in.base.quant = str_match(htrc$combination$constr.all.str.expr.input[[index]], "^-\\S+ [+] \\((.*)\\)$")[,2]
    cat(" ", "=", gamma.in.base.quant); cat("\n")
  })

##
## PDG 2015 correlations for recent data
##

## cat("##\n")
## cat("## PDG 2015 correlations for recent data\n")
## cat("##\n")

papers = c(
  "AUBERT.07AP",
  "AUBERT.08",
  "AUBERT.10F",
  "LEES.12X",
  "LEES.12Y"
  )

rc = lapply(papers,
  function(paper) {
    grep(paper,paste(pdg15$meas.in$author, pdg15$meas.in$year, sep="."))
  })
meas.in.sel = do.call(c, rc)

rc = lapply(pdg15$meas.in[meas.in.sel, ],
  function(meas.in) {
  })

##
## save
##
cat("##\n")
cat("## save HFAG cards\n")
cat("##\n")
##--- save pdg measurements
save.pdg.measurements(pdg15, ht, cmp.meas.pdg15.ht$pdg.match, file.dir="pdg")
##--- save HFAG measurements
save.hfag.measurements(ht$measurements, ht, "hfag")
