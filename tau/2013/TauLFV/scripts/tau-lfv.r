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

##
## order list containing gamma by type and then gamma
##
tau.lfv.br.order = function(list, info) {
  exps = sapply(list, function(el) el$exp)
  gammas = sapply(list, function(el) el$gamma)
  type.nums = sapply(info[as.character(gammas)], function(el) el$type.num)
  list[order(type.nums, gammas, exps)]
}

## /////////////////////////////////////////////////////////////////////////////
##
##	code
##
## /////////////////////////////////////////////////////////////////////////////

tau.lfv.data = yaml.load_file("tau-lfv-data.yaml")
ofname = "tau-lfv-data.tex"

##
## get all types (categories) of tau LFV modes
## each category gets for sorting purposes the smallest
## BR gamma number that belongs to the category
##
gammas = sapply(tau.lfv.data$gamma, function(el) el$gamma)
types = sapply(tau.lfv.data$gamma, function(el) el$type)
types.uniq = unique(types)
types.uniq.num = unname(sapply(types.uniq, function(type) min(gammas[types == type])))
tau.lfv.data$gamma = lapply(tau.lfv.data$gamma, function(el) {
  el$type.num = types.uniq.num[types.uniq == el$type]
  el
})

##--- get gamma info in data.frame
gamma.df = list.to.df(tau.lfv.data$gamma)
names(tau.lfv.data$gamma) = gamma.df$gamma

##
## remove elements according to references
##
## remove preliminary references
## - Hayasaka:2011zz Belle 2011 ell pi0, ell eta, ell eta'
## - Hayasaka:2012pj Belle 2012 pi/K lambda(bar)
## - Lafferty:2007zz BaBar 2007 pi/K lambda(bar)
##
## remove limits from CLEO
##

refs.prelim = c(
  "Hayasaka:2011zz",
  "Hayasaka:2012pj",
  "Lafferty:2007zz"
)

##--- order limits by gamma
tau.lfv.data$limits = tau.lfv.br.order(tau.lfv.data$limits, tau.lfv.data$gamma)

##--- get list of prelim limits, select them and remove them
rc = sapply(tau.lfv.data$limits, function(br) {
  ifelse(br$ref %in% refs.prelim, FALSE, TRUE)
})
tau.lfv.data$limits.prelim = tau.lfv.data$limits[!rc]
tau.lfv.data$limits = tau.lfv.data$limits[rc]

##--- get list of CLEO limits and remove
rc = sapply(tau.lfv.data$limits, function(br) {
  ifelse(br$exp == "CLEO", FALSE, rc)
})
tau.lfv.data$limits = tau.lfv.data$limits[rc]


##--- order combs by gamma
tau.lfv.data$combs = tau.lfv.br.order(tau.lfv.data$combs, tau.lfv.data$gamma)

##--- remove combs using prelim limits
rc = sapply(tau.lfv.data$combs, function(br) {
  refs = as.vector(strsplit(br$refs, ",", fixed=TRUE)[[1]])
  ifelse(any(refs %in% refs.prelim), FALSE, TRUE)
})
tau.lfv.data$combs = tau.lfv.data$combs[rc]

##--- order extra info for combinations by gamma
tau.lfv.data$combs.extra = tau.lfv.br.order(tau.lfv.data$combs.extra, tau.lfv.data$gamma)

##--- remove extra info for combinations of preliminary limits
rc = sapply(tau.lfv.data$combs.extra, function(br) {
  ifelse(br$ref %in% refs.prelim, FALSE, TRUE)
})
tau.lfv.data$combs.extra = tau.lfv.data$combs.extra[rc]

##
## prepare LaTeX source
##

out.txt = character(0)

rc = mapply(function(gamma, descr) {
  key = paste0("g", gamma)
  alurep.tex.cmd.short(paste0(key, ".descr"), paste0("\\ensuremath{", descr, "}"))
}, gamma.df$gamma, gamma.df$descr)
out.txt = c(out.txt, rc)

##
## extra info for combinations
##

rc = lapply(tau.lfv.data$combs.extra, function(br) {
  gamma.info = tau.lfv.data$gamma[[as.character(br$gamma)]]
  rc = paste0(
    "\\htCombExtraLine%\n",
    "  {\\ensuremath{\\Gamma_{", br$gamma, "} = ", gamma.info$descr, "}}%\n",
    "  {", br$exp, "}%\n",
    "  {", br$ref, "}%\n",
    ## "  {", br$lumi, "}%\n",
    ## "  {", br$cross.section, "}%\n",
    "  {", sprintf("%.0f", br$num.tau), "}%\n",    
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

##
## limits
##

descr.last = ""
type.last = ""
sep.first = TRUE
rc = lapply(tau.lfv.data$limits, function(br) {
  gamma.info = tau.lfv.data$gamma[[as.character(br$gamma)]]
  rc = character(0)
  descr = ifelse(gamma.info$descr != descr.last, gamma.info$descr, "")
  descr.last <<- gamma.info$descr
  type = ifelse(gamma.info$type != type.last, gamma.info$type, "")
  type.last <<- gamma.info$type

  if (!sep.first && type != "") {
    rc = c(rc, "\\midrule")
  } else {
    sep.first <<- FALSE
  }
  
  rc = c(
    rc,
    paste0(
      "\\htLimitLine%\n",
      "  {\\ensuremath{\\Gamma_{", br$gamma, "} = ", gamma.info$descr, "}}%\n",
      "  {\\ensuremath{", type, "}}%\n",
      "  {\\ensuremath{", sub("^(.*)e([+-]*)0*(\\d*)$", "\\1 \\\\cdot 10^{\\2\\3}", sprintf("%.1e", br$limit)), "}}%\n",
      "  {", br$exp, "}%\n",
      "  {", br$ref, "}"
    ))
})

rc = alurep.tex.cmd(
  "LimitLines",
  paste0(unlist(rc), collapse="%\n")
)
out.txt = c(out.txt, rc)

##
## preliminary limits
##

descr.last = ""
type.last = ""
sep.first = TRUE
rc = lapply(tau.lfv.data$limits.prelim, function(br) {
  gamma.info = tau.lfv.data$gamma[[as.character(br$gamma)]]
  rc = character(0)
  descr = ifelse(gamma.info$descr != descr.last, gamma.info$descr, "")
  descr.last <<- gamma.info$descr
  type = ifelse(gamma.info$type != type.last, gamma.info$type, "")
  type.last <<- gamma.info$type

  if (!sep.first && type != "") {
    rc = c(rc, "\\midrule")
  } else {
    sep.first <<- FALSE
  }
  
  rc = c(
    rc,
    paste0(
      "\\htLimitLine%\n",
      "  {\\ensuremath{\\Gamma_{", br$gamma, "} = ", gamma.info$descr, "}}%\n",
      "  {\\ensuremath{", type, "}}%\n",
      "  {\\ensuremath{", sub("^(.*)e([+-]*)0*(\\d*)$", "\\1 \\\\cdot 10^{\\2\\3}", sprintf("%.1e", br$limit)), "}}%\n",
      "  {", br$exp, "}%\n",
      "  {", br$ref, "}"
    ))
})

rc = alurep.tex.cmd(
  "PrelimLimitLines",
  paste0(unlist(rc), collapse="%\n")
)
out.txt = c(out.txt, rc)

##
## combinations
##

descr.last = ""
type.last = ""
sep.first = TRUE
rc = lapply(tau.lfv.data$combs, function(br) {
  gamma.info = tau.lfv.data$gamma[[as.character(br$gamma)]]
  rc = character(0)
  descr = ifelse(gamma.info$descr != descr.last, gamma.info$descr, "")
  descr.last <<- gamma.info$descr
  type = ifelse(gamma.info$type != type.last, gamma.info$type, "")
  type.last <<- gamma.info$type

  if (!sep.first && type != "") {
    rc = c(rc, "\\midrule")
  } else {
    sep.first <<- FALSE
  }
  
  rc = c(
    rc,
    paste0(
      "\\htCombLimitLine%\n",
      "  {\\ensuremath{\\Gamma_{", br$gamma, "} = ", gamma.info$descr, "}}%\n",
      "  {\\ensuremath{", type, "}}%\n",
      "  {\\ensuremath{", sub("^(.*)e([+-]*)0*(\\d*)$", "\\1 \\\\cdot 10^{\\2\\3}", sprintf("%.1e", br$limit)), "}}%\n",
      "  {\\cite{", br$refs, "}}"
    ))
})

rc = alurep.tex.cmd(
  "CombLines",
  paste0(unlist(rc), collapse="%\n")
)
out.txt = c(out.txt, rc)

cat(paste(out.txt, collapse="\n"), "\n", file=ofname, sep="")
cat("file '", ofname, "' created\n", sep="")
