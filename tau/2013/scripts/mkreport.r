#!/usr/bin/env Rscript

##
## mkreport.r
##

require("optparse", quietly=TRUE)
require("stringr", quietly=TRUE)
source("../../../Common/bin/alucomb2-utils.r")
source("../../../Common/bin/alucomb2-fit.r")
source("../../../Common/bin/alureport.r")
source("../../../Common/bin/aluelab3.r")

## ////////////////////////////////////////////////////////////////////////////
## functions

##
## example function to use vars of an env in the function env
##
mkreport.rc = function(alucomb.env) {
  attach(alucomb.env)
  rc = try({
  })
  detach(alucomb.env)
  return(invisible(NULL))
}

##
## for the spec. measurement return list with:
## - experiment
## - latex \cite{} reference
## - type: prelim or pub
## - formatted val +- stat +- syst
##
get.tex.meas = function(meas, precision, order) {
  meas.item = list()
  meas.item$exp = meas$tags[1]
  meas.item$exp = sub("BaBar", "\\\\babar", meas.item$exp, ignore.case=TRUE)
  meas.item$ref = paste("\\cite{", get.reference(meas$tags), "}", sep="")
  meas.item$type = meas$tags[3]
  meas.item$value = alurep.tex.meas.val(meas, precision, order)
  return(meas.item)
}

##
## latex defs for all measurements
##
get.tex.meas.defs = function() {
  meas.used = measurements[combination$measurements]
  meas.used.names = names(meas.used)
  meas.used.names = meas.used.names[order(meas.used.names)]

  rc = mapply(function(meas.name, meas) {
    exp = meas$tags[1]
    exp = sub("BaBar", "\\\\babar", exp)
    ref = get.reference(meas$tags)
    type = meas$tags[3]
    quant = meas$quant
    rc2 = alurep.tex.meas.val.card.fields(meas)
    paste("\\htmeasdef{", meas.name, "}{",
          quant, "}{", exp, "}{", ref, "}{",
          rc2$quant, "}{",
          rc2$val, "}{", rc2$stat, "}{", rc2$syst, "}%",
          sep="")
  }, meas.used.names, measurements[meas.used.names])
  return(paste(rc, collapse="\n"))
}

##
## return TeX defs for listing a quantity and its measurements
## with the same precision and order
##
## - create TeX defs for the quantiy and all its measurements
## - create a
## return body of latex tabular environment for the requested quantities
## include
## - quantity description, HFAG average
## - list of experimental measurements, with reference
##
get.tex.table = function(quant.names, perc=FALSE) {
  quant.order = order(alucomb2.gamma.num.id(quant.names))
  quant.names = quant.names[quant.order]

  ##--- get TeX defs for quantities and their measurements
  qm.defs = mapply(function(quant.name, quant) {
    ##--- quantity value according to precision / order found for quantity and its measurements
    rc = alurep.precision.order.quant(quant.name, perc=perc, with.meas=TRUE)
    precision = rc$precision
    order = rc$order

    ##--- TeX def for the quantity
    rc = paste(
      "\\htdef{", quant.name, ".qt}{\\ensuremath{",
      alurep.tex.val.err.prec.ord(quant.val[quant.name], quant.err[quant.name], precision, order, perc=FALSE),
      "}}%", sep=""
      )
    
    ##--- TeX defs for the quantity measurements
    meas.names = alurep.meas.quant(quant.name, delta)
    meas.lines = sapply(meas.names, function(meas.name) {
      meas.item = get.tex.meas(measurements[[meas.name]], precision, order)
      rc = paste(
        "\\htdef{", meas.name, ",qt}{", meas.item$value, "}%",
        sep=""
        )
    })
    
    ##--- join all defs into one string
    if (length(meas.lines) > 0) {
      meas.lines = paste(meas.lines, collapse="\n")
      rc = paste(rc, meas.lines, sep=" \n")
    }

    return(rc)
  },
    quant.names, combination$quantities[quant.names])
  qm.defs = paste(qm.defs, collapse=" \n")

  ##--- get TeX defs for a tabular body listing a quantity and its measurements (using previous defs)
  qm.defs2 = mapply(function(quant.name, quant) {
    ##--- TeX expr for quantity
    quant.descr = paste("\\htuse{", quant.name, ".gn}", " = ", "\\htuse{", quant.name, ".td}", sep="")
    quant.descr = paste("\\begin{ensuredisplaymath}\n", quant.descr, "\n\\end{ensuredisplaymath}\n", sep="")

    ##--- get quantity value with proper precision and order
    rc = paste(quant.descr, " & \\htuse{", quant.name, ".qt} & \\hfagFitLabel", sep="")

    ##--- get measurements value, experiment and citation
    meas.names = alurep.meas.quant(quant.name, delta)
    meas.lines = sapply(meas.names, function(meas.name) {
      rc = paste(
        paste("\\htuse{", meas.name, ",qt}", sep=""),
        paste("\\htuse{", meas.name, ",exp}", sep=""),
        paste("\\htuse{", meas.name, ",ref}", sep=""),
        sep=" & ")
      return(rc)
    })
    meas.lines = paste(meas.lines, collapse=" \\\\\n")
    
    ##--- define macro for tabular body
    quant.meas.tex = paste("\\htdef{", quant.name, ".qm}{%\n", rc, sep="")
    ##--- add measurements if there are any
    if (nchar(meas.lines) > 0) {
      quant.meas.tex = paste(quant.meas.tex, "\\\\\n", meas.lines, "\n", sep="")
    }
    quant.meas.tex = paste(quant.meas.tex, "}%", sep="")
    return(quant.meas.tex)
  },
    quant.names, combination$quantities[quant.names])
  qm.defs2 = paste(qm.defs2, collapse=" \n")

  ##--- join the two previous blocks of definitions
  qm.defs = paste(qm.defs, qm.defs2, sep="\n")
  
  ##--- build complete tabular body of all quantities and measurements
  qm.table = mapply(function(quant.name, quant) {
    rc = paste("\\htuse{", quant.name, ".qm}", sep="")
    return(rc)
  },
    quant.names, combination$quantities[quant.names])
  
  with.meas = TRUE
  if (with.meas) {
    qm.table = paste(qm.table, collapse=" \\\\\n\\midrule\n")
  } else {
    qm.table = paste(qm.table, collapse=" \\\\\n")
  }

  return(list(defs=qm.defs, table=qm.table))
}

##
## return body of latex tabular environment for the requested quantities
## include
## - quantity names
## - constraint coefficients
## must define two macros to expand the two fields
## \newcommand{\htQuantLabel}[1]{\htuse{#1.td}}
## \newcommand{\htQuantValue}[2]{#1}
##
get.tex.table.constr = function(comb.constr, precision, order, constr.prec, constr.order) {
  ##--- order by gamma number
  comb.constr = comb.constr[order(alucomb2.gamma.num.id(names(comb.constr)))]
  rc = mapply(
    function(quant.name, coeff.val) {
      rc = paste(
        "\\htConstrLine{", quant.name, "}{",
        alurep.tex.val.err.prec.ord.noee(quant.val[quant.name], quant.err[quant.name], precision, order), "}{",
        alurep.tex.val.prec.ord(coeff.val, constr.prec, constr.order, ), "}{",
        sprintf("%d", order), "}{",
        sprintf("%d", constr.order), "}",
        sep="")
    }, names(comb.constr), comb.constr)
  return(paste(rc, collapse=" \n"))
}

##
## return body of latex tabular environment for the requested quantities
## include
## - quantity description and HFAG average
## must define two macros to expand the two fields
## \newcommand{\htQuantLabel}[1]{\htuse{#1.td}}
## \newcommand{\htQuantValue}[2]{#1}
## +++ unused quant argument
## +++ uses quant.val from outside
##
get.tex.table.simple = function(quant.names, precision, order) {
  quant.order = order(alucomb2.gamma.num.id(quant.names))
  quant.names = quant.names[quant.order]
  rc = mapply(
    function(quant.name, quant) {
      rc = paste(
        "\\htQuantLine{", quant.name, "}{",
        alurep.tex.val.err.prec.ord.noee(quant.val[quant.name], quant.err[quant.name], precision, order), "}{",
        sprintf("%d", order), "}",
        sep="")
    }, quant.names, combination$quantities[quant.names])
  return(paste(rc, collapse=" \n"))
}

##
## return measurements by collaboration
##
get.tex.meas.by.collab = function() {
  toTex = TrStr.num2tex$new()
  collab.meas = sapply(measurements[combination$measurements], function(meas) meas$tags[1])
  collabs = sort(unique(collab.meas))
  collab.nmeas = sapply(collabs, function(collab) sum(collab.meas == collab))
  return(alurep.tex.cmd.short(paste("NumMeas", collabs, sep=""), collab.nmeas))
}

##
## return table body with measurements by reference
##
get.tex.meas.by.ref = function() {
  meas.paper = list()
  rc = mapply(function(meas.name, meas) {
    quant.name = meas$quant

    pub.flag = (regexpr("^pub", meas$tags[3], ignore.case=TRUE) != -1)
    if (pub.flag) {
      pub = "pub"
      cit = paste(meas$tags[-(1:3)], collapse=" ")
    } else {
      pub = ""
      cit = paste(meas$tags[1], "prelim.", meas$tags[-(1:3)], sep=" ")
    }
    expnt = meas$tags[1]
    ref = get.reference(meas$tags)
    
    return(list(cit=cit, pub.flag=pub, expnt=expnt, ref=ref, quant=quant.name, meas=meas.name))
  }, combination$measurements, measurements[combination$measurements])

  ##--- convert to data.frame
  rc = do.call(rbind.data.frame, apply(rc, 2, function(x) {as.data.frame(t(unlist(x)))}))

  ##--- get citations grouped by experiment
  expnt = unique(rc$expnt)
  ##--- sort experiments alphabetically
  expnt = expnt[order(expnt)]
  cits = lapply(expnt, function(x) {
    df = subset(rc, expnt == x & pub.flag == "pub")
    cits = as.character(unique(df$cit))
    cits = cits[order(cits)]
    df = subset(rc, expnt == x & pub.flag != "pub")
    cits.prelim = as.character(unique((df$cit)))
    cits.prelim = cits.prelim[order(cits.prelim)]
    return(c(cits, cits.prelim))
  })
  cits = unlist(cits, recursive = TRUE)
  
  meas.by.ref = lapply(cits, function(x) {
    df = subset(rc, cit == x)
    ##--- reorder measurements according to quantity
    quant.order = order(alucomb2.gamma.num.id(df$quant))
    df = df[quant.order,]

    expnt = df$expnt[1]
    cite.tex = paste("\\cite{", df$ref[1], "}", sep="")
    if (df$pub.flag[1] == "pub") {
      expnt = sub("^BaBar$", "\\\\babar", expnt, perl=TRUE, ignore.case=TRUE)
      ref.tex = paste(df$cit[1], " (", expnt, ") ", cite.tex, sep="")
    } else {
      cit = df$cit[1]
      cit = sub("^BaBar", "\\\\babar", cit, perl=TRUE, ignore.case=TRUE)
      ref.tex = paste(cit, cite.tex, sep=" ")
    }

    quant.descr = paste("\\htuse{", df$quant, ".gn}", " = ", "\\htuse{", df$quant, ".td}", sep="")
    quant.descr = paste("\\begin{ensuredisplaymath}\n", quant.descr, "\n\\end{ensuredisplaymath}", sep="")

    val.lines = paste(quant.descr, " & ", "\\htuse{", df$meas, "}", sep="")
    val.tex = paste(val.lines, collapse="\n\\\\\n")

    list(cit=x, ref.tex=ref.tex, val.tex=val.tex, cite.tex=cite.tex, collab.tex=expnt)
  })

  meas.by.ref.defs = lapply(meas.by.ref, function(x) {
    def.cite = paste("\\htdef{", x$cit, ".cite", "}{", x$cite.tex, "}%", sep="")
    def.collab = paste("\\htdef{", x$cit, ".collab", "}{", x$collab.tex, "}%", sep="")
    def.ref = paste("\\htdef{", x$cit, ".ref", "}{", x$ref.tex, "}%", sep="")
    def.meas = paste("\\htdef{", x$cit, ".meas", "}{%\n", x$val.tex, "}%", sep="")
    paste(def.cite, def.collab, def.ref, def.meas, sep="\n")
  })
  meas.by.ref.defs.tex = paste(meas.by.ref.defs, collapse="\n")

  meas.by.ref.tex = lapply(meas.by.ref, function(x) {
    rc = paste(
      "\\multicolumn{2}{l}{\\htuse{", x$cit, ".ref", "}} \\\\\n",
      "\\htuse{", x$cit, ".meas", "}",
      sep="")
  })
  meas.by.ref.tex = paste(meas.by.ref.tex, collapse=" \\\\\\hline\n")

  list(defs=meas.by.ref.defs.tex, table=meas.by.ref.tex)
}

##
## return latex code with the basis quantities correlation coefficients
##
get.tex.base.nodes.corr = function() {
  tex.all.tau.br.corr = NULL
  quant.names = alucomb2.base.quant(combination)

  ##--- tex code preceding the correlation table content
  corr.pre = c(
    "%%",
    "%% basis quantities correlation, @@num@@",
    "%%",
    "\\ifhevea\\begin{table}\\fi%% otherwise cannot have normalsize caption",
    "\\begin{center}",
    "\\ifhevea",
    "\\caption{Basis quantities correlation coefficients in percent, section @@num@@.\\label{tab:tau:br-fit-corr@@num@@}}%",
    "\\else",
    "\\begin{minipage}{\\linewidth}",
    "\\begin{center}",
    "\\captionof{table}{Basis quantities correlation coefficients in percent, section @@num@@.}\\label{tab:tau:br-fit-corr@@num@@}%",
    "\\fi",
    "\\begin{envsmall}",
    "\\begin{center}",
    "\\renewcommand*{\\arraystretch}{1.1}%",
    "\\begin{tabular}{@@tabcols@@}",
    "\\hline")

  ##--- tex code following the correlation table content
  corr.post = c(
    "\\\\\\hline",
    "\\end{tabular}",
    "\\end{center}",
    "\\end{envsmall}",
    "\\ifhevea\\else",
    "\\end{center}",
    "\\end{minipage}",
    "\\fi",
    "\\end{center}",
    "\\ifhevea\\end{table}\\fi")

  ##--- correlation of basis quantities in percent
  quant.corr.base = round(quant.corr[quant.names, quant.names] * 100)
  ##--- fix negative zero output by sprintf
  quant.corr.base = ifelse(quant.corr.base == 0, abs(quant.corr.base), quant.corr.base)

  ##--- write data in tables with at most spec. rows and columns
  coeff.per.row = 14
  coeff.per.col = 14
  inum = 1
  for (j in seq(1, length(quant.names), by=coeff.per.col)) {
    for (i in seq(1, length(quant.names), by=coeff.per.row)) {
      submat.txt = NULL
      for(ii in i:min(i+coeff.per.row-1, length(quant.names))) {
        if (ii<=j) next
        maxcol = min(ii-1,j+coeff.per.col-1)
        row.txt = paste("\\(", alurep.gamma.texlabel(quant.names[ii]), "\\)")
        row.txt = c(row.txt, sprintf("%4.0f", quant.corr.base[ii, j:maxcol]))
        row.txt = c(row.txt, rep("", length.out=min(j+coeff.per.col-1, length(quant.names))-j+1-(maxcol-j+1)))
        submat.txt = c(submat.txt, paste(row.txt, collapse=" & "))
      }
      if (is.null(submat.txt)) next
      label.txt = alurep.gamma.texlabel(quant.names[j:min(j+coeff.per.col-1, length(quant.names))])
      label.num = min(j+coeff.per.col-1, length(quant.names)) - j + 1
      label.txt = paste("\\(", label.txt, "\\)", sep=" ", collapse=" & ")
      label.txt = paste("", label.txt, sep=" & ")
      submat.txt = paste(submat.txt, collapse=" \\\\\n")
      submat.txt = paste(submat.txt, label.txt, sep=" \\\\\n")
      corr.pre.mod = gsub("@@num@@", as.character(inum), corr.pre)
      corr.pre.mod = gsub("@@tabcols@@", paste(rep("r", length.out=label.num+1), collapse=""), corr.pre.mod)
      ## corr.pre.mod = paste(corr.pre.mod, collapse="\n")
      tex.all.tau.br.corr = c(tex.all.tau.br.corr, corr.pre.mod, submat.txt, corr.post)
      inum = inum + 1
    }
  }
  return(paste(tex.all.tau.br.corr, collapse="\n"))
}

##--- get all constraints names that were used in the fit
get.constraints.used.names = function() {
  constr.used = combination$constr.all.lin | combination$constr.all.nl
  constr.used.names = names(constr.used[constr.used])
  constr.order = order(alucomb2.gamma.num.id(constr.used.names))
  constr.used.names = constr.used.names[constr.order]
  return(constr.used.names)
}

##--- get all used constraint equations
get.tex.constraints.used = function() {
  constr.used.names = get.constraints.used.names()
  comb.str = combination$constr.all.str.expr[constr.used.names]
  comb.val = combination$constr.all.val[constr.used.names]
  comb.nl = intersect(constr.used.names, names(combination$constr.nl.str.expr))
  if (length(comb.nl) > 0) {
    comb.str[comb.nl] = combination$constr.nl.str.expr[comb.nl]
    comb.val[comb.nl] = combination$constr.nl.str.val[comb.nl]
  }
  return(alurep.tex.constraint(unlist(comb.val), unlist(comb.str)))
}

##
## return latex code with constraint equations
## (only constraints not corresponding to ratio of BRs)
##
get.tex.constraint.defs = function() {
  rc = get.tex.constraints.used()
  return(paste("\\htconstrdef{", names(rc$left), "}{", rc$left, "}{", rc$right, "}{", rc$right.split, "}%", sep=""))
}

##
## return latex code with constraint equations
## (only constraints not corresponding to ratio of BRs)
##
get.tex.constraint.equations = function() {
  rc = get.tex.constraints.used()
  constr.used.names = names(rc$left)

  constr.used.names = get.constraints.used.names()
  
  ##
  ## selection of constraints related to BRs that are ratios of 2 BRs
  ## in the report, these constraints are listed in the definition of the BRs not in the list of constraints
  ## also select the dummy constraint to compute the Unitarity discrepancy
  ##
  sel = grepl("Gamma\\d+by\\d+.*|Unitarity", constr.used.names, perl=TRUE)
  constr.names = constr.used.names[!sel]

  return(paste("\\begin{align*}\n",
               "\\htuse{", constr.names, ".left}",
               " ={}& ",
               "\\htuse{", constr.names, ".right.split}",
               "\n\\end{align*}", sep="", collapse="\n"))
}

##
## create .tex files for HFAG report
##
mkreport = function(fname) {
  load(fname, .GlobalEnv)

  ##--- compute from saved data
  assign("quant.err", sqrt(diag(quant.cov)), envir = .GlobalEnv)
  quant.err.nonzero = (quant.err != 0)
  quant.corr = quant.cov
  quant.corr[quant.err.nonzero, quant.err.nonzero] = cov2cor(quant.cov[quant.err.nonzero, quant.err.nonzero])
  assign("quant.corr", quant.corr, envir = .GlobalEnv)

  ##--- special rounding to cope with rounding errors
  eval(parse(text="quant.val[\"GammaAll\"]=round(quant.val[\"GammaAll\"], digits=12)"), .GlobalEnv)
  eval(parse(text="quant.err[\"GammaAll\"]=round(quant.err[\"GammaAll\"], digits=12)"), .GlobalEnv)
  eval(parse(text="quant.val[\"Gamma998\"]=round(quant.val[\"Gamma998\"], digits=12)"), .GlobalEnv)
  eval(parse(text="quant.err[\"Gamma998\"]=round(quant.err[\"Gamma998\"], digits=12)"), .GlobalEnv)

  ##--- assemble file name
  fname = basename(fname)
  fname = sub("[.][^.]*$", "", fname, perl=TRUE)
  fname = sub("average[^-]*-*", "", fname)
  if (fname != "") fname = paste0("-", fname)
  fname = file.path("../report", paste0("tau-br-fit", fname))

  fname = paste0(fname, ".tex")
  cat("", file=fname)
  cat("file '", fname, "' created\n", sep="")

  ##
  ## prepare inputs for TeX defs
  ##

  quant.names.nonratio = alucomb2.nonratio.quant(combination)
  quant.names.ratio = alucomb2.ratio.quant(combination)

  ##--- all fitted quantities but GammaAll and Gamma998, and which are not ratios of BRs
  quant.num.nonratio = length(quant.names.nonratio)
  ##--- all fitted quantities but GammaAll and Gamma998, and which are ratios of BRs
  quant.num.ratio = length(quant.names.ratio)
  ##--- all fitted quantities but GammaAll and Gamma998
  quant.num.nonunitarity = quant.num.nonratio + quant.num.ratio

  ##--- get number of measurements for each quantity
  quant.meas.num = apply(delta, 2, sum)

  ##--- number of quantities with at least a meaurement
  quant.num.nonratio.with.meas = sum(quant.meas.num[quant.names.nonratio] > 0)
  quant.num.ratio.with.meas  = sum(quant.meas.num[quant.names.ratio] > 0)
  quant.num.nonunitarity.with.meas = quant.num.nonratio.with.meas + quant.num.ratio.with.meas

  ##
  ## PDG fit special determinations
  ##

  ##
  ## decay modes used in HFAG but referred to a single node in the PDG
  ##
  quant.names.merged.in.pdg = c(
    "Gamma168", # K phi(K+K-) 
    "Gamma169", # K phi(KSKL)
    "Gamma910", # G(2pi- pi+ eta(3pi0) nu(tau) (ex. K0)) / G(total)
    "Gamma911", # G(pi- 2pi0 eta(pi+ pi- pi0) nu(tau)) / G(total)
    "Gamma930", # G(2pi- pi+ eta(pi+ pi- pi0) nu(tau) (ex. K0)) / G(total)
    "Gamma944", # G(2pi- pi+ eta(gamma gamma) nu(tau) (ex. K0)) / G(total)
    NULL
    )

  quant.names.nonratio.pdg = setdiff(alucomb2.nonratio.quant(combination), quant.names.merged.in.pdg)
  quant.names.ratio.pdg = setdiff(alucomb2.ratio.quant(combination), quant.names.merged.in.pdg)

  quant.num.nonratio.pdg = length(quant.names.nonratio.pdg)
  quant.num.ratio.pdg = length(quant.names.ratio.pdg)
  quant.num.nonunitarity.pdg = quant.num.nonratio.pdg + quant.num.ratio.pdg

  ##
  ## numbers of PDG quantities in the fit with at least one measurement
  ## simulate clustering of modes from HFAG to PDG
  ##

  quant.names.nonratio.pdg = alucomb2.nonratio.quant(combination)
  quant.names.nonratio.pdg = quant.names.nonratio.pdg[quant.meas.num[quant.names.nonratio.pdg] > 0]
  quant.names.nonratio.pdg = gsub("Gamma910", "3pi eta", quant.names.nonratio.pdg)
  quant.names.nonratio.pdg = gsub("Gamma930", "3pi eta", quant.names.nonratio.pdg)
  quant.names.nonratio.pdg = gsub("Gamma944", "3pi eta", quant.names.nonratio.pdg)
  quant.names.nonratio.pdg = gsub("Gamma911", "pi2pi0 eta", quant.names.nonratio.pdg)
  quant.names.nonratio.pdg = unique(quant.names.nonratio.pdg)
  
  quant.names.ratio.pdg = alucomb2.ratio.quant(combination)
  quant.names.ratio.pdg = quant.names.ratio.pdg[quant.meas.num[quant.names.ratio.pdg] > 0]
  quant.names.ratio.pdg = gsub("Gamma910", "3pi eta", quant.names.ratio.pdg)
  quant.names.ratio.pdg = gsub("Gamma930", "3pi eta", quant.names.ratio.pdg)
  quant.names.ratio.pdg = gsub("Gamma944", "3pi eta", quant.names.ratio.pdg)
  quant.names.ratio.pdg = gsub("Gamma911", "pi2pi0 eta", quant.names.ratio.pdg)
  quant.names.ratio.pdg = unique(quant.names.ratio.pdg)

  quant.num.nonratio.with.meas.pdg = length(quant.names.nonratio.pdg)
  quant.num.ratio.with.meas.pdg  = length(quant.names.ratio.pdg)
  quant.num.nonunitarity.with.meas.pdg = quant.num.nonratio.with.meas.pdg + quant.num.ratio.with.meas.pdg

  ##
  ## write tex defs of some quantities
  ##
  tex.defs = c(
    ##--- unitarity check
    alurep.tex.cmd.short("UnitarityResid",
                         alurep.tex.val.err.auto(quant.val["Gamma998"], quant.err["Gamma998"], perc=TRUE)),
    ##--- measurements
    alurep.tex.cmd.short("MeasNum", as.character(meas.num)),

    alurep.tex.cmd.short("QuantNum", as.character(quant.num.nonunitarity)),
    alurep.tex.cmd.short("QuantNumNonRatio", as.character(quant.num.nonratio)),
    alurep.tex.cmd.short("QuantNumRatio", as.character(quant.num.ratio)),

    alurep.tex.cmd.short("QuantNumWithMeas", as.character(quant.num.nonunitarity.with.meas)),
    alurep.tex.cmd.short("QuantNumNonRatioWithMeas", as.character(quant.num.nonratio.with.meas)),
    alurep.tex.cmd.short("QuantNumRatioWithMeas", as.character(quant.num.ratio.with.meas)),

    alurep.tex.cmd.short("QuantNumPdg", as.character(quant.num.nonunitarity.pdg)),
    alurep.tex.cmd.short("QuantNumNonRatioPdg", as.character(quant.num.nonratio.pdg)),
    alurep.tex.cmd.short("QuantNumRatioPdg", as.character(quant.num.ratio.pdg)),

    alurep.tex.cmd.short("QuantNumWithMeasPdg", as.character(quant.num.nonunitarity.with.meas.pdg)),
    alurep.tex.cmd.short("QuantNumNonRatioWithMeasPdg", as.character(quant.num.nonratio.with.meas.pdg)),
    alurep.tex.cmd.short("QuantNumRatioWithMeasPdg", as.character(quant.num.ratio.with.meas.pdg)),

    ##--- basis quantities
    alurep.tex.cmd.short("IndepQuantNum", as.character(quant.num - constr.num)),
    alurep.tex.cmd.short("BaseQuantNum", as.character(length(alucomb2.base.quant(combination)))),

    ##--- quantities in the unitarity constraint sum
    alurep.tex.cmd.short("UnitarityQuantNum", as.character(length(alucomb2.unitarity.constr(combination)))),

    ##
    ## constraints used to relate measurements to basis quantities
    ## traditionally, the constraints used for unitarity and unitarity residual are not counted
    ## this number is consistent with "QuantNum", "Dof" and "MeasNum", both "QuantNum" and "ConstrNum" are two less
    ##
    alurep.tex.cmd.short("ConstrNum", as.character(constr.num - (quant.num - quant.num.nonunitarity))),
    alurep.tex.cmd.short("ConstrNumPdg", as.character(constr.num - (quant.num - quant.num.nonunitarity.pdg))),

    ##--- chisq, dof, chisqprob
    alurep.tex.cmd.short("Chisq", sprintf("%.1f", chisq)),
    alurep.tex.cmd.short("Dof", as.character(dof)),
    alurep.tex.cmd.short("ChisqProb", alurep.tex.val.auto(chisq.prob, perc=TRUE)),
    alurep.tex.cmd.short("ChisqProbRound", sprintf("%.0f\\%%", round(chisq.prob*100)))
    )
  cat(tex.defs, sep="\n", file=fname, append=TRUE)
  cat("file '", fname, "', initial defs\n", sep="")

  ##
  ## define measurements
  ##
  rc = get.tex.meas.defs()
  cat(rc, file=fname, append=TRUE)
  cat("\n", file=fname, append=TRUE)
  cat("file '", fname, "', measurement description defs\n", sep="")

  ##
  ## write tex macro containing BR fit data
  ##
  quant.names = combination$combine
  quant.names = setdiff(quant.names, "GammaAll")
  rc = get.tex.table(quant.names)

  cat(rc$defs, sep="\n", file=fname, append=TRUE)
  cat("file '", fname, "', quantity - measurements defs\n", sep="")

  tex.all.tau.br.val = alurep.tex.cmd("BrVal", rc$table)
  cat(tex.all.tau.br.val, sep="\n", file=fname, append=TRUE)
  cat("file '", fname, "', quantity - measurements table\n", sep="")

  ##
  ## write tex macro containing all measurements by reference
  ##
  rc = get.tex.meas.by.ref()
  cat(rc$defs, sep="\n", file=fname, append=TRUE)
  cat("file '", fname, "', BR meas by ref definitions\n", sep="")
  ##
  tex.meas.paper = alurep.tex.cmd("MeasPaper", rc$table)
  cat(tex.meas.paper, sep="\n", file=fname, append=TRUE)
  cat("file '", fname, "', BR meas by ref table content\n", sep="")
  rm(rc)

  ##
  ## write text macro containing all strange BR values and refs
  ##
  gamma110.names = names(combination$constr.all.comb$Gamma110.c)
  gamma110.names = setdiff(gamma110.names, "Gamma110")
  tex.tau.br.strange.val = alurep.tex.cmd("BrStrangeVal", get.tex.table.simple(gamma110.names, 4, -2))
  cat(tex.tau.br.strange.val, sep="\n", file=fname, append=TRUE)
  cat("file '", fname, "', BR strange table content\n", sep="")

  ##
  ## write text macro containing total strange BR
  ##
  gamma110.names = "Gamma110"
  tex.tau.br.strange.tot.val = alurep.tex.cmd("BrStrangeTotVal", get.tex.table.simple(gamma110.names, 4, -2))
  cat(tex.tau.br.strange.tot.val, sep="\n", file=fname, append=TRUE)
  cat("file '", fname, "', BR strange tot table content\n", sep="")

  ##
  ## write text macro containing all modes in unitarity constraint
  ##
  constr.comb = alucomb2.unitarity.constr(combination)
  tex.tau.unitarity.quants = alurep.tex.cmd("UnitarityQuants", get.tex.table.constr(constr.comb, 4, -2, 4, 0))
  cat(tex.tau.unitarity.quants, sep="\n", file=fname, append=TRUE)
  cat("file '", fname, "', unitarity quantities\n", sep="")

  ##
  ## write text macro containing all basis modes
  ##
  gammaAll.names = alucomb2.base.quant(combination)
  tex.tau.base.quants = alurep.tex.cmd("BaseQuants", get.tex.table.simple(gammaAll.names, 4, -2))
  cat(tex.tau.base.quants, sep="\n", file=fname, append=TRUE)
  cat("file '", fname, "', basis quantities\n", sep="")

  ##
  ## write text macro containing correlation of basis quantities
  ##
  tex.all.tau.br.corr = alurep.tex.cmd("BrCorr", get.tex.base.nodes.corr())
  cat(tex.all.tau.br.corr, sep="\n", file=fname, append=TRUE)
  cat("file '", fname, "', BR correlations table content\n", sep="")

  ##
  ## write TeX definitions for all constraint equations
  ##
  tex.constr.defs = get.tex.constraint.defs()
  cat(tex.constr.defs, sep="\n", file=fname, append=TRUE)
  cat("file '", fname, "', constraint definitions\n", sep="")
  
  ##
  ## write text macro containing constraint equations
  ##
  tex.constr.val = alurep.tex.cmd("ConstrEqs", get.tex.constraint.equations())
  cat(tex.constr.val, sep="\n", file=fname, append=TRUE)
  cat("file '", fname, "', constraint table content\n", sep="")

  ##
  ## measurements by collaboration
  ##
  tex.num.meas.per.collab = get.tex.meas.by.collab()
  cat(tex.num.meas.per.collab, sep="\n", file=fname, append=TRUE)
  cat("file '", fname, "', measurements per collaboration\n", sep="")
}

## ////////////////////////////////////////////////////////////////////////////
## code

args = commandArgs(TRUE)
if (length(args) == 1) {
  rc = mkreport(fname = args[1])
} else {
  mkreport(fname = "average.rdata")
}
