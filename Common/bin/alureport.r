## ////////////////////////////////////////////////////////////////////////////
##
## utility functions for HFAG-tau
##

require(methods, quietly=TRUE)
require(stringr, quietly=TRUE)

## ////////////////////////////////////////////////////////////////////////////
##
## object for doing string translations
##
TrStr = setRefClass("TrStr",
  fields = list(
    .table = "numeric"
    ),
  methods=list(
    initialize = function(set1=character(), set2=character()) {
      ## callSuper(...)
      table = 1:256
      if (length(set1) != 0 || length(set2) != 0) {
        table[as.integer(charToRaw(set1))+1] = as.integer(charToRaw(set2))+1
      }
      .self$.table = table
      .self
    })
  )

##--- translate a string using stored table
rc = TrStr$methods(
  tr = function(str) {
    rawToChar(as.raw(.table[as.integer(charToRaw(str))+1]-1))
  })

##--- translate a string using stored table, delimit numbers by "N"
rc = TrStr$methods(
  trN = function(str) {
    str = gsub("(\\d+)", "N\\1N", str)
    rawToChar(as.raw(.table[as.integer(charToRaw(str))+1]-1))
  })

##--- specialized string translation class for conversion of alphanumeric strings to TeX keywords
TrStr.num2tex = setRefClass("TrStr.num2tex", contains="TrStr",
  methods=list(
    initialize = function() {
      callSuper("0123456789_", "zothfvsneiU")
    })
  )

##
## round data frame to print it in a reproductible way
##
round.data.frame <- function(x, digits = 0) {
  x.names = rownames(x)
  rc = data.frame(lapply(x, function(y) if(is.numeric(y)) round(y, digits) else y))
  rownames(rc) = x.names
  rc
}

##
## reporting functions
##

##
## return tex label such as \Gamma_1 or \frac{\Gamma_1}{\Gamma_2}
##
alurep.gamma.texlabel = function(str) {
  str = gsub("GammaAll", "\\\\Gamma_{\\\\text{All}}", str, perl=TRUE)
  str = gsub("Gamma(\\d+)by(\\d+)", "\\\\frac{\\\\Gamma_{\\1}}{\\\\Gamma_{\\2}}", str, perl=TRUE)
  str = gsub("Gamma(\\d+)", "\\\\Gamma_{\\1}", str, perl=TRUE)
  return(str)
}

##
## get tex description of a quantity in alucomb2 structure
##
alurep.get.texdescr.nv = function(quant, descr, texdescr) {
  repeat {
    if (quant == "Gamma998" || quant == "GammaAll") {
      break
    }

    if (!is.null(texdescr) && texdescr != "") {
      ##+++ check for ratio of BR
      texdescr = paste("\\BRF{\\tau^-}{", texdescr, "}", sep="")
      texdescr = gsub("mathrm", "text", texdescr, fixed=TRUE)
      break
    }
    
    if (is.null(descr) || descr == "") {
      texdescr = ""
      break
    }

    ##--- transform text description into TeX
    descr = sub("\\s+/\\s+G\\(total\\)\\s*", "", descr)
    descr = gsub("ex[.]\\s+K\\(S\\)0 --> pi- pi[+]", "ex. K0", descr, perl=TRUE)
    descr = gsub("\\s*\\(``\\d-prong''\\)", "", descr, perl=TRUE)
    descr = gsub("-->", "\\to", descr, fixed=TRUE)
    descr = gsub("([+-])", "^\\1", descr)
    descr = gsub("\\^-prong", "\\\\text{-prong}", descr)
    descr = gsub("\\s+or\\s+", "\\\\,\\\\text{or}\\\\,", descr)

    ##--- convert <decay mode> into \BRF{\tau}{<decay mode>}
    descr = gsub("G(\\(((?:[^()]++|(?1))+)*+\\))", "\\\\BRF{tau^-}{\\2}", descr, perl=TRUE)

    ##--- transform text description into TeX
    descr = gsub("([<=>]+)\\s*(\\d+\\s*)(neutrals|K0|K\\(L\\)0|K\\(S\\)0|pi0)\\s*nu\\((mu|tau|e)\\)",
                 "\\1 \\2 \\3\\\\, nu(\\4)", descr)
    descr = gsub("([<=>]+)\\s*(\\d+)", "\\1 \\2\\\\,", descr)
    descr = gsub(">=", "\\\\ge{}", descr)
    descr = gsub("nubar\\((mu|tau|e)\\)", "\\\\bar{nu}_\\1", descr)
    descr = gsub("nu\\((mu|tau|e)\\)", "nu_\\1", descr)
    descr = gsub("pi0", "pi^0", descr)
    descr = gsub("Kstar", "K^*", descr)
    descr = gsub("Kbar0", "\\\\bar{K}^0", descr)
    descr = gsub("K0", "K^0", descr)
    descr = gsub("K\\(S\\)0", "K_S^0", descr)
    descr = gsub("K\\(L\\)0", "K_L^0", descr)
    descr = gsub("(nu|tau|mu|pi|omega|eta|gamma)", "\\\\\\1", descr)
    descr = gsub("\\s+\\(ex[.]", "\\\\;(\\\\text{ex.~}", descr)
    descr = gsub("(neutrals)", "\\\\text{\\1}", descr)
    descr = gsub("(particles|strange|total)", "\\\\text{\\1}", descr)

    ##--- transform division of BR into a TeX fraction
    descr = gsub("^(.*\\S)\\s*/\\s*(\\S.*)$", "\\\\frac{\\1}{\\2}", descr, perl=TRUE)

    texdescr = descr
    break
  }

  return(texdescr)
}
alurep.get.texdescr = Vectorize(alurep.get.texdescr.nv, USE.NAMES=TRUE)

##
## get tex description of a quantity in alucomb2 structure
##
alurep.get.texdescr.quantities = function(quantities) {
  quant.names = names(quantities)
  quant.texdescr = sapply(quantities, function(x) {if (is.null(x$texdescr)) {""} else {x$texdescr}})
  quant.descr = sapply(quantities, function(x) {if (is.null(x$descr)) {""} else {x$descr}})
  alurep.get.texdescr(quant.names, quant.descr, quant.texdescr)
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
## return measurement value, stat, syst original values
## all properly formatted for latex printing, separately for
## - quant = val +- stat +- syst
## - val
## - stat
## - syst
##
alurep.tex.meas.val.card.fields = function(meas) {
  quant.tex = character(0)
  ee.data = character(0)

  val.txt = attr(meas$value.orig, "input")
  quant.tex = c(quant.tex, val.txt)
  ee.data = c(ee.data, val.txt)

  if (attr(meas$stat, "input") != "") {
    stat.txt = paste("\\pm", attr(meas$stat, "input"))
    ee.data = c(ee.data, attr(meas$stat, "input"))
  } else {
    stat.txt = paste("{}^{", attr(meas$stat.p, "input"), "}_{", attr(meas$stat.n, "input"), "}", sep="")
    ee.data = c(ee.data, attr(meas$stat.p, "input"), attr(meas$stat.n, "input"))
  }
  quant.tex = c(quant.tex, stat.txt)

  if (attr(meas$syst, "input") != "") {
    syst.txt = attr(meas$syst, "input")
    if (syst.txt != "0") {
      quant.tex = c(quant.tex, paste("\\pm", syst.txt))
      ee.data = c(ee.data, syst.txt)
    }
  } else {
    syst.txt = paste("{}^{", attr(meas$syst.p, "input"), "}_{", attr(meas$syst.n, "input"), "}", sep="")
    quant.tex = c(quant.tex, syst.txt)
    ee.data = c(ee.data, attr(meas$syst.p, "input"), attr(meas$syst.n, "input"))
  }

  order = gsub("^[^e]+(|e[+]?([-])?0*(\\d+))$", "\\2\\3", ee.data, perl=TRUE)
  order = ifelse(order == "", 0, as.numeric(order))  
  
  quant.tex = paste(quant.tex, collapse=" ")
  if (max(order) == min(order)) {
    if (max(order) != 0) {
      quant.tex = gsub("(\\S+)e[+-]*\\d+", "\\1", quant.tex, perl=TRUE, ignore.case=TRUE)
      quant.tex = paste("(", quant.tex, ") \\cdot 10^{", max(order), "}")
    }
  }
  
  val.tex = gsub("e[+]?([-])?0*(\\d+)", "\\\\cdot 10^{\\1\\2}", val.txt, ignore.case=TRUE)
  stat.tex = gsub("e[+]?([-])?0*(\\d+)", "\\\\cdot 10^{\\1\\2}", stat.txt, ignore.case=TRUE)
  syst.tex = gsub("e[+]?([-])?0*(\\d+)", "\\\\cdot 10^{\\1\\2}", syst.txt, ignore.case=TRUE)
  quant.tex = gsub("e[+]?([-])?0*(\\d+)", "\\\\cdot 10^{\\1\\2}", quant.tex, ignore.case=TRUE)

  val.tex = gsub("%", "\\%", val.txt, fixed=TRUE)
  stat.tex = gsub("%", "\\%", stat.txt, fixed=TRUE)
  syst.tex = gsub("%", "\\%", syst.txt, fixed=TRUE)
  quant.tex = gsub("%", "\\%", quant.tex, fixed=TRUE)

  return(list(quant=quant.tex, val=val.tex, stat=stat.tex, syst=syst.tex))
}

##
## return value +- stat +- syst original values
## convert exponential format for latex
##
alurep.tex.meas.val.card = function(meas) {
  rc = alurep.tex.meas.val.card.fields(meas)
  return(paste("\\ensuremath{", rc$quant, "}", sep=""))
}

##
## for the spec vector of measurement, compute the
## appropriate precision and power-of-ten order
##

alurep.precision.order = function(vals, perc=FALSE, signif=4, signif.min=2) {
  vals = vals[vals != 0]
  if (length(vals) == 0) {
    vals = 1
  }
  vals.str = sprintf("%.4g", vals)
  order.str = gsub("^[^eE]+(|[eE]((-[0-9]+)+|[+]([0-9]+)))$", "\\2",
    sprintf("%.4e", vals), perl=TRUE)
  order = ifelse(order.str == "", 0, as.numeric(order.str))
  order.max = max(order)
  order.min = min(order)  
  precision = signif.min - 1 + 1 * max(signif-signif.min, order.max - order.min)
  ## print(cbind(order.max, order.min, precision))
  if (order.max == -1) {
    if (!perc) {
      order.max = 0
      precision = precision+1
    } else {
      order.max = -2
      precision = precision-1
    }
  } else if (order.max == -3) {
    order.max = -2
    precision = precision+1
  } else if (order.max == 1) {
    order.max = 0
    precision = precision-1
  } else if (order.max == 2) {
    order.max = 0
    precision = precision-2
  }
  else {
  }
  return(list(precision=precision, order=order.max))
}

##
## return numeric value formatted val in a string
## according to the specified precision and power-of-ten order
##
alurep.tex.val.prec.ord = function(val, precision, order, width=0, perc=FALSE) {
  val = val/10^order
  if (order == 0) {
    rc = sprintf(paste("%", width, ".", precision, "f", sep=""), val)
  } else if (perc && order == -2) {
    ##--- will act depending on  alurep.precision.order()
    rc = sprintf(paste("%", width, ".", precision, "f\\%%", sep=""), val)
  } else {
    rc = sprintf(paste("%", width, ".", precision, "f\\cdot 10^{%d}", sep=""), val, order)
  }
  return(rc)
}

##
## return numeric value formatted val in a string
## use automatic precision and order of magnitude
##
alurep.tex.val.auto = function(val, width=0, perc=FALSE) {
  rc = alurep.precision.order(val, perc=perc)
  precision = rc$precision
  order = rc$order
  return(alurep.tex.val.prec.ord(val, precision, order, width=width, perc=perc))
}

##
## return quantity formatted val +- stat in a string
## according to the specified precision and power-of-ten order
##
alurep.tex.val.err.prec.ord = function(val, err, precision, order, width=0, perc=FALSE) {
  val = val/10^order
  err = err/10^order
  if (order == 0) {
    rc = sprintf(paste("%", width, ".", precision, "f \\pm %", width, ".", precision, "f", sep=""),
      val, err)
  } else if (perc && order == -2) {
    rc = sprintf(paste("(%", width, ".", precision, "f \\pm %", width, ".", precision, "f)\\%%", sep=""),
      val, err)
  } else {
    rc = sprintf(paste("(%", width, ".", precision, "f \\pm %", width, ".", precision, "f) \\cdot 10^{%d}", sep=""),
      val, err, order)
  }
  return(rc)
}

##
## return quantity formatted val +- stat in a string
## according to the specified precision and power-of-ten order
## do not add the x10^order part
##
alurep.tex.val.err.prec.ord.noee = function(val, err, precision, order, width=0, perc=FALSE) {
  val = val/10^order
  err = err/10^order
  rc = sprintf(paste("%", width, ".", precision, "f \\pm %", width, ".", precision, "f", sep=""), val, err)
  return(rc)
}

##
## return quantity formatted val +- stat in a string
## according to the self-determined optimal precision and power-of-ten order
##
alurep.tex.val.err.auto = function(val, err, width=0, perc=FALSE) {
  rc = alurep.precision.order(c(val, err), perc=perc)
  precision = rc$precision
  order = rc$order
  return(alurep.tex.val.err.prec.ord(val, err, precision, order, width=width, perc=perc))
}

##
## return measurement formatted val +- stat +- syst in a string
## according to the specified precision and power-of-ten order
##
alurep.tex.meas.val = function(meas, precision, order, width=0, perc=FALSE) {
  norm = 10^order
  str.val = sprintf(paste("%", width, ".", precision, "f", sep=""), meas$value.orig/norm)
  if (meas$stat.p == -meas$stat.n) {
    str.stat = sprintf(paste("\\pm %", width, ".", precision, "f", sep=""), meas$stat/norm)
  } else {
    str.stat = sprintf(paste("{}^{%+", width, ".", precision, "f}_{%+", width, ".", precision, "f}", sep=""),
      meas$stat.p/norm, meas$stat.n/norm)
  }
  if (meas$syst.p == -meas$syst.n) {
    str.syst = sprintf(paste("\\pm %", width, ".", precision, "f", sep=""), meas$syst.orig/norm)
  } else {
    str.syst = sprintf(paste("{}^{%+", width, ".", precision, "f}_{%+", width, ".", precision, "f}", sep=""),
      meas$syst.p/norm, meas$syst.n/norm)
  }
  str.meas = paste(str.val, str.stat, str.syst)
  if (order == 0) {
    rc = paste("", str.meas, "", sep="")
  } else {
    rc = sprintf(paste("(", str.meas, ") \\cdot 10^{%d} ", sep=""), order)
  }
  return(paste("\\ensuremath{", rc, "}", sep=""))
}

##
## for the spec. quantity return the associated measurements
##
alurep.meas.quant = function(quant.name, delta) {
  delta.names = delta[, quant.name]
  meas.names = delta.names[delta.names!=0]
  return(names(meas.names))
}

##
## compute appropriate precision and order of magnitude for printing a measurement
##
alurep.precision.order.meas = function(meas, perc=FALSE) {
  vals = c(meas$value.orig, meas$stat.p, meas$stat.n, meas$syst.p, meas$syst.n)
  return(alurep.precision.order(vals, perc=perc))
}

##
## return measurement formatted val +- stat +- syst in a string
## according to the specified precision and power-of-ten order
##
alurep.tex.meas.val.auto = function(meas, width=0, perc=FALSE) {
  rc = alurep.precision.order.meas(meas, perc=perc)
  precision = rc$precision
  order = rc$order
  return(alurep.tex.meas.val(meas, precision, order, width, perc=perc))
}

##
## get all measurements related to a quantity
## compute appropriate precision and order of magnitude
##
alurep.precision.order.quant = function(quant.name, perc=FALSE, with.meas=FALSE) {
  vals = c(quant.val[quant.name], quant.err[quant.name])
  if (with.meas) {
    meas.names = alurep.meas.quant(quant.name, delta)
    meas.vals = unlist(lapply(measurements[meas.names], function(m) c(m$value.orig, m$stat.p, m$stat.n, m$syst.p, m$syst.n)))
    vals = c(vals, meas.vals)
  }
  return(alurep.precision.order(vals, perc=perc))
}

##
## substitute parameter labels with TeX printable names
##
alurep.subst.params = function(str) {
  str = gsub("BR_eta_2gam", "\\Gamma_{\\eta\\to\\gamma\\gamma}", str, fixed=TRUE)
  str = gsub("BR_eta_neutral", "\\Gamma_{\\eta\\to\\text{neutral}}", str, fixed=TRUE)
  str = gsub("BR_eta_3piz", "\\Gamma_{\\eta\\to3\\pi^0}", str, fixed=TRUE)
  str = gsub("BR_eta_pimpippiz", "\\Gamma_{\\eta\\to\\pi^+\\pi^-\\pi^0}", str, fixed=TRUE)
  str = gsub("BR_eta_charged", "\\Gamma_{\\eta\\to\\text{charged}}", str, fixed=TRUE)
  str = gsub("BR_KS_2piz", "\\Gamma_{K_S\\to\\pi^0\\pi^0}", str, fixed=TRUE)
  str = gsub("BR_KS_pimpip", "\\Gamma_{K_S\\to\\pi^+\\pi^-}", str, fixed=TRUE)
  str = gsub("BR_om_pimpippiz", "\\Gamma_{\\omega\\to\\pi^+\\pi^-\\pi^0}", str, fixed=TRUE)
  str = gsub("BR_om_pimpip", "\\Gamma_{\\omega\\to\\pi^+\\pi^-}", str, fixed=TRUE)
  str = gsub("BR_om_pizgamma", "\\Gamma_{\\omega\\to\\pi^0\\gamma}", str, fixed=TRUE)
  str = gsub("BR_phi_KmKp", "\\Gamma_{\\phi\\to K^+K^-}", str, fixed=TRUE)
  str = gsub("BR_phi_KSKL", "\\Gamma_{\\phi\\to K_S K_L}", str, fixed=TRUE)
  str = gsub("BRA_Kz_KS_KET", "\\Gamma_{<K^0|K_S>}", str, fixed=TRUE)
  str = gsub("BRA_Kz_KL_KET", "\\Gamma_{<K^0|K_L>}", str, fixed=TRUE)
  str = gsub("BRA_Kzbar_KS_KET", "\\Gamma_{<\\bar{K}^0|K_S>}", str, fixed=TRUE)
  str = gsub("BRA_Kzbar_KL_KET", "\\Gamma_{<\\bar{K}^0|K_L>}", str, fixed=TRUE)
  str = gsub("BRA_KzKzbar_KLKL_KET_by_BRA_KzKzbar_KSKS_KET", "R_{0\\bar{0}SS/LL}", str, fixed=TRUE)

  str = gsub("BR_f1_2pizpippim", "\\Gamma_{f_1\\to\\pi^+\\pi^-2\\pi^0}", str, fixed=TRUE)
  str = gsub("BR_f1_2pip2pim", "\\Gamma_{f_1\\to2\\pi^+2\\pi^-}", str, fixed=TRUE)
}

##
## get left and right elements of a constraint equation as TeX expressions
##
alurep.tex.constraint.zeroval = function(constr.str) {
  ##--- replace "0 = -x + y" into "x = y"
  rc = str_match(unlist(constr.str), "-([[:alnum:]]+)\\s+[+]\\s+(.*\\S)\\s*$")
  constr.left = ifelse(!is.na(rc[,2]), rc[,2], "")
  constr.right = ifelse(!is.na(rc[,3]), rc[,3], "")
  ##--- remove outer braces
  constr.right = gsub("^\\s*(\\(((?:[^()]++|(?1))+)*+\\))\\s*$", "\\2", constr.right, perl=TRUE)

  ##--- force GammaXbyY constraint to be GammaX/GammaY
  constr.right = ifelse(regexpr("^Gamma(\\d+)by(\\d+)$", constr.left, perl=TRUE) != -1, constr.left, constr.right)

  ##--- convert Gamma# notation for TeX
  constr.left = alurep.gamma.texlabel(constr.left)
  constr.right = alurep.gamma.texlabel(constr.right)
  return(list(left=constr.left, right=constr.right))
}
alurep.tex.constraint.nonzeroval = function(constr.val, constr.str) {
  constr.left = as.character(constr.val)
  constr.right = constr.str
  ##--- convert Gamma# notation for TeX
  constr.right = alurep.gamma.texlabel(constr.right)
  return(list(left=constr.left, right=constr.right))
}
alurep.tex.constraint = function(constr.val, constr.str) {
  sel.zero = which(constr.val == 0)
  sel.nonzero = which(constr.val != 0)
  rc.zero = alurep.tex.constraint.zeroval(constr.str[sel.zero])
  rc.nonzero = alurep.tex.constraint.nonzeroval(constr.val[sel.nonzero], constr.str[sel.nonzero])

  constr.left = character(0)
  constr.right = character(0)
  constr.left[c(sel.zero, sel.nonzero)] = c(rc.zero$left, rc.nonzero$left)
  constr.right[c(sel.zero, sel.nonzero)] = c(rc.zero$right, rc.nonzero$right)
  names(constr.left) = names(constr.val)
  names(constr.right) = names(constr.val)
  
  ##--- split long eqs
  constr.right.split = gsub("(([^+*]+[+*]){6}[^+]+)[+]", "\\1 \\\\\\\\ \n  {}& +", constr.right, perl=TRUE)
  
  constr.right = gsub("*", "\\cdot{}", constr.right, fixed=TRUE)
  constr.right = alurep.subst.params(constr.right)
  constr.right.split = gsub("*", "\\cdot{}", constr.right.split, fixed=TRUE)
  constr.right.split = alurep.subst.params(constr.right.split)
  
  return(list(left = constr.left, right=constr.right, right.split=constr.right.split))
}

##
## get the inspire bibtex key from the measurement tags
##
get.reference = function(tags) {
  bib.pdg.table = list(
    "BARTELT.96" = "Bartelt:1996iv",
    "DEL-AMO-SANCHEZ.11E" = "delAmoSanchez:2010pc",
    "AUBERT.10B" = "Aubert:2009ag",
    "AUBERT.10F" = "Aubert:2009qj",
    "HAYASAKA.10" = "Hayasaka:2010np",
    "LEE.10" = "Lee:2010tc",
    "LEES.10A" = "not found",
    "MIYAZAKI.10" = "Miyazaki:2009wc",
    "MIYAZAKI.10A" = "Miyazaki:2010qb",
    "AUBERT.09AK" = "Aubert:2009ra",
    "AUBERT.09D" = "Aubert:2009ys",
    "AUBERT.09W" = "Aubert:2009ap",
    "GROZIN.09A" = "Grozin:2008nw",
    "INAMI.09" = "Inami:2008ar",
    "MIYAZAKI.09" = "Miyazaki:2008mw",
    "AUBERT.08" = "Aubert:2007mh",
    "AUBERT.08AE" = "Aubert:2008nj",
    "AUBERT.08K" = "Aubert:2007kx",
    "FUJIKAWA.08" = "Fujikawa:2008ma",
    "HAYASAKA.08" = "Hayasaka:2007vc",
    "MIYAZAKI.08" = "Miyazaki:2007zw",
    "NISHIO.08" = "Nishio:2008zx",
    "ANASHIN.07" = "Anashin:2007zz",
    "AUBERT.07AP" = "Aubert:2007jh",
    "AUBERT.07BK" = "Aubert:2007pw",
    "AUBERT.07I" = "Aubert:2006cz",
    "BELOUS.07" = "Abe:2006vf",
    "EPIFANOV.07" = "Epifanov:2007rf",
    "MIYAZAKI.07" = "Miyazaki:2007jp",
    "ABDALLAH.06A" = "Abdallah:2003cw",
    "AUBERT.06C" = "Aubert:2005wa",
    "AUBERT,B.06" = "Aubert:2006hw",
    "INAMI.06" = "Inami:2006vd",
    "MIYAZAKI.06" = "Miyazaki:2005ng",
    "MIYAZAKI.06A" = "Miyazaki:2006sx",
    "YUSA.06" = "Yusa:2006qq",
    "ARMS.05" = "Arms:2005qg",
    "AUBERT,B.05A" = "Aubert:2005ye",
    "AUBERT,B.05F" = "Aubert:2005nx",
    "AUBERT,B.05W" = "Aubert:2005waa",
    "AUBERT,BE.05D" = "Aubert:2005tp",
    "ENARI.05" = "Enari:2005gc",
    "HAYASAKA.05" = "Hayasaka:2005xw",
    "SCHAEL.05C" = "Schael:2005am",
    "ABBIENDI.04J" = "Abbiendi:2004xa",
    "ABDALLAH.04K" = "Abdallah:2003xd",
    "ABDALLAH.04T" = "Abdallah:2003yq",
    "ABE.04B" = "Abe:2003sx",
    "ACHARD.04G" = "Achard:2004jj",
    "AUBERT.04J" = "Aubert:2003pc",
    "ENARI.04" = "Enari:2004ax",
    "YUSA.04" = "Yusa:2004gm",
    "ABBIENDI.03" = "Abbiendi:2002jw",
    "BRIERE.03" = "Briere:2003fr",
    "HEISTER.03F" = "Heister:2002ik",
    "INAMI.03" = "Inami:2002ah",
    "CHEN.02C" = "not found",
    "ABBIENDI.01J" = "Abbiendi:2000ee",
    "ABREU.01M" = "Abreu:2001wc",
    "ACCIARRI.01F" = "Acciarri:2001sg",
    "ACHARD.01D" = "Achard:2001pk",
    "ANASTASSOV.01" = "Anastassov:2000xu",
    "HEISTER.01E" = "Heister:2001me",
    "ABBIENDI.00A" = "Abbiendi:2000kz",
    "ABBIENDI.00C" = "Abbiendi:1999pm",
    "ABBIENDI.00D" = "Abbiendi:1999cq",
    "ABREU.00L" = "Abreu:2000sg",
    "ACCIARRI.00B" = "Acciarri:2000vq",
    "AHMED.00" = "not found",
    "ALBRECHT.00" = "Albrecht:2000yg",
    "ASNER.00" = "Asner:1999kj",
    "ASNER.00B" = "Asner:2000nx",
    "BERGFELD.00" = "Bergfeld:1999yh",
    "BROWDER.00" = "Browder:1999fr",
    "EDWARDS.00A" = "Edwards:1999fj",
    "GONZALEZ-SPRINBERG.00" = "GonzalezSprinberg:2000mk",
    "ABBIENDI.99H" = "Abbiendi:1998cx",
    "ABREU.99X" = "Abreu:1999rb",
    "ACKERSTAFF.99D" = "Ackerstaff:1998yk",
    "ACKERSTAFF.99E" = "Ackerstaff:1998ia",
    "BARATE.99K" = "Barate:1999hi",
    "BARATE.99R" = "Barate:1999hj",
    "BISHAI.99" = "Bishai:1998gf",
    "GODANG.99" = "Godang:1999ge",
    "RICHICHI.99" = "Richichi:1998bc",
    "ACCIARRI.98C" = "Acciarri:1998zc",
    "ACCIARRI.98E" = "Acciarri:1998iv",
    "ACCIARRI.98R" = "Acciarri:1998as",
    "ACKERSTAFF.98M" = "Ackerstaff:1997tx",
    "ACKERSTAFF.98N" = "Ackerstaff:1998mt",
    "ALBRECHT.98" = "Albrecht:1997gn",
    "BARATE.98" = "Barate:1997ma",
    "BARATE.98E" = "Barate:1997tt",
    "BLISS.98" = "Bliss:1997iq",
    "ABE.97O" = "Abe:1997dy",
    "ACKERSTAFF.97J" = "Ackerstaff:1997is",
    "ACKERSTAFF.97L" = "Ackerstaff:1996gy",
    "ACKERSTAFF.97R" = "Ackerstaff:1997dv",
    "ALEXANDER.97F" = "Alexander:1997bv",
    "AMMAR.97B" = "Ammar:1996xh",
    "ANASTASSOV.97" = "Anastassov:1996tc",
    "ANDERSON.97" = "Anderson:1997kb",
    "AVERY.97" = "not found",
    "BARATE.97I" = "Barate:1997hw",
    "BARATE.97R" = "Barate:1997be",
    "BERGFELD.97" = "Bergfeld:1997zt",
    "BONVICINI.97" = "Bonvicini:1997bw",
    "BUSKULIC.97C" = "Buskulic:1996qs",
    "COAN.97" = "Coan:1997am",
    "EDWARDS.97" = "not found",
    "EDWARDS.97B" = "not found",
    "ESCRIBANO.97" = "Escribano:1996wp",
    "ABREU.96B" = "Abreu:1995xt",
    "ACCIARRI.96H" = "Acciarri:1996ht",
    "ACCIARRI.96K" = "Acciarri:1996kc",
    "ALAM.96" = "Alam:1995mt",
    "ALBRECHT.96E" = "Albrecht:1996gr",
    "ALEXANDER.96D" = "Alexander:1995tw",
    "ALEXANDER.96E" = "Alexander:1996nc",
    "ALEXANDER.96S" = "Alexander:1996um",
    "BAI.96" = "Bai:1995hf",
    "BALEST.96" = "Balest:1996cr",
    "BARTELT.96" = "Bartelt:1996iv",
    "BUSKULIC.96" = "Buskulic:1995ty",
    "BUSKULIC.96C" = "Buskulic:1995rh",
    "COAN.96" = "Coan:1996iu",
    "ABE.95Y" = "Abe:1995yj",
    "ABREU.95T" = "Abreu:1995gr",
    "ABREU.95U" = "Abreu:1995gs",
    "ACCIARRI.95" = "Acciarri:1994vr",
    "ACCIARRI.95F" = "Acciarri:1995kx",
    "AKERS.95F" = "Akers:1994xf",
    "AKERS.95I" = "Akers:1995nh",
    "AKERS.95P" = "Akers:1995vy",
    "AKERS.95Y" = "Akers:1995ry",
    "ALBRECHT.95" = "Albrecht:1994nm",
    "ALBRECHT.95C" = "Albrecht:1995un",
    "ALBRECHT.95G" = "Albrecht:1995ht",
    "ALBRECHT.95H" = "Albrecht:1995wx",
    "BALEST.95C" = "Balest:1995kq",
    "BUSKULIC.95C" = "Buskulic:1994yh",
    "BUSKULIC.95D" = "Buskulic:1994hi",
    "ABREU.94K" = "Abreu:1994fi",
    "AKERS.94E" = "Akers:1994ex",
    "AKERS.94G" = "Akers:1994td",
    "ALBRECHT.94E" = "Albrecht:1994es",
    "ARTUSO.94" = "Artuso:1994ii",
    "BARTELT.94" = "Bartelt:1994hn",
    "BATTLE.94" = "Battle:1994by",
    "BAUER.94" = "Bauer:1993wn",
    "BUSKULIC.94D" = "Buskulic:1993mn",
    "BUSKULIC.94E" = "Buskulic:1994je",
    "BUSKULIC.94F" = "Buskulic:1994jf",
    "GIBAUT.94B" = "Gibaut:1994ik",
    "ADRIANI.93M" = "Adriani:1993gk",
    "ALBRECHT.93C" = "Albrecht:1992ka",
    "ALBRECHT.93G" = "Albrecht:1993fr",
    "BALEST.93" = "not found",
    "BEAN.93" = "Bean:1992hh",
    "BORTOLETTO.93" = "Bortoletto:1993px",
    "ESCRIBANO.93" = "Escribano:1993pq",
    "PROCARIO.93" = "Procario:1992hd",
    "ABREU.92N" = "Abreu:1992gn",
    "ACTON.92F" = "Acton:1992ff",
    "ACTON.92H" = "Acton:1992by",
    "AKERIB.92" = "Akerib:1992hb",
    "ALBRECHT.92D" = "Albrecht:1991rh",
    "ALBRECHT.92K" = "Albrecht:1992uba",
    "ALBRECHT.92M" = "Albrecht:1992td",
    "ALBRECHT.92Q" = "Albrecht:1992pe",
    "AMMAR.92" = "Ammar:1991jc",
    "ARTUSO.92" = "Artuso:1992qu",
    "BAI.92" = "Bai:1992bu",
    "BATTLE.92" = "Battle:1992qv",
    "BUSKULIC.92J" = "Buskulic:1992jj",
    "DECAMP.92C" = "Decamp:1991jp",
    "ADEVA.91F" = "Adeva:1991qq",
    "ALBRECHT.91D" = "Albrecht:1990nc",
    "ALEXANDER.91D" = "Alexander:1991am",
    "ANTREASYAN.91" = "Antreasyan:1990fq",
    "GRIFOLS.91" = "Grifols:1990ha",
    "ABACHI.90" = "Abachi:1989vr",
    "ALBRECHT.90E" = "Albrecht:1990zj",
    "BEHREND.90" = "Behrend:1989qx",
    "BOWCOCK.90" = "Bowcock:1989mq",
    "DELAGUILA.90" = "delAguila:1990jg",
    "GOLDBERG.90" = "Goldberg:1990ep",
    "WU.90" = "Wu:1989hx",
    "ABACHI.89B" = "Abachi:1988gx",
    "BEHREND.89B" = "Behrend:1989wc",
    "JANSSEN.89" = "Janssen:1989wg",
    "KLEINWORT.89" = "Kleinwort:1988sb",
    "ADEVA.88" = "Adeva:1988eq",
    "ALBRECHT.88B" = "Albrecht:1987zf",
    "ALBRECHT.88L" = "Albrecht:1988bt",
    "ALBRECHT.88M" = "Albrecht:1988ik",
    "AMIDEI.88" = "Amidei:1987zj",
    "BEHREND.88" = "Behrend:1987fa",
    "BRAUNSCHWEIG.88C" = "Braunschweig:1988yy",
    "KEH.88" = "Keh:1988gs",
    "TSCHIRHART.88" = "Tschirhart:1987uh",
    "ABACHI.87B" = "Abachi:1987hw",
    "ABACHI.87C" = "Abachi:1987nh",
    "ADLER.87B" = "Adler:1987bf",
    "AIHARA.87B" = "Aihara:1986mw",
    "AIHARA.87C" = "Aihara:1987zm",
    "ALBRECHT.87L" = "not found",
    "ALBRECHT.87P" = "Albrecht:1987fb",
    "BAND.87" = "Band:1987nv",
    "BAND.87B" = "Band:1987bm",
    "BARINGER.87" = "Baringer:1987tr",
    "BEBEK.87C" = "Bebek:1987nh",
    "BURCHAT.87" = "Burchat:1986yp",
    "BYLSMA.87" = "Bylsma:1986zy",
    "COFFMAN.87" = "Coffman:1987da",
    "DERRICK.87" = "Derrick:1987sp",
    "FORD.87" = "Ford:1986zk",
    "FORD.87B" = "Ford:1987ha",
    "GAN.87" = "Gan:1987mr",
    "GAN.87B" = "Gan:1987xs",
    "AIHARA.86E" = "Aihara:1986nj",
    "BARTEL.86D" = "Bartel:1986un",
    "RUCKSTUHL.86" = "Ruckstuhl:1986qg",
    "SCHMIDKE.86" = "Schmidke:1986gp",
    "YELTON.86" = "Yelton:1985sx",
    "ALTHOFF.85" = "Althoff:1984ja",
    "ASH.85B" = "Ash:1985hp",
    "BALTRUSAITIS.85" = "Baltrusaitis:1985fh",
    "BARTEL.85F" = "not found",
    "BEHRENDS.85" = "Behrends:1985pm",
    "BELTRAMI.85" = "Beltrami:1985vb",
    "BERGER.85" = "Berger:1985in",
    "BURCHAT.85" = "Burchat:1985mv",
    "FERNANDEZ.85" = "Fernandez:1984if",
    "MILLS.85" = "Mills:1985mh",
    "AIHARA.84C" = "Aihara:1984hr",
    "BEHREND.84" = "Behrend:1984up",
    "MILLS.84" = "Mills:1984dn",
    "BEHREND.83C" = "not found",
    "SILVERMAN.83" = "Silverman:1982ft",
    "BEHREND.82" = "not found",
    "BLOCKER.82B" = "Blocker:1982rw",
    "BLOCKER.82D" = "not found",
    "FELDMAN.82" = "Feldman:1981md",
    "HAYES.82" = "Hayes:1981bn",
    "BERGER.81B" = "not found",
    "DORFAN.81" = "Dorfan:1980mk",
    "BRANDELIK.80" = "not found",
    "BACINO.79B" = "Bacino:1979fz",
    "BACINO.78B" = "Bacino:1978gb",
    "BRANDELIK.78" = "not found",
    "JAROS.78" = "Jaros:1978um",
    "DAVIER.06" = "Davier:2005xq",
    "RAHAL-CALLOT.98" = "RahalCallot:1998ki",
    "GENTILE.96" = "Gentile:1995ue",
    "WEINSTEIN.93" = "not found",
    "PERL.92" = "Perl:1991gd",
    "PICH.90" = "Pich:1990gh",
    "BARISH.88" = "Barish:1987nj",
    "GAN.88" = "Gan:1987zi",
    "HAYES.88" = "Hayes:1988ut",
    "PERL.80" = "not found",
    "LEES.12X" = "Lees:2012ks",
    "LEES.12Y" = "Lees:2012de",
    "Belous:2013dba" = "Belous:2013dba",
    "RYU.14vpc" = "Ryu:2014vpc"
    )

  bib.conf.table = list(
    "BaBar.ICHEP08" = "Aubert:2008an",
    "BaBar.DPF09" = "Paramesvaran:2009ec",
    "Belle.PHIPSI11" = "SooRyu:2011aa"
    )

  ref.conf = bib.conf.table[[paste(tags[1], tail(tags, -3), collapse=".", sep=".")]]
  if (!is.null(ref.conf)) return(ref.conf)

  tag = paste(tail(tags, -3), collapse=".")
  ref.insp = bib.pdg.table[[tag]]
  if (is.null(ref.insp) || ref.insp=="not found") {
    ref.insp = paste("not found:", paste(tags, collapse="."))
  }
  return(ref.insp)
}
