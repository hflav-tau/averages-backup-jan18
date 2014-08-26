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
## return numeric id for sorting labels like "Gamma5", "Gamma3by5"
## <n>by<m> are sorted after <n> in ascending order ov <m>
##
alurep.gamma.num.id = function(gamma.name) {
  gamma.name = sub("Unitarity", "Gamma1000", gamma.name, fixed=TRUE)
  gamma.name = sub("GammaAll", "Gamma999", gamma.name, fixed=TRUE)
  num1 = as.numeric(str_match(gamma.name, ("^(\\D+)(\\d+[.]?\\d*)"))[,3])
  tmp = str_match(gamma.name, ("^\\D+\\d+(by|)(\\d*)"))[,3]
  num2 = ifelse(tmp=="", 0, as.numeric(tmp))
  return(1000*num1 + num2)
}

##
## reporting functions
##

##
## return tex label such as \Gamma_1 or \frac{\Gamma_1}{\Gamma_2}
##
alurep.gamma.texlabel.nv = function(gamma.name) {
  if (gamma.name == "GammaAll") return("1")
  if (regexpr("^Gamma\\d+(|by\\d+)$", gamma.name, perl=TRUE) == -1) return("")
  str = str_match(gamma.name, "Gamma(\\d+)(by(\\d+))?")[1,]
  if (is.na(str[1])) return(gamma.name)
  if (str[4] == "") {
    return(paste("\\Gamma_{", str[2], "}", sep=""))
  }
  return(paste("\\frac{\\Gamma_{", str[2], "}}{\\Gamma_{", str[4], "}}", sep=""))
}
alurep.gamma.texlabel = Vectorize(alurep.gamma.texlabel.nv)

##
## get tex description of a quantity in alucomb2 structure
##
alurep.get.texdescr.nv = function(descr, texdescr) {
  repeat {
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
    descr = sub("\\s+/\\s+G\\(total\\)\\s*", "", descr)
    descr = gsub("ex[.]\\s+K\\(S\\)0 --> pi- pi[+]", "ex. K0", descr, perl=TRUE)
    descr = gsub("\\s*\\(``\\d-prong''\\)", "", descr, perl=TRUE)
    descr = gsub("-->", "\\to", descr, fixed=TRUE)
    descr = gsub("([+-])", "^\\1", descr)
    descr = gsub("\\^-prong", "\\\\text{-prong}", descr)
    descr = gsub("G(\\(((?:[^()]++|(?1))+)*+\\))", "\\\\BRF{tau^-}{\\2}", descr, perl=TRUE)
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
    descr = gsub("\\s+\\(ex[.]", "\\\\;(\\\\text{ex.}", descr)
    descr = gsub("(neutrals)", "\\\\text{\\1}", descr)
    descr = gsub("(particles|strange|total)", "\\\\text{\\1}", descr)
    descr = gsub("(.*\\S)\\s*/\\s*(\\S.*)$", "\\\\frac{\\1}{\\2}", descr, perl=TRUE)
    texdescr = descr
    break
  }

  return(texdescr)
}
alurep.get.texdescr = Vectorize(alurep.get.texdescr.nv)

##
## get tex description of a quantity in alucomb2 structure
##
alurep.tex.quant.descr = function(quant) {
  alurep.get.texdescr.nv(quant$descr, quant$texdescr)
}

##--- return latex command def with specified multi-line body
alurep.tex.cmd = function(cmd, body) {
  paste("\\newcommand{\\", cmd, "}{%\n", body, "%\n}\n", sep="")
}

##--- return latex command def with specified one-line body
alurep.tex.cmd.short = function(cmd, body) {
  paste("\\newcommand{\\", cmd, "}{", body, "\\xspace}\n", sep="")
}

##
## return measurement value, stat, syst original values
## all properly formatted for latex printing, separately for
## - quant = val +- stat +- syst
## - val
## - stat
## - syst
##
alurep.tex.meas.val.card.fields = function(meas) {
  if (attr(meas$stat, "input") != "") {
    stat.txt = paste("\\pm", attr(meas$stat, "input"))
  } else {
    stat.txt = paste("{}^{", attr(meas$stat.p, "input"), "}_{", attr(meas$stat.n, "input"), "}", sep="")
  }
  if (attr(meas$syst, "input") != "") {
    syst.txt = paste("\\pm", attr(meas$syst, "input"))
  } else {
    syst.txt = paste("{}^{", attr(meas$syst.p, "input"), "}_{", attr(meas$syst.n, "input"), "}", sep="")
  }
  val.txt = attr(meas$value.orig, "input")

  order = gsub("^[^e]+(|e[+]?([-])?0*(\\d+))$", "\\2\\3", c(val.txt, stat.txt, syst.txt), perl=TRUE)
  order = ifelse(order == "", 0, as.numeric(order))

  if(max(order) == -2) stop()
  
  val.tex = gsub("e[+]?([-])?0*(\\d+)", "\\\\cdot 10^{\\1\\2}", val.txt, ignore.case=TRUE)
  stat.tex = gsub("e[+]?([-])?0*(\\d+)", "\\\\cdot 10^{\\1\\2}", stat.txt, ignore.case=TRUE)
  syst.tex = gsub("e[+]?([-])?0*(\\d+)", "\\\\cdot 10^{\\1\\2}", syst.txt, ignore.case=TRUE)
  
  if (max(order) == min(order)) {
    if (max(order) == 0) {
      if (syst.txt == "\\pm 0") {
        quant = paste(val.tex, stat.tex)
      } else {
        quant = paste(val.tex, stat.tex, syst.tex)
      }
    } else {
      if (syst.txt == "\\pm 0") {
        rc = c(val.txt, stat.txt)
      } else {
        rc = c(val.txt, stat.txt, syst.txt)
      }
      quant = paste(gsub("^([^e]+)(|e[+]?([-])?0*(\\d+))$", "\\1", rc, perl=TRUE), collapse=" ")
      quant = paste("(", quant, ") \\cdot 10^{", max(order), "}")
    }
  } else {
    if (syst.txt == "\\pm 0") {
      quant = paste(val.txt, stat.txt)
    } else {
      quant = paste(val.txt, stat.txt, syst.txt)
    }
  }
  return(list(quant=quant, val=val.txt, stat=stat.txt, syst=syst.txt))
}

##
## return value +- stat +- syst original values
## convert exponential format for latex
##
alurep.tex.meas.val.card = function(meas) {
  rc = alurep.tex.meas.val.card.fields(meas)
  return(rc$quant)
}

##
## for the spec vector of measurement, compute the
## appropriate precision and power-of-ten order
##

alurep.precision.order = function(vals, perc=FALSE, signif=4, signif.min=2) {
  vals = vals[vals != 0]
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
  return(c(precision=precision, order=order.max))
}

##
## return numeric value formatted val +- stat in a string
## use automatic precision and order of magnitude
##
alurep.tex.val.auto = function(vals, width=0, perc=FALSE) {
  rc = alurep.precision.order(vals, perc)
  precision = rc[1]
  order = rc[2]
  vals = vals / 10^order
  if (order == 0) {
    rc = sprintf(paste("%", width, ".", precision, "f", sep=""), vals)
  } else if (perc && order == -2) {
    ##--- will act depending on  alurep.precision.order()
    rc = sprintf(paste("%", width, ".", precision, "f\\%%", sep=""), vals)
  } else {
    rc = sprintf(paste("%", width, ".", precision, "f\\cdot 10^{%d}", sep=""), vals, order)
  }
  return(rc)
}

##
## return numeric value formatted val +- stat in a string
## according to the specified precision and power-of-ten order
##
alurep.tex.val.prec.ord = function(quant.val, precision, order, width=0, perc=FALSE) {
  quant.val = quant.val/10^order
  if (order == 0) {
    rc = sprintf(paste("%", width, ".", precision, "f", sep=""), quant.val)
  } else if (perc && order == -2) {
    ##--- will act depending on  alurep.precision.order()
    rc = sprintf(paste("%", width, ".", precision, "f%%", sep=""), quant.val)
  } else {
    rc = sprintf(paste("%", width, ".", precision, "f\\cdot 10^{%d}", sep=""), quant.val, order)
  }
  return(rc)
}

##
## return quantity formatted val +- stat in a string
## according to the specified precision and power-of-ten order
##
alurep.tex.val.err.prec.ord = function(quant.val, quant.err, precision, order, width=0, perc=FALSE) {
  quant.val = quant.val/10^order
  quant.err = quant.err/10^order
  if (order == 0) {
    rc = sprintf(paste("%", width, ".", precision, "f \\pm %", width, ".", precision, "f", sep=""),
      quant.val, quant.err)
  } else {
    rc = sprintf(paste("(%", width, ".", precision, "f \\pm %", width, ".", precision, "f) \\cdot 10^{%d}", sep=""),
      quant.val, quant.err, order)
  }
  return(rc)
}

##
## return quantity formatted val +- stat in a string
## according to the self-determined optimal precision and power-of-ten order
##
alurep.tex.val.err.prec.ord.auto = function(quant.val, quant.err, width=0, perc=FALSE) {
  rc = alurep.precision.order(c(quant.val, quant.err), perc)
  precision = rc[1]
  order = rc[2]
  return(alurep.tex.val.err.prec.ord(quant.val, quant.err, precision, order, width, perc))
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
  return(rc)
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
  return(alurep.precision.order(vals, perc))
}

##
## return measurement formatted val +- stat +- syst in a string
## according to the specified precision and power-of-ten order
##
alurep.tex.meas.val.auto = function(meas, width=0, perc=FALSE) {
  rc = alurep.precision.order.meas(meas, perc)
  precision = rc[1]
  order = rc[2]
  return(alurep.tex.meas.val(meas, precision, order, width, perc))
}

##
## get all measurements related to a quantity
## compute appropriate precision and order of magnitude
##
alurep.precision.order.quant = function(quant.name, perc=FALSE) {
  meas.names = alurep.meas.quant(quant.name, delta)
  vals = unlist(lapply(measurements[meas.names], function(m)
    c(m$value.orig, m$stat.p, m$stat.n, m$syst.p, m$syst.n)))
  vals = c(vals, quant.val[quant.name], quant.err[quant.name])
  return(alurep.precision.order(vals, perc))
}
