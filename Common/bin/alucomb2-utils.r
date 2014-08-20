## ////////////////////////////////////////////////////////////////////////////
##
## utility functions for alucomb2.r
##

## ////////////////////////////////////////////////////////////////////////////
## functions

##
## create diagonal matrix also for vectors of length one
##
diag.m = function(vec) {
  if (length(vec) <= 1) {
    rc = as.matrix(vec)
  } else {
    rc = diag(vec)
  }
  rc
}

##
## substitute after evaluating arg
##
esub = function(expr, sublist=NULL) do.call("substitute", list(expr, sublist))
esub.expr = function(expr, sublist=NULL) {
  sapply(as.expression(expr), function(call) as.expression(esub(call, sublist)))
}

##
## deparse expression and produce single line
##
deparse.one.line = function(expr) {
  paste(gsub("^\\s+|\\s+$", "", sapply(as.expression(expr), function(x) deparse(x)), perl=TRUE), collapse="")
}

##
## print matrix with even formatting
##
alu.rbind.print = function(x, width=13, num.columns=NULL) {
  number.width = width
  row.names = rownames(x)
  col.names = colnames(x)
  if (is.null(row.names)) row.names = rep("", nrow(x))
  if (is.null(col.names)) col.names = rep("", ncol(x))
  rn.max = max(nchar(row.names))
  cn.max = max(nchar(col.names))
  cn.max = max(number.width, cn.max)
  fmt = paste("% -", number.width, ".", number.width-7, "g", sep="")
  width = getOption("width")
  digits = getOption("digits")
  if (is.null(num.columns)) {
    items.per.row = max(1, floor((width-rn.max) / (cn.max+1)))
  } else {
    items.per.row = num.columns
  }
  for (i.first in seq(1, ncol(x), by=items.per.row)) {
    ## print(format(x), quote=FALSE)
    i.last = min(i.first + items.per.row - 1, ncol(x))
    cat(format("", width=rn.max+1), paste(format(col.names[i.first:i.last], width=cn.max), collapse=" "), "\n", sep="")
    mapply(function(label, vec) {
      cat(format(label, width=rn.max), " ", paste(format(sprintf(fmt, unlist(vec)), width=cn.max), collapse=" "), "\n", sep="")
    }, row.names, apply(x[,i.first:i.last, drop=FALSE], 1, list))
  }
  return(invisible(NULL))
}

##
## test if pattern matches string irrespective of letters case
##
match.nocase = function(pattern, str) {
  return(regexpr(pattern, str, ignore.case=TRUE) != -1)
}

##--- measurement names matching a specific name
meas.match = function(gamma) {
  return( regexpr(paste("[.]", gamma, "[.]", sep=""), meas.names) != -1 )
}

##--- sum in quadrature
quadrature = function(x) {
  return(sqrt(sum(x^2)))
}

##--- normalize publication state
alu.norm.pubstate = function(str) {
  str = sub("^pub[^.]*", "pub", str, perl=TRUE)
  str = sub("^sub[^.]*", "sub", str, perl=TRUE)
  str = sub("^acc[^.]*", "acc", str, perl=TRUE)
  str = sub("^prelim[^.]*", "prelim", str, perl=TRUE)
  str = sub("^PUB[^.]*", "PUB", str, perl=TRUE)
  str = sub("^SUB[^.]*", "SUB", str, perl=TRUE)
  str = sub("^ACC[^.]*", "ACC", str, perl=TRUE)
  str = sub("^PRELIM[^.]*", "PRELIM", str, perl=TRUE)
  return(str)
}

##
## return asymmetric excursions, possibly as percent or relative
##
alucomb2.handle.asymm.excurs = function(value, delta) {
  delta.p = as.vector(delta)
  if (is.null(attr(delta, "negval"))) {
    ##--- no separate negative value
    if (!is.null(attr(delta, "%"))) {
      delta.p = value * delta.p /100
    } else if (!is.null(attr(delta, "&"))) {
      delta.p = value * delta.p
    }
    delta.n = -delta.p
    delta.pn = delta.p
    attr(delta.p, "input") = ""
    attr(delta.n, "input") = ""
    attr(delta.pn, "input") = attr(delta, "input")
  } else {
    ##--- there is a separate negative value
    delta.n = as.vector(attr(delta, "negval"))
    if (!is.null(attr(delta, "%"))) {
      delta.p = value * delta.p /100
      delta.n = value * delta.n /100
    } else if (!is.null(attr(delta, "&"))) {
      delta.p = value * delta.p
      delta.n = value * delta.n
    }
    ##--- symmetrize according to combos prescription
    delta.pn = sqrt( (delta.p^2 + delta.n^2)/2 )
    attr(delta.p, "input") = attr(delta, "input")
    attr(delta.n, "input") = attr(attr(delta, "negval"), "input")
    attr(delta.pn, "input") = ""
  }
  return(list(delta.pn, delta.p, delta.n))
}

##--- utility function to only have an empty output line when needed
print.empty.line.first.time = function() {
  if (flag.empty.line.first.time) {
    cat("\n")
    flag.empty.line.first.time <<- FALSE
  }
}

alucomb2.def.global.flag = function() {
  flag = FALSE
  list(
       set = function(x) flag <<- x,
       get = function() flag
       )}

alucomb2.eol.first.time = list(
  flag = alucomb2.def.global.flag(),
  reset = function() alucomb2.eol.first.time$flag$set(TRUE),
  print = function() {
    if (alucomb2.eol.first.time$flag$get()) {
      alucomb2.eol.first.time$flag$set(FALSE)
      cat("\n")
    }
  }
)

##--- print single parameter
alucomb2.format.param = function(label, value) {
  input = attr(value, "input")
  rc.out = paste(format(label, width=16), format(input[1], width=12))
  if (input[2] == "" && value[3] == -value[4]) {
    input[2] = substr(input[4], 2, nchar(input[4]))
  }
  if (input[2] == "") {
    rc.out = paste(rc.out, format(paste(" ", input[3], sep=""), width=12),  " ", input[4], sep="")
  } else {
    rc.out = paste(rc.out, " +-", input[2], sep="")
  }
}

##--- print parameters
alucomb2.print.params = function(params) {
  if (length(params) > 0) {
    cat("\nPARAMETERS\n")
    mapply(function(label, value) {
      cat("  ", alucomb2.format.param(label, value), "\n", sep="")
    }, names(params), params)
  }
}

##--- print systematic terms, the original input string and the resulting absolute value
alucomb2.print.meas.syst.terms = function(syst.label, syst.terms, mask) {
  syst.terms.input = attr(syst.terms, "input")[mask]
  syst.terms = syst.terms[mask]
  if (length(syst.terms) > 0) {
    cat("\n", syst.label, "\n", sep="")
    rc= mapply(function(label, val.abs, val.input) {
      cat("  ", format(val.input, width=16), " ", format(label, width=16), sep="")
      if (val.input != val.abs) cat(" # ", val.abs, sep="")
      cat("\n")
    },
      names(syst.terms.input),
      paste(ifelse(syst.terms>=0, "+", ""), as.character(syst.terms), sep=""),
      syst.terms.input)
  }
}

##--- print correlation terms, using the cards input values
alucomb2.print.meas.correlation = function(corr.label, corr.terms) {
  if (length(corr.terms) > 0) {
    cat("\n")
    rc = mapply(function(label, value, input) {
      cat(corr.label, " ", format(input, width=12), " ", gsub(".", " ", label, fixed=TRUE), sep="")
      if (value != input) cat(" # ", value, sep="")
      cat("\n")
    },
      names(corr.terms),
      paste(ifelse(corr.terms>=0, "+", ""), as.character(corr.terms), sep=""),
      attr(corr.terms, "input"))
  }
}

##
## print a measurement
##
alucomb2.print.meas = function(meas, quantities) {
  node = quantities[[meas$quant]]$node
  cat("#\n#", meas$quant)
  if (!is.null(node) && node != "") {
    cat(" (PDG node = ", node, ")\n", sep="")
  } else {
    cat("\n")
  }
  descr = quantities[[meas$quant]]$descr
  if (!is.null(descr) && descr != "") cat("# ", descr, "\n", sep="")
  cat("#\nBEGIN MEASUREMENT ", paste(meas$tags, collapse=" "), "\n", sep="")

  cat("\n")
  cat("VALUE", attr(meas$value, "input"))

  flag.stat.mismatch = FALSE
  if (attr(meas$stat, "input") != "") {
    cat(" +-", attr(meas$stat, "input"), sep="")
    if (attr(meas$stat, "input") != meas$stat) flag.stat.mismatch = TRUE
  } else {
    cat(" ", attr(meas$stat.p, "input"), " ", attr(meas$stat.n, "input"), sep="")
    if (attr(meas$stat.p, "input") != meas$stat.p ||
        attr(meas$stat.n, "input") != meas$stat.n) flag.stat.mismatch = TRUE
  }

  flag.syst.mismatch = FALSE
  if (attr(meas$syst, "input") != "") {
    cat(" +-", attr(meas$syst, "input"), sep="")
    if (attr(meas$syst, "input") != meas$syst) flag.syst.mismatch = TRUE
  } else {
    cat(" ", attr(meas$syst.p, "input"), " ", attr(meas$syst.n, "input"), sep="")
    if (attr(meas$syst.p, "input") != meas$syst.p ||
        attr(meas$syst.n, "input") != meas$syst.n) flag.syst.mismatch = TRUE
  }

  if (flag.stat.mismatch || flag.syst.mismatch) {
    cat("  # ", meas$value)
    if (meas$stat.p == -meas$stat.n) {
      cat(" +-", meas$stat, sep="")
    } else {
      cat(" +", meas$stat.p, " -", -meas$stat.n, sep="")
    }
    if (meas$syst.p == -meas$syst.n) {
      cat(" +-", meas$syst, sep="")
    } else {
      cat(" +", meas$syst.p, " -", -meas$syst.n, sep="")
    }
  }
  cat("\n")

  alucomb2.print.meas.correlation("STAT_CORR_WITH", meas$corr.terms.stat)
  alucomb2.print.meas.correlation("ERROR_CORR_WITH", meas$corr.terms.tot)

  syst.local.mask = attr(meas$syst.terms, "type") == "l"
  syst.paper.mask = attr(meas$syst.terms, "type") == "p"
  alucomb2.print.meas.syst.terms("SYSTEMATICS", meas$syst.terms, !(syst.local.mask | syst.paper.mask))
  alucomb2.print.meas.syst.terms("SYSTPAPER", meas$syst.terms, syst.paper.mask)
  alucomb2.print.meas.syst.terms("SYSTLOCAL", meas$syst.terms, syst.local.mask)

  alucomb2.print.params(meas$params)

  cat("\n")
  cat("END\n")
}

##
## handle:
## - numbers preceded by +, -, +-
## - sequences of two numbers, +a -b
##   - replaced by a single number, equal to a, with attribute "negval" equal to b
## - numbers followed by %, &
## - labels followed by %, &
## returns:
## - vector of labels (what cannot be converted to a number)
## - vector of numbers (numbers or +a -b sequences, possibly followed by %, &)
##
alucomb2.get.labels.values = function(clause.fields) {
  ##
  ## values are numbers possibly preceded by "+", "-", "+-"
  ## sequences "+num1 -num2" are replaced by a +num1 with attribute "negval" = -num2
  ## the word +-num is replaced by +num with attribute "pm"
  ##
  
  ##--- for each words get possible +/-/+- at beginning
  clause.fields.signtag = sub("^([+]|-|[+]-|).*$", "\\1", clause.fields)
  ##--- remove "+-" sequence at begin of word, which would stop conversion to numeric
  clause.fields.convert = sub("^[+]-", "", clause.fields)
  
  ##--- for each words get possible ^/% (relative and percent value)
  clause.fields.relperc = sub("^[^&%]*(&|%|)$", "\\1", clause.fields)
  ##--- remove "&" or "%" sequence at end of wrd, which would stop conversion to numeric
  clause.fields.convert = sub("(&|%)$", "", clause.fields.convert)
  
  ##--- convert if possible to numbers, else NA
  clause.fields.val = suppressWarnings(as.numeric(clause.fields.convert))
  ##--- take note wich fields were converted to numbers
  clause.fields.numeric = which(!is.na(clause.fields.val))
  ##--- labels are all non-numeric words
  clause.labels = clause.fields.convert[is.na(clause.fields.val)]
  clause.labels.relperc = clause.fields.relperc[is.na(clause.fields.val)]
  
  ##--- indices to numeric values preceded by plus sign
  plus.ind = which(clause.fields.signtag == "+")
  plus.ind = plus.ind[!is.na(clause.fields.val[plus.ind])]
  ##--- indices to numeric values preceded by minus sign
  minus.ind = which(clause.fields.signtag == "-")
  minus.ind = minus.ind[!is.na(clause.fields.val[minus.ind])]
  
  ##--- indices to sequences +<numeric value 1> -<numeric value 2>
  plus.minus.ind = intersect(plus.ind, minus.ind-1)
  
  ##--- transform to a list in order to store per-item attributes
  clause.fields.val = mapply(function(val, pm, relperc, input) {
    if (is.na(val)) return(NA)
    if (pm == "+-") attr(val, "+-") = TRUE
    if (relperc != "") attr(val, relperc) = TRUE
    ##--- save input card string without %|&
    attr(val, "input") = input
    val
  }, clause.fields.val, clause.fields.signtag,
    clause.fields.relperc, clause.fields.convert,
    SIMPLIFY=FALSE)
  
  plus.minus.val = lapply(plus.minus.ind, function(i) {
    second.relperc = clause.fields.relperc[i+1]
    if (clause.fields.relperc[i] != second.relperc) {
      stop("+num -num sequence with different suffixes...\n  ",
           paste(clause.fields[i:(i+1)], collapse=" "))
    }
    first = clause.fields.val[[i]]
    second = clause.fields.val[[i+1]]
    if (second.relperc != "") {
      attr(second, second.relperc) = NULL
    }
    attr(first, "negval") = second
    return(first)
  })
  
  ##--- replace +<numeric value 1> -<numeric value 2> with a single number
  clause.fields.val[plus.minus.ind] = plus.minus.val
  clause.fields.val[plus.minus.ind+1] = NA
  clause.values = clause.fields.val[!is.na(clause.fields.val)]
  
  ##--- apply labels attributes to corresponding values
  if (length(clause.values) > 0) {
    clause.values = mapply(function(val, relperc) {
      if (!is.na(relperc) && relperc != "") attr(val, relperc) = TRUE
      if (!is.null(attr(val, "&")) && !is.null(attr(val, "%"))) {
        stop("both % and & suffixes specified for same value in...\n  ",
             paste(clause.fields, sep=" "))
      }
      ##--- append %|& attribute to card input value (even if it was originally on its label)
      for (tag in c("&", "%")) {
        if (!is.null(attr(val, tag))) {
          attr(val, "input") = paste(attr(val, "input"), tag, sep="")
          ##--- append to "negval" attribute as well if existing
          negval = attr(val, "negval")
          if (!is.null(negval)) {
            attr(negval, "input") = paste(attr(negval, "input"), tag, sep="")
            attr(val, "negval") = negval
          }
        }
      }
      return(val)
    }, clause.values, clause.labels.relperc[1:length(clause.values)], SIMPLIFY=FALSE)
  }
  return(list(labels=clause.labels, values=clause.values, fields.numeric=clause.fields.numeric))
}

##
## get matrix deriv(meas_i, syst.term_k) * delta(syst.term_k)
##
alucomb2.meas.by.syst.term = function(measurements, syst.terms) {
  rc = sapply(measurements, function(m) {
    if (is.null(m$syst.terms)) {
      st = rep(0, length(syst.terms))
    } else {
      st = m$syst.terms[syst.terms]
      st = ifelse(is.na(st), 0, st)
    }
    names(st) = names(syst.terms)
    return(st)
  })
  return(t(rc))
}

##////////////////////////////////
##
## cards documentation
##

##
## measurement block:
##
## old style
##   BEGIN <exp> <quantity> <pub|prelim> <reference>.[<other.tags>]
##   MEASUREMENT    m_<quantity>  statistical  systematic
##   DATA           m_<quantity>  statistical  systematic
##                  <val>         <stat>       <syst>
##   ...
##   STAT_CORR_WITH 0.235119 BaBar Gamma9by5 pub AUBERT.10P
##   ...
##   ##--- systematics terms
##   DATA
##     term1 <signed val>
##     term2 <signed val>
##   END
## new style
##   BEGIN MEASUREMENT <exp> <quantity> <pub|prelim> <reference> [<other tags>]
##   VALUE <val> <stat> <syst>
##   ...
##   SYSTEMATICS
##     term1 <signed val>
##     term2 <signed val>
##   ...
##   END
##
## a "measurement" is stored in a list with fields
##
## - tags = character array: measurement name fields
## - quant = character: measured quantity
## - value = numeric,
## - stat = numeric,
## - syst = numeric,
## - params = list of numeric arrays: parameter value, pos/neg uncertainty
## - syst.terms = numeric array: syst. contributions related to external parameters
## - corr.terms.stat = numeric array: stat. corr. with other measurements
## - corr.terms.tot = numeric array: error corr. with other measurements
##

##
## combination block:
##
## old style
##   BEGIN <tag> [<tag 2> ...]
##   COMBINE * * *
##   PARAMETERS !      assumed_value  pos_excurs   neg_excurs
##     par1            0.919          +0.003       -0.003
##   ...
##   MEASUREMENT m_<quantity> statistical systematic [seed <value>] [scale <value>]
##   ...
##   END
## new style
##   BEGIN COMBINATION <tag> [<tag 2> ...]
##   PARAMETERS !      assumed_value  pos_excurs   neg_excurs
##     par1            0.919          +0.003       -0.003
##   ...
##   QUANTITY <quantity> [seed <value>] [scale <value>]
##   ...
##   END
##
## a "combination" is stored as a list with fields:
##
## - tags = character array: combination name fields
## - quantities = character array: measured quantities to be fitted
## - params = list of numeric arrays: parameter value, pos/neg uncertainty
## - constr.lin.comb = list of numeric arrays: coeffs of quantities in linear constraints
## - constr.lin.val = list or numeric: values of linear constraint equations
## - constr.nl.str.expr = list of character: non-linear constraint expression
## - constr.nl.str.val = list or character: expression of non-linear constraint
##

## ////////////////////////////////////////////////////////////////////////////
## code

##
## alucomb.read
##

alucomb.read = function(file = "") {
  if (match.nocase("^\\s*$", file)) {
    stop("alucomb: please provide as argument the card input file\n")
  }

  flag.replace.corr = FALSE
  flag.build.delta = FALSE

  if (!file.exists(file)) {
    stop("cannot find file ", file, "\n")
  }
  dir.base = dirname(file)
  lines = readLines(file)
  cat("read file", file, "\n")

  ##
  ## deal with INCLUDE lines
  ##
  iline = 1
  lines.len = length(lines)
  repeat {
    if (iline > lines.len) break
    fields = unlist(strsplit(lines[iline], "\\s+", perl=TRUE))
    if (length(fields) > 0 &&
        match.nocase("^INCLUDE$", fields[1])) {
      file.inc = paste(dir.base, fields[2], sep="/")
      if (!file.exists(file.inc)) {
        stop("cannot find included file ", file.inc, "\n")
      }
      lines.inc = readLines(file.inc)
      cat("read file", file.inc, "\n")
      lines = c(
        if (iline-1>=1) {lines[1:(iline-1)]} else {NULL},
        lines.inc,
        if (iline+1<=lines.len) {lines[(iline+1):(lines.len)]} else {NULL}
        )
      lines.len = length(lines)
      next
    }
    iline = iline + 1
  }

  ##
  ## init global storage for cards
  ##
  measurements = list()
  combination = list()

  s.inblock = FALSE
  s.inclause = FALSE

  t.inblock = FALSE
  t.inclause = FALSE
  t.continue = FALSE
  t.endclause = FALSE
  t.endblock = FALSE
  t.endfile = FALSE

  measurement.in.data = FALSE

  iline = 1
  lines.len = length(lines)
  repeat {
    ##
    ## read new line
    ##
    if (iline > lines.len) {
      t.endfile = TRUE
    } else {
      line = lines[iline]
      iline = iline+1

      ##--- skip comment and empty lines
      if (regexpr("^\\s*$", line, perl=TRUE) != -1 ||
          regexpr("^\\s*[\\#!*;]", line, perl=TRUE) != -1) {
        next
      }
      ##--- remove trailing whitespace
      line = sub("\\s+$", "", line, perl=TRUE)
      ##--- remove end-of-line comments
      line = sub("\\s*[\\#!].*$", "", line, perl=TRUE)

      ##
      ## tokenize input
      ##

      ##
      ## split line in quoted strings segments
      ## matched substrings will be
      ## - a string for any unquoted sequence between line-begin, line-end or double-quotes
      ## - a string for every quoted string
      ## - an empty string both before and after any quoted string
      ##
      matches = gregexpr("([^\"]*(\\.[^\"]*)*)", line)[[1]]
      if (length(matches) > 1) {
        ##--- special treatment for lines with quoted strings
        matches.len = attr(matches, "match.length")
        matches.zerolen = which(matches.len == 0)
        matches.qstr = matches.zerolen[c(TRUE, FALSE)] + 1
        matches.str = setdiff(seq(1, length(matches)), c(matches.zerolen, matches.qstr))
        matches.all.ind = sort(c(matches.str, matches.qstr))
        matches.qstr.flag = matches.all.ind %in% matches.qstr
        matches.all = matches[matches.all.ind]
        attr(matches.all, "match.length") = matches.len[matches.all.ind]
        line.qstrings = regmatches(line, list(matches.all))[[1]]
        ##--- split
        fields = unlist(mapply(function(word, qstr.flag) {
          if (qstr.flag) {return(word)}
          ##--- trim and split unquoted strings
          word = sub("^\\s+|\\s+$", "", word)
          unlist(strsplit(word, "\\s+", perl=TRUE))
	}, line.qstrings, matches.qstr.flag, USE.NAMES=FALSE))
        ##--- restore empty string that strsplit would return for lines that begin with whitespace
        if (regexpr("^\\s+", line, perl=TRUE) != -1) {
          fields = c("", fields)
        }
      } else {
        fields = unlist(strsplit(line, "\\s+", perl=TRUE))
      }
      ##--- block begin
      t.inblock = (fields[1] == "BEGIN")
      ##--- block end
      t.endblock = (fields[1] == "END")
      ##--- a line beginning with space is the continuation of the last clause
      t.continue = (fields[1] == "")
      ##--- begin a clause if: 1) does not begin/end a block and 2) does not begin with blank
      t.inclause = !(t.inblock || t.endblock || t.continue)
    }

    ##--- when continuing a clause collect all its fields
    if (t.continue) {
      if (s.inclause) {
        ##--- continuation of previous clause
        clause.fields = c(clause.fields, fields[-1])
      } else {
        t.continue = FALSE
        cat("warning, continuation line outside of clause, ignored, line ...\n")
        cat("  '", line, "'\n", sep="")
      }
      next
    }

    ##
    ## end of a clause -> interpret its fields
    ##
    if (s.inclause && !t.continue && !is.null(clause.fields)) {
      ##--- if no continuation then the current clause ends
      s.inclause = FALSE
      if (t.inblock) {
        cat("warning, clause ended by BEGIN block, line ...\n")
        cat("  '", line, "'\n", sep="")
      }
      if (t.endfile) {
        cat("warning, clause ended by end of file, clause...\n")
        cat("  ", paste(clause.fields, collapse=" "), "\n", sep="")
      }
      ##
      ## handle a complete clause here
      ##
      ## cat("END clause", unlist(clause.fields), "\n\n")
      clause.keyw = toupper(clause.fields[1])
      clause.fields = clause.fields[-1]
      clause.line = paste(c(clause.keyw, clause.fields), collapse=" ")
      
      if (!(clause.keyw %in% c("QUANTITY", "STAT_CORR_WITH", "ERROR_CORR_WITH"))) {
        ##--- separate numeric values from string identifiers
        clause = alucomb2.get.labels.values(clause.fields)
      }
      
      if (clause.keyw == "COMBINE") {
        ##
        ## COMBINE
        ##

        ##+++ combos
        block.type = "COMBINATION"
        ##+++ combos
        clause$labels = clause$labels[clause$labels != "*"]
        ##--- collect quantities to combine
        block$combine = c(block$combine, clause$labels)

      } else if (clause.keyw == "PARAMETERS") {
        ##
        ## PARAMETERS (same treatment in MEASUREMENT and COMBINATION block)
        ##

        ##--- after collapsing +devp -devm into devp with attr negval expect 2x values as labels
        if (length(clause$values) != 2*length(clause$labels)) {
          stop("wrong parameter data in line...\n  ", clause.line)
        }

        param.values = clause$values[seq(1, 2*length(clause$labels), by=2)]
        names(param.values) = clause$labels
        param.deltas = clause$values[seq(2, 2*length(clause$labels), by=2)]

        ##--- add defined parameters, each is an array with value and excursions
        params = mapply(
          function(value, delta) {
            rc = alucomb2.handle.asymm.excurs(value, delta)
            param.list = list(value=value, delta=rc[[1]], delta.p=rc[[2]], delta.n=rc[[3]])
            param = unlist(param.list)
            attr(param, "input") = sapply(param.list, function(el) attr(el, "input"))
            return(param)
          }, param.values, param.deltas, SIMPLIFY=FALSE)
        
        ##--- get which parameters where already input
        matched = clause$labels %in% names(block$params)

        ##--- warn about overridden parameters
        rc = mapply(
          function(label, old, new) {
            cat("\nwarning, override PARAMETER\n")
            cat("old: ", alucomb2.format.param(label, old), "\n", sep="")
            cat("new: ", alucomb2.format.param(label, new), "\n", sep="")
          }, clause$labels[matched], params[matched], block$params[clause$labels[matched]])

        ##--- override matched parameters and add the unmatched ones
        block$params[clause$labels] = params
        
      } else if (clause.keyw == "QUANTITY" || (block.type == "COMBINATION" && clause.keyw == "MEASUREMENT")) {
        ##
        ## QUANTITY or (+++combos) MEASUREMENT in COMBINE block
        ##

        ##--- get first label, quantity name
        ##+++ combos, remove leading "m_"
        meas.name = sub("^m_", "", clause.fields[1], ignore.case=TRUE)
        if (is.null(block$quantities[[meas.name]])) {
          block$quantities[[meas.name]] = list()
        }
        clause.fields = clause.fields[-1]

        ##+++ combos, remove labels named "stat*", "syst*"
        clause.fields = clause.fields[!match.nocase("^(syst|stat)", clause.fields)]

        ##--- if "descr" or "texdescr", use as argument all what follows
        descr.pos = match(c("descr", "texdescr"), tolower(clause.fields))
        descr.pos = descr.pos[which(!is.na(descr.pos))[1]]
        if (!is.na(descr.pos)) {
          descr.label = tolower(clause.fields[descr.pos])
          descr.text = paste(tail(clause.fields, -descr.pos), collapse=" ")
          ##--- remove possible double quotes
          descr.text = sub("^\"", "", sub("\"$", "", descr.text))
          clause.fields = head(clause.fields, -(descr.pos-1))
        }

        ##
        ## get label/values where values are strings
        ## for QUANTITY it is required that each value follows its label
        ##
        str.labels.mask = tolower(clause.fields) %in% c("node")
        if (length(str.labels.mask) > 0 && tail(str.labels.mask, 1)) {
          stop("label without a value in line...\n  ", clause.line)
        }
        str.values.mask = c(FALSE, head(str.labels.mask, -1))
        clause.labels = clause.fields[str.labels.mask]
        clause.values = as.list(clause.fields[str.values.mask])
        clause.fields = clause.fields[!(str.labels.mask | str.values.mask)]

        ##--- deal with remaining keywords that require numeric values
        clause = alucomb2.get.labels.values(clause.fields)
        
        ##--- add string-valued keywords
        clause$labels = c(clause$labels, clause.labels)
        clause$values = c(clause$values, clause.values)

        ##--- add description if present
        if (!is.na(descr.pos)) {
          clause$labels = c(clause$labels, descr.label)
          clause$values = c(clause$values, descr.text)
        }

        if (length(clause$labels) != length(clause$values)) {
          stop("mismatch between labels and numeric values in line...\n  ", clause.line)
        }

        if (length(clause$values) > 0) {
          labels.exist = names(block$quantities[[meas.name]])
          labels.override = clause$labels[clause$labels %in% labels.exist]
          if (length(labels.override) > 0) {
            cat("\nwarning, override QUANTITY", meas.name, "\n")
            rc = print(rbind(
              unlist(block$quantities[[meas.name]][labels.override]),
              unlist(clause$values[clause$labels %in% labels.override])
              ))
          }
          block$quantities[[meas.name]][clause$labels] = unlist(clause$values)
        }

      } else if (clause.keyw == "MEASUREMENT") {
        ##
        ##+++ MEASUREMENT clause, combos compatibility
        ##
        block.type = "MEASUREMENT"
        meas.name = sub("^m_", "", clause.fields[1], ignore.case=TRUE)
        if (meas.name != block.fields[2]) {
          stop("MEASUREMENT quantity does not match BLOCK quantity...\n  ", clause.line)
        }
        measurement.in.data = TRUE
      } else if (clause.keyw == "VALUE" || (clause.keyw == "DATA" && measurement.in.data)) {
        ##
        ## VALUE --> enter measurement value
        ##+++ DATA after MEASUREMENT clause, combos compatibility
        ##

        if (length(clause$values) == 3) {
          if (!is.null(attr(clause$values[[1]], "negval"))) {
            stop("\"+val -val\" used for measurements value in line...\n  ", clause.line)
          }
          ##--- use values if existing
          block.meas$value = clause$values[[1]]

          ##--- get possibly asymmetric stat error, possibly percent or relative
          rc = alucomb2.handle.asymm.excurs(clause$values[[1]], clause$values[[2]])
          block.meas$stat = rc[[1]]
          block.meas$stat.p = rc[[2]]
          block.meas$stat.n = rc[[3]]

          ##--- get possibly asymmetric syst error, possibly percent or relative
          rc = alucomb2.handle.asymm.excurs(clause$values[[1]], clause$values[[3]])
          block.meas$syst = rc[[1]]
          block.meas$syst.p = rc[[2]]
          block.meas$syst.n = rc[[3]]

          ##+++ set end of MEASEREMENT DATA (combos compatibility)
          measurement.in.data = FALSE
        } else {
          stop("wrong number of numeric data in line...\n  ", clause.line)
        }
      } else if (clause.keyw == "DATA" ||
                 clause.keyw == "SYSTEMATICS" || clause.keyw == "SYSTLOCAL" || clause.keyw == "SYSTPAPER" ) {
        ##
        ## SYSTEMATICS or (+++combos) DATA
        ## list systematic terms with a label that possibly refers
        ## to a global parameter or to the same effect in another
        ## measurement
        ##
        ## SYSTLOCAL like SYSTEMATICS but add measurement name to make it local
        ## SYSTPAPER like SYSTEMATICS but add reference name to make it paper-wide
        ##

        if (length(clause$values) < 1) {
          stop("at least one numeric value required in line...\n  ", clause.line)
        }

        ##
        ## - "%" on either name or value means percent of measurement value
        ## - "&" on either name or value means relative of measurement value
        ##
        clause$values = lapply(clause$values, function(val) {
          if (!is.null(attr(val, "&"))) {
            val = val*block.meas$value
          } else if (!is.null(attr(val, "%"))) {
            val = val*block.meas$value/100
          }
          return(val)
        })

        clause$labels.input = clause$labels
        type.attr = "g"
        if (clause.keyw == "SYSTLOCAL") {
          type.attr = "l"
          clause$labels = paste(paste(block.fields, collapse="."), clause$labels, sep=".")
        } else if (clause.keyw == "SYSTPAPER") {
          type.attr = "p"
          clause$labels = paste(paste(block.fields[-2], collapse="."), clause$labels, sep=".")
        }
        if (any(clause$labels %in% names(block.meas$syst.terms))) {
          stop("already used systematics label in line...\n  ", clause.line)
        }
        names(clause$values) = clause$labels

        input.attr = sapply(clause$values, function(el) attr(el, "input"))
        input.attr = sub("^[+-]", "", input.attr)
        input.attr = paste(ifelse(clause$values >= 0, "+", "-"), input.attr, sep="")
        names(input.attr) = clause$labels.input

        type.attr = sapply(clause$values, function(el) type.attr)

        input.attr = c(attr(block.meas$syst.terms, "input"), input.attr)
        type.attr = c(attr(block.meas$syst.terms, "type"), type.attr)
        block.meas$syst.terms = c(block.meas$syst.terms, unlist(clause$values))
        attr(block.meas$syst.terms, "input") = input.attr
        attr(block.meas$syst.terms, "type") = type.attr

      } else if (clause.keyw == "STAT_CORR_WITH" || clause.keyw == "ERROR_CORR_WITH") {
        ##
        ## STAT_CORR_WITH, ERROR_CORR_WITH
        ##

        clause = alucomb2.get.labels.values(clause.fields[1])
        if (length(clause$values) == 1) {
          ##--- alucomb2.r format (has a number in 1st word)
          corr = clause$values[[1]]
          clause$labels = clause.fields[-1]
        } else {
          ##--- if first word is not a number, check for old format, 3 strings followed by a number
          corr = NA
          if (length(clause.fields) == 4) {
            clause = alucomb2.get.labels.values(clause.fields[4])
            if (length(clause$values) == 1) {
              corr = clause$values[[1]]
              clause$labels = clause.fields[-4]
            }
          }
          if (length(clause.fields) != 4 || is.na(corr)) {
            stop("invalid ", clause.keyw, " clause in line...  ...\n  ", clause.line)
          }
        }

        if (!is.null(attr(corr, "+-"))) {
          stop("invalid ", clause.keyw, " clause with \"+-\" in line...  ...\n  ", clause.line)
        }
          
        ##--- fix the format of pub status to pub|prelim
        clause$labels[3] = alu.norm.pubstate(clause$labels[3])
        ##--- set name to correlated measurement
        names(corr) =  paste(clause$labels, collapse=".")
        
        if (!is.null(attr(corr, "%"))) {
          ##--- percent correlation coefficient
          corr = corr/100
        }
        
        ##--- add + if not present and positive value
        attr(corr, "input") = sub("^([^+-])", "+\\1", attr(corr, "input"))

        ##--- update meas block
        input.attr = attr(corr, "input")
        names(input.attr) = names(corr)
        if (clause.keyw == "STAT_CORR_WITH") {
          if (any(names(corr) %in% names(block.meas$corr.terms.stat))) {
            stop("correlation specified twice in line...\n  ", clause.line)
          }
          input.attr = c(attr(block.meas$corr.terms.stat, "input"), input.attr)
          block.meas$corr.terms.stat = c(block.meas$corr.terms.stat, corr)
          attr(block.meas$corr.terms.stat, "input") = input.attr
        } else {
          if (any(names(corr) %in% names(block.meas$corr.terms.tot))) {
            stop("correlation specified twice in line...\n  ", clause.line)
          }
          input.attr = c(attr(block.meas$corr.terms.tot, "input"), input.attr)
          block.meas$corr.terms.tot = c(block.meas$corr.terms.tot, corr)
          attr(block.meas$corr.terms.tot, "input") = input.attr
        }

      } else if (clause.keyw == "SUMOFQUANT" || clause.keyw == "COMBOFQUANT" || clause.keyw == "CONSTRAINT") {
        ##
        ## SUMOFQUANT, COMBOFQUANT, CONSTRAINT
        ##

        if (clause.keyw != "CONSTRAINT") {
          constr.name = paste(clause.fields[1], "coq", sep=".")
          constr.val = 0
        }

        if (clause.keyw == "SUMOFQUANT") {
          constr.comb = c(-1, rep(1, length(clause.fields)-1))
          names(constr.comb) = clause.fields
        } else if (clause.keyw == "COMBOFQUANT") {
          if (length(clause$labels)-1 != length(clause$values)) {
            stop("error: mismatch between labels and data values in line...\n  ", clause.line)
          }
          constr.comb = c(-1, unlist(clause$values))
          names(constr.comb) = clause$labels
        } else if (clause.keyw == "CONSTRAINT") {
          if (length(clause$labels) != length(clause$values)) {
            stop("error: mismatch between labels and data values in line...\n  ", clause.line)
          }
          constr.name = clause$labels[1]
          constr.val = clause$values[[1]]
          constr.comb = unlist(clause$values[-1])
          names(constr.comb) = clause$labels[-1]
        }

        block$constr.lin.val[[constr.name]] = as.vector(constr.val)
        block$constr.lin.comb[[constr.name]] = constr.comb

      } else if (clause.keyw == "NLCONSTRAINT") {
        ##
        ## NLCONSTRAINT = non-linear constraint equation
        ##

        if (is.na(suppressWarnings(as.numeric(clause.fields[2])))) {
          stop("error: missing numeric value in non-linear constraint in...\n  ", clause.line)
        }
        if (FALSE && length(clause$labels) != 2) { #+++ fields preprocessing removed, disable check
          stop("error: NLCONSTRAINT needs one numeric value and one expression in line...\n  ", clause.line)
        }
        constr.name = clause$labels[1]
        constr.val = as.vector(clause$values[[1]])
        constr.expr = paste(clause$labels[-1], collapse=" ")
        rc = suppressWarnings(try(parse(text=constr.expr), silent=TRUE))
        if (inherits(rc, "try-error")) {
          stop("error, malformed NLCONSTRAINT\n", constr.name, " ", constr.val, " \"", constr.expr, "\"")
        }
        ##--- warn about overridden NLCONSTRAINT
        if (!is.null(block$constr.nl.str.val[[constr.name]])) {
          cat("\nwarning, override NLCONSTRAINT\n")
          old.str =
            paste("old: ", block$constr.nl.str.val[[constr.name]], " = ", block$constr.nl.str.expr[[constr.name]], "\n", sep="")
          new.str = paste("new: ", constr.val, " = ", constr.expr, "\n", sep="")
          if (nchar(old.str) > 78 || nchar(new.str) > 78) {
            old.str = paste("\n", old.str, sep="")
            new.str = paste("\n", new.str, sep="")
          }
          cat(old.str)
          cat(new.str)
        }
        block$constr.nl.str.val[[constr.name]] = constr.val
        block$constr.nl.str.expr[[constr.name]] = constr.expr

      } else if (clause.keyw == "MODMEAS") {
        ##
        ## MODMEAS drop|keep <measurement tags>
        ##
        keyw = toupper(clause.fields[1])
        if (keyw == "KEEP" || keyw == "DROP") {
          ##--- fix the format of pub status to pub|prelim
          meas.name = clause.fields[-1]
          meas.name[3] = alu.norm.pubstate(meas.name[3])
          meas.name = paste(meas.name, collapse=".")
          if (keyw == "DROP") {
            if (!meas.name %in% block$meas.drop.cards) {
              block$meas.drop.cards = c(block$meas.drop.cards, meas.name)
            }
          } else {
            block$meas.drop.cards = setdiff(block$meas.drop.cards, meas.name)
          }
        } else {
          stop("error, invalid MODMEAS keyword in line...\n  ", clause.line)
        }

      } else {
        stop("error, invalid keyword in line...\n  ", clause.line)
      }
    }

    if (t.inclause) {
      ##--- new clause (previous one has been handled already)
      t.inclause = FALSE
      if (s.inblock) {
        s.inclause = TRUE
        clause.fields = fields
      } else {
        cat("warning, begin clause outside of block, ignored, line ...\n")
        cat("  '", line, "'\n", sep="")
      }
    }

    if (t.inblock) {
      t.inblock = FALSE
      if (s.inblock) {
        cat("error, nested BEGIN block, ignored, line ...\n")
        cat("  '", line, "'\n", sep="")
        next
      }
      s.inblock = TRUE
      clause.fields = NULL

      ##--- remove "BEGIN"
      block.fields = fields[-1]

      block.type = ""
      if (length(block.fields) > 0) {
        block.type = toupper(block.fields[1])
        if (block.type == "COMBINATION" || block.type == "MEASUREMENT") {
          block.fields = block.fields[-1]
        } else {
          ##+++ temporarily just warning for COMBOS-like cards compatibility
          if (TRUE) {
            cat("warning, unknown BEGIN block in line...\n")
            cat(" ", fields, "\n")
          }
        }
      }
      ##
      ## initialization both for MEASUREMENT and COMBINATION
      ## when combos compatibility is abandoned, one can init just the relevant list
      ## block$params is used for parameters both of MEASUREMENT and COMBINATION
      ##

      ##--- special treatment in case this is a MEASUREMENT block
      meas.fields = block.fields
      ##+++combos: separate tags joint with dots
      if (length(meas.fields) == 3) {
        meas.fields = c(meas.fields[-3], strsplit(meas.fields[3], "[.]", perl=TRUE)[[1]])
      }
      meas.fields[3] = alu.norm.pubstate(meas.fields[3])
      meas.tag = paste(meas.fields, collapse=".")
      rm(meas.fields)

      block.meas = list()

      block = list()
      block$params = list()
      block$quantities = list()
      block$constr.lin.val = list()
      block$constr.lin.comb = list()
      block$constr.nl.str.val = list()
      block$constr.nl.str.expr = list()
    }

    if (s.inblock && (t.endblock || t.endfile)) {
      ##--- end block
      t.endblock = FALSE
      s.inblock = FALSE
      ## cat("END block", unlist(block.fields), "\n\n")
      if (block.type == "MEASUREMENT") {
        ##
        ## MEASUREMENT BLOCK
        ##

        if (length(block.fields) < 3) {
          stop("too few fields in line...\n  ", paste(c("BEGIN", block.type, block.fields), collapse=" "))
        }

        block.meas$quant = block.fields[2]
        block.meas$tags = strsplit(meas.tag, "[.]", perl=TRUE)[[1]]
        block.meas$params = block$params

        if (!is.null(measurements[[meas.tag]])) {
          old.meas = measurements[[meas.tag]]
          cat("warning, BEGIN update measurement", meas.tag, "\n")
          el.check = intersect(names(old.meas), names(block.meas))
          for (el in el.check) {
            if (any(unlist(old.meas[[el]]) != unlist(block.meas[[el]]))) {
              old.val = unlist(old.meas[[el]])
              new.val = unlist(block.meas[[el]])
              if (is.null(names(old.val))) {
                names(old.val) = el
                names(new.val) = el
              } else {
                cat(el, "\n")
                val.names = unique(c(names(old.val), names(new.val)))
                old.val = old.val[val.names]
                new.val = new.val[val.names]
              }
              alu.rbind.print(rbind(old=old.val, new=new.val))
            } else if (!identical(old.meas[[el]], block.meas[[el]])) {
              cat("\n<<<<<<<< old", el, "\n")
              print(old.meas[[el]])
              cat(">>>>>>>> new", el, "\n")
              print(block.meas[[el]])
            }
          }
          el.removed = setdiff(names(old.meas), names(block.meas))
          if (length(el.removed) > 0) {
            cat("<<<<<<<< removed elements\n")
            for (el in el.removed) {
              print(old.meas[[el]])
            }
          }
          el.added = setdiff(names(block.meas), names(old.meas))
          if (length(el.added) > 0) {
            cat(">>>>>>>> added elements\n")
            for (el in el.added) {
              print(block.meas[[el]])
            }
          }
          cat("warning, END update measurement", meas.tag, "\n")
        }

        ##
        ## measured quantity = 2nd tag in the BEGIN MEASUREMENT statement
        ## BEGIN [MEASUREMENT] <esperiment> <quantity> <pub|prelim> <reference> [...]
        ##
        measurements[[meas.tag]] = block.meas

      } else if (block.type == "COMBINATION") {
        ##
        ## COMBINE BLOCK
        ##
        block$tags = unlist(block.fields)
        if (is.null(block$combine) || length(block$combine) == 0) {
          ##--- average all quantities mentioned in block if no explicit selection
          block$combine = names(block$quantities)
        }

        ##--- add "USE = 1" quantities, drop "USE = 0" quantities
        quant.use = unlist(lapply(block$quantities, function(el) { unname(el["use"]) }))
        quant.drop = names(quant.use[quant.use == 0])
        quant.use = names(quant.use[quant.use != 0])
        block$combine = setdiff(block$combine, quant.drop)
        block$combine = c(block$combine, setdiff(quant.use, block$combine))

        ##--- insure the descr field has at least an empty string
        quant.descr = names(unlist(lapply(block$quantities, function(el) { unname(el["descr"]) })))
        quant.descr.wo = setdiff(names(block$quantities), quant.descr)
        block$quantities = lapply(block$quantities,
          function(quant) {
            if (is.null(quant$descr)) quant$descr = ""
            quant
          })

        combination = block

      } else {
        stop("unknown block type in line...\n  ", paste(c("BEGIN", block.type, block.fields[-1]), collapse=" "))
      }
    }

    if (t.endblock && !s.inblock) {
      t.endblock = FALSE
      cat("error, END block without BEGIN, ignored, line ...\n")
      cat("  '", line, "'\n", sep="")
    }

    if (t.endfile) {
      t.endfile = FALSE
      break
    }
  }
  invisible(list(combination=combination, measurements=measurements))
} ##--- end of alucomb.read

##
## transform the name of a measurement in the HFAG-tau format into a format usable by Root
##
hfag.to.root = function(str) {
  if (str == "TauMass") {
    return("m_{#tau} [MeV]")
  }

  str.orig = str

  str = gsub("Numbar", "#bar{#nu}_{#mu}",str)

  str = gsub("(^|[^A-Z])([A-Z][a-y]*)([z+-]*)b", "\\1#bar{\\2}\\3", str, perl=TRUE)
  str = gsub("F1", "f_{1}",str)
  str = gsub("Pi", "#pi",str)

  str = gsub("Nue", "#nu_{e}",str)
  str = gsub("Num", "#nu_{#mu}",str)
  str = gsub("Nut", "#nu_{#tau}",str)
  str = gsub("Nu", "#nu",str)

  str = gsub("M", "#mu",str)
  str = gsub("H", "h",str)
  str = gsub("m($|[#}A-Zh])", "^{-}\\1", str, perl=TRUE)
  str = gsub("p($|[#}A-Zh])", "^{+}\\1", str, perl=TRUE)
  str = gsub("z($|[#}A-Zh])", "^{0}\\1", str, perl=TRUE)

  if (str.orig %in% c("HmHmHpNu", "PimKmPipNu", "PimPimPipNu")) {
    str = paste(str, "(ex.K^{0})")
  }

  str = paste("B(#tau^{-} #rightarrow ", str, ")", sep="")
  return(str)
}

##
## compute root of symmetric real positive semi-definite matrix
##
alu.matr.sqrt.symm.semipos = function(X, tol = sqrt(.Machine$double.eps)) {
  if (!is.matrix(X))
    X <- as.matrix(X)
  if (length(dim(X)) > 2L || !(is.numeric(X)))
    stop("'X' must be a numeric or complex matrix")

  X = (X + t(X))/2
  rc.e = eigen(X)
  comp.positive = rc.e$values > tol*rc.e$values[1]
  rc.vals = ifelse(comp.positive, rc.e$values, 0)
  rc = rc.e$vectors %*% diag.m(sqrt(rc.vals)) %*% t(rc.e$vectors)
  rc = (rc + t(rc))/2
  rownames(rc) = rownames(X)
  colnames(rc) = colnames(X)

  return(rc)
}

##
## compute pseudo-inverse square root of symmetric real positive semi-definite matrix
## zero eigenvalues produce zero eigenvalues in the inverse-square matrix
##
alu.matr.inv.sqrt.symm.semipos = function(X, tol = sqrt(.Machine$double.eps)) {
  if (!is.matrix(X))
    X <- as.matrix(X)
  if (length(dim(X)) > 2L || !(is.numeric(X)))
    stop("'X' must be a numeric or complex matrix")

  X = (X + t(X))/2
  rc.e = eigen(X)
  ##--- force to zero computationally zero or negative eigenvalues
  comp.positive = rc.e$values > tol*rc.e$values[1]
  rc.vals = ifelse(comp.positive, rc.e$values, 0)
  rc.vecs.inv = ifelse(rc.vals > 0, 1/rc.vals, 0)

  rc.inv = rc.e$vectors %*% diag.m(rc.vecs.inv) %*% t(rc.e$vectors)
  rc.inv = (rc.inv + t(rc.inv))/2
  rownames(rc.inv) = rownames(X)
  colnames(rc.inv) = colnames(X)
  rc = rc.e$vectors %*% diag.m(sqrt(rc.vecs.inv)) %*% t(rc.e$vectors)
  rc = (rc + t(rc))/2
  rownames(rc) = rownames(X)
  colnames(rc) = colnames(X)

  attr(rc, "pos.eigen.num") = sum(comp.positive)
  attr(rc, "inv") = rc.inv

  attr(rc, "dof.eff") = apply(rc.e$vectors, 1, function(x) sum(x^2 * (rc.vals != 0)))
  ## attr(rc, "dof.factors") = ifelse(attr(rc, "dof.factors") != 0, 1/sqrt(attr(rc, "dof.factors")), 0)
  names(attr(rc, "dof.eff")) = rownames(X)

  return(rc)
}

##
## compute pseudo-inverse square root of symmetric real positive semi-definite matrix
## zero eigenvalues produce zero eigenvalues in the inverse-square matrix
## apply provided weight vector to matrix before inversion in order to reduce computational errors
## the resulting matrix Y satisfies: t(Y) %*% Y = X^-1
##
alu.matr.inv.sqrt.symm.semipos.norm = function(X, X.norm, tol = sqrt(.Machine$double.eps)) {
  if (!is.matrix(X))
    X <- as.matrix(X)
  if (length(dim(X)) > 2L || !is.numeric(X))
    stop("'X' must be a numeric matrix")

  if (dim(X)[1] != dim(X)[2])
    stop("'X' must be a square matrix")

  if (!is.numeric(X.norm))
    stop("'X.norm' must be a numeric vector/matrix")

  if (length(X.norm) != dim(X)[1])
    stop("X.norm vector dimension must match X matrix rank")

  XX = diag.m(1/X.norm) %*% X %*% diag.m(1/X.norm)
  rownames(XX) = rownames(X)
  colnames(XX) = colnames(X)
  XX.inv.sqrt.attr = alu.matr.inv.sqrt.symm.semipos(XX, tol)
  XX.inv.sqrt = XX.inv.sqrt.attr %*% diag.m(1/X.norm)
  rownames(XX.inv.sqrt) = rownames(X)
  colnames(XX.inv.sqrt) = colnames(X)

  attr(XX.inv.sqrt, "inv") = diag.m(1/X.norm) %*% attr(XX.inv.sqrt.attr, "inv") %*% diag.m(1/X.norm)
  attr(XX.inv.sqrt, "inv") = (attr(XX.inv.sqrt, "inv")+ t(attr(XX.inv.sqrt, "inv")))/2
  rownames(attr(XX.inv.sqrt, "inv")) = rownames(X)
  colnames(attr(XX.inv.sqrt, "inv")) = colnames(X)

  attr(XX.inv.sqrt, "pos.eigen.num") = attr(XX.inv.sqrt.attr, "pos.eigen.num")
  attr(XX.inv.sqrt, "dof.eff") = attr(XX.inv.sqrt.attr, "dof.eff")

  return(XX.inv.sqrt)
}
