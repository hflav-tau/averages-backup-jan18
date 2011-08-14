#!/usr/bin/env Rscript

## ////////////////////////////////////////////////////////////////////////////
## definitions

##
## create diagonal matrix also for vectors of length one
##
diag.m <- function(vec) {
  if (length(vec) <= 1) {
    rc = diag(as.matrix(vec))
  } else {
    rc = diag(vec)
  }
  rc
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

##
## representation of Combos input card files data
##

##
## a "measurement" is stored in a list with fields
## - tag = "character",
## - value = "numeric",
## - stat = "numeric",
## - syst = "numeric",
## - bibitem = "character",
## - params = "list",
## - syst.terms = "numeric",
## - corr.terms = "list",
## - corr.terms.tot = "list"
## 

##
## a "combination" is stored as a list with fields:
##
## - tag = "character",
## - bibitem = "character",
## - quantities = "character",
## - params = "list",
## - meas.lin.combs = "list"
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

iline = 1
repeat {
  lines.len = length(lines)
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
    next
  }
  iline = iline + 1
}

##
## define storage for measurements and one combination
##
measurements = list()
combination = list()
meas.options = list()

flag.in.meas = FALSE
flag.in.data = FALSE
flag.in.params = FALSE
flag.in.combine = FALSE
flag.in.sumofquant = FALSE
flag.in.combofquant = FALSE
flag.in.constraint = FALSE
flag.in.combine.meas = FALSE

data.labels = character()
data.values = numeric()
meas.labels = character()
meas.values = numeric()

for (line in lines) {
  ## cat(line,"\n")
  if (regexpr("^\\s*$", line, perl=TRUE) != -1 ||
      regexpr("^[*#;]", line, perl=TRUE) != -1) {
    next
  }
  line = gsub("\\s*!.*", "", line, perl=TRUE)
  fields = unlist(strsplit(line, "\\s+", perl=TRUE))

  ##--- first field == whitespace means it is a continuation line
  flag.continuation = match.nocase("^\\s*$", fields[1])
  
  ##--- begin of measurement cards
  if (match.nocase("^BEGIN$", fields[1])) {
    if (flag.in.meas) {
      cat("error, BEGIN keyword inside definition of", meas$tag, "\n")
    } else {
      ##--- init measurement list
      meas = list()
      meas$tag = paste(as.character(fields[2:4]), collapse=".")
      meas$bibitem = as.character(fields[-1])
      meas$value = numeric()
      meas$stat = numeric()
      meas$syst = numeric()
      meas$params = list()
      meas$syst.terms = numeric()
      meas$corr.terms = list()
      meas$corr.terms.tot = list()

      sumofquant.values = numeric()
      combofquant.labels = character()
      combofquant.values = numeric()
      measlincombs.list = list()

      constraint.labels = character()
      constraint.values = numeric()

      constraints.list.comb = list()
      constraints.list.val = list()

      nlconstr.comb = list()
      nlconstr.val = list()

      flag.in.meas = TRUE
    }
    next
  }
  if (!flag.in.meas) {
    cat("error, ", fields[1], "outside a measurement definition (BEGIN..END)\n")
    next
  }
  ##
  ## field[1] != <white space> means end of data for current data card field
  ##
  if (!flag.continuation) {
    ##
    ## store collected data of SUMOFQUANT cards
    ##
    if (flag.in.sumofquant) {
      ##--- vector with one for each quantity whose sum corresponds to the measurement
      val = rep(1, length(sumofquant.values)-1)
      names(val) = tail(sumofquant.values, -1)
      ##--- list of measurements with the coefficients corresponding to the related quantities
      measlincombs.list[[sumofquant.values[1]]] = val
      ##--- reset vector of SUMOFQUANT parameter
      sumofquant.values = numeric()
    }
    ##
    ## store collected data of COMBOFQUANT cards
    ##
    if (flag.in.combofquant) {
      if (length(combofquant.labels)-1 != length(combofquant.values)) {
        stop("error: comb.lin. coeffs ", length(combofquant.values), " while expecting ", length(combofquant.labels)-1)
      }
      names(combofquant.values) = tail(combofquant.labels, -1)
      measlincombs.list[[combofquant.labels[1]]] = combofquant.values
      ##--- reset vectors of COMBOFQUANT parameter
      combofquant.labels = character()
      combofquant.values = numeric()
    }
    ##
    ## store collected data of CONSTRAINT cards
    ##
    if (flag.in.constraint) {
      if (length(constraint.labels) != length(constraint.values)) {
        stop("error: comb.lin. coeffs ", length(constraint.values), " while expecting ", length(constraint.labels))
      }
      if (length(constraint.labels) <= 1) {
        stop("error: no quantities listed for constraint ", constraint.labels[1])
      }
      if (!is.null(constraints.list.comb[[constraint.labels[1]]]) ||
          !is.null(constraints.list.val[[constraint.labels[1]]])) {
        stop("error: constraint entered twice: ", constraint.labels[1])
      }
      names(constraint.values) = constraint.labels
      constraints.list.comb[[constraint.labels[1]]] = constraint.values[-1]
      tmp = constraint.values[1]
      names(tmp) = NULL
      constraints.list.val[[constraint.labels[1]]] = tmp
      ##--- reset vectors of CONSTRAINT parameter
      constraint.labels = character()
      constraint.values = numeric()
    }
    ##
    ## store collected data of MEASUREMENT card in COMBINE section
    ##
    if (flag.in.combine.meas) {
      ##--- add this measurement to the list of measurments to combine
      meas.label = sub("^m_","", data.labels[1])
      meas.labels = c(meas.labels, meas.label)
      if ((length(data.labels) - 3) != length(data.values)) {
        stop("error: MEASUREMENT ", meas.label, " mismatch between labels and data values")
      }
      names(data.values) = data.labels[-(1:3)]
      if (length(data.values) > 0) {
        meas.options[[meas.label]] = data.values
      }
      data.labels = character()
      data.values = numeric()
    }
    flag.in.data = FALSE
    flag.in.params = FALSE
    flag.in.sumofquant = FALSE
    flag.in.combofquant = FALSE
    flag.in.constraint = FALSE
    flag.in.combine.meas = FALSE
  }
  if (match.nocase("^END$", fields[1])) {
    if (!flag.in.combine) {
      ##
      ## measurement cards
      ##
      if (length(meas.labels) != 3) {
        stop("wrong number of MEASUREMENT labels: ", length(meas.labels), "instead of 3\n")
      }
      if (length(data.labels) != length(data.values)) {
        stop("number of numeric items does not match number of labels\n")
      }
      ##
      ## deal with stat & syst errors expressed as percentage of value
      ##++ assume that first data value is the measurement
      ##
      patt.perc = "([[:alnum:]]+[^[:alnum:]]*)(%)([^[:alnum:]]*)$"
      for (i in 2:length(data.values)) {
        if (regexpr(patt.perc, data.labels[i]) != -1) {
          data.labels[i] = gsub(patt.perc, "\\1\\3", data.labels[i])
          data.values[i] = data.values[1] * data.values[i] /100
        }
      }
      names(data.values) = data.labels
      
      ##--- get measurement values from following DATA values
      meas.values = numeric(length(meas.labels))
      names(meas.values) = meas.labels
      for (meas.label in meas.labels) {
        if (is.na(data.values[meas.label])) {
          stop("missing MEASUREMENT data for ", meas.label, "\n")
        }
        meas.values[meas.label] = data.values[meas.label]
      }

      ##--- edit measurement value, stat, syst labels, preserve meas.labels
      labels = meas.labels
      labels = gsub("^m_", "", labels, ignore.case=TRUE)
      labels = gsub("statistical", "stat", labels, ignore.case=TRUE)
      labels = gsub("systematic", "syst", labels, ignore.case=TRUE)

      if (meas$bibitem[2] != labels[1]) {
        ##
        ## The Combos field "method" was probably intended to indicate the method with which
        ## a quantity was measured.  However, when averaging measurements of different quantities
        ## possibly statistically correlated, the STAT_CORR_WITH card requires the indication
        ## of "experiment" "method" "where".  In this case "method" should actually indicate what
        ## quantity is measured. Different methods labels should be merged with the "where" label.
        ## 
        ## Here we suggest that "method" must actually indicate what quantity is measured.
        ## Due to Combos, the MEASUREMENT label must be surreptitiously set to a quantity other than the measured one
        ## in order to include a measurement that does not measure any of the quantities that are to be combined
        ## but rather a linear combination of them. A warning is issued when that happens, and the MEASUREMENT
        ## label is set to the "method", which actually should indicate what is being measured.
        ##
        cat("warning: label '", labels[1], "' != method '", paste(meas$bibitem[1:3], collapse="."), "': replaced\n", sep="")
        labels[1] = meas$bibitem[2]
      }

      ##--- set labels for measurement and stat./syst. errors
      names(meas.values) = labels

      ##--- store MEASUREMENT DATA
      meas$value = meas.values[1]
      meas$stat = meas.values[2]
      meas$syst = meas.values[3]
      meas$syst.terms = data.values[!data.labels %in% meas.labels]
      if (!is.null(measurements[[meas$tag]])) {
        ##--- two measurements have the same tag
        stop(paste("error: two measurements with tag", meas$tag))
      }
      measurements[[meas$tag]] = meas

      ##--- reset temporary storage for data values and labels
      data.labels = character()
      data.values = numeric()
      meas.labels = character()
      meas.values = numeric()
    } else {
      ##
      ## combination cards
      ##
      combination$tag = meas$tag
      combination$bibitem = meas$bibitem
      combination$quantities = meas.labels
      combination$quantities.options = meas.options
      combination$params = meas$params
      combination$meas.lin.combs = measlincombs.list
      combination$constr.comb = constraints.list.comb
      combination$constr.val = constraints.list.val
      combination$nlconstr.comb = nlconstr.comb
      combination$nlconstr.val = nlconstr.val
    }
    flag.in.meas = FALSE
    flag.in.combine = FALSE
    next
  }
  ##--- MEASUREMENT in measurement cards
  if (match.nocase("^MEASUREMENT$", fields[1]) && !flag.in.combine) {
    ##--- get labels for value, stat.error, syst.error
    meas.labels = fields[2:4]
    next
  }
  if (match.nocase("^COMBINE$", fields[1])) {
    flag.in.combine = TRUE
    next
  }
  if (match.nocase("^DATA$", fields[1]) ||
      (flag.in.data && flag.continuation)) {
    flag.in.data = TRUE
    data.labels = c(data.labels,
      unlist(lapply(fields[-1], function(elem) {if (is.na(suppressWarnings(as.numeric(elem)))) {elem}})))
    data.values = c(data.values,
      unlist(lapply(fields[-1], function(elem) {if (!is.na(suppressWarnings(as.numeric(elem)))) {as.numeric(elem)}})))
    next
  }
  if (match.nocase("^STAT_CORR_WITH$", fields[1])) {
    corr = as.numeric(fields[5])
    names(corr) =  paste(as.character(fields[2:4]), collapse=".")
    meas$corr.terms = c(meas$corr.terms,corr)
    next
  }
  if (match.nocase("^ERROR_CORR_WITH$", fields[1])) {
    corr = as.numeric(fields[5])
    names(corr) =  paste(as.character(fields[2:4]), collapse=".")
    meas$corr.terms.tot = c(meas$corr.terms.tot, corr)
    next
  }
  if (match.nocase("^PARAMETERS$", fields[1]) ||
      (flag.in.params && flag.continuation)) {
    flag.in.params = TRUE
    params.data = lapply(fields[-1], function(x) type.convert(x, as.is=TRUE))
    if (length(params.data) == 0) next
    if (!is.character(params.data[[1]]) ||
        !is.numeric(params.data[[2]]) ||
        !is.numeric(params.data[[3]]) ||
        !is.numeric(params.data[[4]])) {
      stop("wrong parameter data ", paste(params.data,sep=","), "\n")
    }
    names(params.data)[2:4] = c("value", "delta_pos", "delta_neg")
    meas$params = c(meas$params, list(unlist(params.data[2:4])))
    names(meas$params)[length(meas$params)] = params.data[[1]]
    next
  }
  ##--- non-linear constraints, label and expression
  if (flag.in.combine && match.nocase("^NLCONSTRAINT$", fields[1])) {
    expr = paste(fields[-(1:3)], collapse="")
    expr = gsub("\"?([^\"]*)\"?.*", "\\1", expr, perl=TRUE)
    nlconstr.comb[[fields[2]]] = expr
    nlconstr.val[[fields[2]]] = ifelse(!is.na(suppressWarnings(as.numeric(fields[3]))), as.numeric(fields[3]), 0)
  }
  ##--- MEASUREMENT in COMBINE cards
  if ((flag.in.combine && match.nocase("^MEASUREMENT$", fields[1])) ||
      (flag.in.combine.meas && flag.continuation)) {
    flag.in.combine.meas = TRUE
    ##--- collect data in the MEASUREMENT line
    data.labels = c(data.labels,
      unlist(lapply(fields[-1], function(elem) {if (is.na(suppressWarnings(as.numeric(elem)))) {elem}})))
    data.values = c(data.values,
      unlist(lapply(fields[-1], function(elem) {if (!is.na(suppressWarnings(as.numeric(elem)))) {as.numeric(elem)}})))
    next
  }
  ##
  ## SUMOFQUANT <quant> <quant 1> [<quant 2> ...]
  ## define a measurement as sum of quantities to fit
  ##
  if (match.nocase("^SUMOFQUANT$", fields[1]) ||
      (flag.in.sumofquant && flag.continuation)) {
    flag.in.sumofquant = TRUE
    sumofquant.values = c(sumofquant.values, as.character(fields[-1]))
  }
  ##
  ## COMBOFQUANT <quant> <quant 1> <coeff 1> [<quant 2> <coeff 2> ...]
  ## define a measurement as linear combination of quantities to fit
  ##
  if (match.nocase("^COMBOFQUANT$", fields[1]) ||
      (flag.in.combofquant && flag.continuation)) {
    flag.in.combofquant = TRUE
    combofquant.labels = c(combofquant.labels,
      unlist(lapply(fields[-1], function(elem) {if (is.na(suppressWarnings(as.numeric(elem)))) {elem}})))
    combofquant.values = c(combofquant.values,
      unlist(lapply(fields[-1], function(elem) {if (!is.na(suppressWarnings(as.numeric(elem)))) {as.numeric(elem)}})))
  }
  ##
  ## CONSTRAINT <name> <value> <quant 1> <coeff 1> [<quant 2> <coeff 2> ...]
  ## define a measurement as linear combination of quantities to fit
  ##
  if (match.nocase("^CONSTRAINT$", fields[1]) ||
      (flag.in.constraint && flag.continuation)) {
    flag.in.constraint = TRUE
    constraint.labels = c(constraint.labels,
      unlist(lapply(fields[-1], function(elem) {if (is.na(suppressWarnings(as.numeric(elem)))) {elem}})))
    constraint.values = c(constraint.values,
      unlist(lapply(fields[-1], function(elem) {if (!is.na(suppressWarnings(as.numeric(elem)))) {as.numeric(elem)}})))
  }
}

list(combination=combination, measurements=measurements)
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
