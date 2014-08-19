## ////////////////////////////////////////////////////////////////////////////
##
## aluelab2.r
##
## utility functions to elaborate alucomb2.r results
##

require(methods, quietly=TRUE)

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

##--- sum in quadrature
quadrature = function(x) {
  return(sqrt(sum(x^2)))
}

##
## substitute after evaluating arg
## copied from alucomb2-utils.r
##
esub = function(expr, sublist=NULL) do.call("substitute", list(expr, sublist))
esub.expr = function(expr, sublist=NULL) {
  sapply(as.expression(expr), function(call) as.expression(esub(call, sublist)))
}

##
## deparse expression and produce single line
## copied from alucomb2-utils.r
##
deparse.one.line = function(expr) {
  paste(gsub("^\\s+|\\s+$", "", sapply(as.expression(expr), function(x) deparse(x)), perl=TRUE), collapse="")
}

## ////////////////////////////////////////////////////////////////////////////
##
## class to do computations on statistically correlated quantities
##

##--- class to store quantities values, their covariance, a parameter list
StatComb = setRefClass("StatComb",
  fields = list(
    .val = "numeric",
    .cov = "matrix",
    .param = "numeric"
    ),
  methods=list(
    initialize = function(quant.val=numeric(), quant.cov=matrix(numeric(0),0,0), parameters=numeric()) {
      ## callSuper(...)
      .self$.val = quant.val
      .self$.cov = quant.cov
      .self$.param = parameters
      .self
    })
  )

rc = StatComb$methods(
  param.add = function(param) {
    .self$.param = c(.param, param)
  })

rc = StatComb$methods(
  param = function(name=NULL) {
    if (is.null(name)) return(.param)
    .param[name]
  })

rc = StatComb$methods(
  param.err = function(name=NULL) {
    if (is.null(name)) return(.param * 0)
    .param[name] * 0
  })

rc = StatComb$methods(
  val = function(name=NULL) {
    if (is.null(name)) return(.val)
    .val[name]
  })

rc = StatComb$methods(
  cov = function(rows=NULL, cols=NULL) {
    if (is.null(cols)) cols=rows
    if (is.null(rows)) return(.cov)
    .cov[rows, cols]
  })

rc = StatComb$methods(
  corr = function(rows=NULL, cols=NULL) {
    if (is.null(cols)) cols=rows
    err = sqrt(diag(.cov))
    corr = .cov / (err %o% err)
    if (is.null(rows)) {
      return(corr)
    }
    corr[rows, cols]
  })

rc = StatComb$methods(
  err = function(name=NULL) {
    if (is.null(name)) return(sqrt(diag(.cov)))
    sqrt(diag(.cov))[name]
  })

rc = StatComb$methods(
  val.err = function(name=NULL) {
    if (is.null(name)) return(NULL)
    c(val=.val[name], err=sqrt(.cov[name, name]))
  })

rc = StatComb$methods(
  vals = function(name=NULL) {
    if (is.null(name)) return(c(.val, .param))
    c(.val, .param)[name]
  })

rc = StatComb$methods(
  errs = function(name=NULL) {
    if (is.null(name)) return(c(err(), param.err()))
    c(err(), param.err())[name]
  })

rc = StatComb$methods(
  vnames = function() {
    c(names(.val), names(.param))
  })

##
## get string expression for a quantity, from its constraint equations
## needs alucomb2 format
##
rc = StatComb$methods(
  str.to.comb = function(str.expr) {
    expr = parse(text=as.character(str.expr))
    vars = intersect(all.vars(expr), names(.val))
    comb = drop(attr(eval(deriv(expr, vars), c(as.list(.val), as.list(.param))), "gradient"))
    return(comb)
  })

##
## add a vector of measurements and their correlation to the list
## to quant.val, quant.err, quant.corr, quant.cov
##
rc = StatComb$methods(
  meas.add = function(add.val, add.err, add.corr=NULL) {
    if (length(add.val) != length(add.err)) stop("mismatch of dimensions of add.val and add.err")
    quant.names = names(add.val)

    if (!is.null(add.corr)) {
      if (length(add.val) != dim(add.corr)[1] || length(add.val) != dim(add.corr)[2]) {
        stop("mismatch of dimensions of add.val and add.corr")
      }
      ##--- don't check for now if names match
    } else {
      add.corr = as.matrix(diag.m(rep(1, length(add.val))))
      rownames(add.corr) = names(add.val)
      colnames(add.corr) = names(add.val)
      ##--- don't check for now if names match
    }
    add.cov = add.corr * (add.err %o% add.err)
    
    ##--- assemble covariance matrix
    cov.right = matrix(0, dim(.cov)[1], dim(add.cov)[2])
    colnames(cov.right) = colnames(add.cov)
    cov.top = cbind(.cov, cov.right)
    cov.left = matrix(0, dim(add.cov)[1], dim(.cov)[2])
    rownames(cov.left) = rownames(add.cov)
    cov.bottom = cbind(cov.left, add.cov)
    .self$.cov = rbind(cov.top, cov.bottom)
    
    ##--- assemble values
    .self$.val = c(.val, add.val)
    return(invisible(NULL))
  })

##
## add a single additional parameter
##
rc = StatComb$methods(
  param.add.single = function(label, val) {
    val = as.numeric(val)
    names(val) = label
    .self$.param = c(.param, val)
    return(val)
  })

##
## add a single additional uncorrelated measurement
## to quant.val, quant.err, quant.corr, quant.cov
##
rc = StatComb$methods(
  meas.add.single = function(label, val, err=0) {
    val = as.numeric(val)
    err = as.numeric(err)
    if (err == 0) return(param.add.single(label, val))
    names(val) = label
    names(err) = label
    meas.add(val, err)
    return(c(val, err))
  })

##
## add a single correlation between two quantities
##
rc = StatComb$methods(
  corr.add.single = function(label.1, label.2, corr) {
    corr = as.numeric(corr)
    cov = corr*sqrt(.cov[label.1, label.1] * .cov[label.2, label.2])
    .self$.cov[label.1, label.2] = cov
    .self$.cov[label.2, label.1] = cov
    return(corr)
  })

##
## compute linear combination given measurements and their variance
##
rc = StatComb$methods(
  linear.comb.with.cov = function(lc, val, cov) {
    lc.val = val %*% lc
    lc.var = drop(lc %*% cov %*% lc)
    lc.err = sign(lc.var)*sqrt(abs(lc.var))
    return(c(val=lc.val, err=lc.err))
  })

##
## compute linear combinations using measurements stored
## in quant.val, quant.cov
##
rc = StatComb$methods(
  linear.comb = function(label, val, err) {
    diff = setdiff(names(lc), names(.val))
    if (length(diff) > 0) {
      stop("error: following quantities were not loaded: ", diff)
    }
    meas.lc = names(lc)
  return(linear.comb.with.cov(lc, .val[meas.lc], .cov[meas.lc, meas.lc]))
  })

##
## add quantity, specifying its gradient to compute the covariance
## - add value, error, covariamce, correlation
##
rc = StatComb$methods(
  meas.val.grad.add = function(add.name, add.val, add.grad) {
    names(add.val) = add.name
    add.comb = drop(add.grad)
    add.comb.full = .val * 0
    add.comb.full[names(add.comb)] = add.grad
    add.cov = add.comb.full %*% .cov %*% add.comb.full
    names(add.cov) = add.name
    .self$.cov = rbind(
      cbind(.cov, matrix(.cov %*% add.comb.full, dimnames=list(NULL, add.name))),
      matrix(c(add.comb.full %*% .cov, add.cov), 1, dim(.cov)[2]+1, dimnames=list(add.name)))
    .self$.val = c(.val, add.val)
    return(invisible())
  })

##
## add quantity that is combination of other quantities
## - add value, error, covariamce, correlation
##
rc = StatComb$methods(
  meas.comb.add = function(add.name, add.comb) {
    add.comb = drop(add.comb)
    add.comb.full = .val * 0
    add.comb.full[names(add.comb)] = add.comb
    add.val = add.comb.full %*% .val
    names(add.val) = add.name
    add.err = sqrt(add.comb.full %*% .cov %*% add.comb.full)
    names(add.err) = add.name
    .self$.val = c(.val, add.val)
    .self$.cov = rbind(
      cbind(.cov, matrix(.cov %*% add.comb.full, dimnames=list(NULL, add.name))),
      matrix(c(add.comb.full %*% .cov, add.err^2), 1, dim(.cov)[2]+1, dimnames=list(add.name)))
    return(invisible())
  })

##
## add quantity defined as expression of existing quantities
##
rc = StatComb$methods(
  meas.expr.add = function(add.name, add.expr) {
    add.expr = substitute(add.expr)
    ##--- substitute parameters
    add.expr = esub.expr(add.expr, list(.param))
    add.deriv.expr = deriv(add.expr, all.vars(add.expr))
    add.val = eval(add.deriv.expr, as.list(.val))
    add.grad = attr(add.val, "gradient")
    meas.val.grad.add(add.name, add.val, add.grad)
    return(invisible())
  })

##
## given a model matrix that describes how a parameter relates to quantities,
## fit for it and return the fitted combination of quantities
##
rc = StatComb$methods(
  model.matrix.fit = function(model.matrix) {
    model.matrix.names = names(model.matrix)
    val = .val[model.matrix.names]
    cov = .cov[model.matrix.names, model.matrix.names]
    model.matrix = matrix(model.matrix, length(model.matrix), 1)
    invcov = solve(cov)
    fit.cov = solve(t(model.matrix) %*% invcov %*% model.matrix)
    fit.comb = fit.cov %*% t(model.matrix) %*% invcov
    return(fit.comb)
  })

##
## fit theory parameter according to model matrix and
## add the result to the quantities
##
rc = StatComb$methods(
  meas.fit.add = function(add.name, add.model.matrix) {
    fit.comb = model.matrix.fit(add.model.matrix)
    meas.comb.add(add.name, fit.comb)
    return(invisible())
  })

##
## compute systematic contribution
##
rc = StatComb$methods(
  syst.contrib = function(quant.name, ...) {
    syst.contrib.list = list(...)
    if (length(syst.contrib.list) == 0) {
      return(0);
    } else if (length(syst.contrib.list) == 1) {
      syst.name = syst.contrib.list[[1]]
      return(.cov[quant.name, syst.name]/err()[syst.name])
    } else {
      return(sqrt(sum(sapply(syst.contrib.list, function(syst.name) {.cov[quant.name, syst.name]/err()[syst.name]})^2)))
    }
  })

rc = StatComb$methods(
  syst.contrib.rel = function(quant.name, ...) {
    syst.contrib(quant.name, ...) / .val[quant.name]
  })

rc = StatComb$methods(
  syst.contrib.perc = function(quant.name, ...) {
    syst.contrib(quant.name, ...) / .val[quant.name] * 100
  })
  
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
