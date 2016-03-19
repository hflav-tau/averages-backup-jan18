## ////////////////////////////////////////////////////////////////////////////
##
## aluelab3.r
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
esub.single = function(expr, sublist=NULL) do.call(substitute, list(expr, sublist))
esub.expr = function(expr, sublist=NULL) {
  sapply(as.expression(expr), function(call) as.expression(esub.single(call, sublist)))
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
    .param = "numeric",
    .gamma = "character",
    .texdescr = "character"
    ),
  methods=list(
    initialize = function(
      quant.val=numeric(),
      quant.cov=matrix(numeric(0),0,0),
      parameters=numeric()
      ) {
      ## callSuper(...)
      .self$.val = quant.val
      .self$.cov = quant.cov
      .self$.param = parameters
      .self$.gamma = character(0)
      .self$.texdescr = character(0)
      .self
    })
  )

rc = StatComb$methods(
  gamma.add = function(gamma) {
    gamma = unlist(gamma, use.names=TRUE)
    .self$.gamma[names(gamma)] = gamma
    return(invisible())
  })

rc = StatComb$methods(
  texdescr.add = function(texdescr) {
    texdescr = unlist(texdescr, use.names=TRUE)
    .self$.texdescr[names(texdescr)] = texdescr
    return(invisible())
  })

rc = StatComb$methods(
  gamma.add.single = function(label, gamma=NULL) {
    if (!is.null(gamma)) {
      .self$.gamma[label] = gamma
    }
    return(invisible())
  })

rc = StatComb$methods(
  texdescr.add.single = function(label, texdescr=NULL) {
    if (!is.null(texdescr)) {
      .self$.texdescr[label] = texdescr
    }
    return(invisible())
  })

rc = StatComb$methods(
  param.add = function(param, gamma = NULL, texdescr=NULL) {
    .self$.param = c(.param, unlist(param, use.names=TRUE))
    if (!is.null(gamma)) {
      gamma = unlist(gamma, use.names=TRUE)
      .self$.gamma[names(gamma)] = gamma
    }
    if (!is.null(texdescr)) {
      texdescr = unlist(texdescr, use.names=TRUE)
      .self$.texdescr[names(texdescr)] = texdescr
    }
  })

rc = StatComb$methods(
  gamma = function(name=NULL) {
    if (is.null(name)) return(.gamma)
    .gamma[name]
  })

rc = StatComb$methods(
  texdescr = function(name=NULL) {
    if (is.null(name)) return(.texdescr)
    .texdescr[name]
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
  vnames = function(name=NULL) {
    if (is.null(name)) return(c(names(.val), names(.param)))
    c(names(.val), names(.param))[name]
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
  quant.add = function(add.val, add.err, add.corr=NULL, gamma=NULL, texdescr=NULL) {
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

    if (!is.null(gamma)) {
      gamma = unlist(gamma, use.names=TRUE)
      .self$.gamma[names(gamma)] = gamma
    }

    if (!is.null(texdescr)) {
      texdescr = unlist(texdescr, use.names=TRUE)
      .self$.texdescr[names(texdescr)] = texdescr
    }
    
    return(invisible(NULL))
  })

##
## add a single additional parameter
##
rc = StatComb$methods(
  param.add.single = function(label, val, gamma=NULL, texdescr=NULL) {
    val = as.numeric(val)
    names(val) = label
    .self$.param = c(.param, val)
    gamma.add.single(label, gamma)
    texdescr.add.single(label, texdescr)
    return(invisible())
  })

##
## add a single additional uncorrelated measurement
## to quant.val, quant.err, quant.corr, quant.cov
##
rc = StatComb$methods(
  quant.add.single = function(label, val, err=0, gamma=NULL, texdescr=NULL) {
    val = as.numeric(val)
    err = as.numeric(err)
    if (err == 0) return(param.add.single(label, val))
    names(val) = label
    names(err) = label
    if (!is.null(gamma)) {
      names(gamma) = label
    }
    if (!is.null(texdescr)) {
      names(texdescr) = label
    }
    quant.add(val, err, gamma, texdescr)
    gamma.add.single(label, gamma)
    texdescr.add.single(label, texdescr)
    return(invisible())
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
    return(invisible())
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
    quant.lc = names(lc)
  return(linear.comb.with.cov(lc, .val[quant.lc], .cov[quant.lc, quant.lc]))
  })

##
## add quantity, specifying its gradient to compute the covariance
## - add value, error, covariamce, correlation
##
rc = StatComb$methods(
  quant.val.grad.add = function(add.name, add.val, add.grad, gamma=NULL, texdescr=NULL) {
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
    gamma.add.single(add.name, gamma)
    texdescr.add.single(add.name, texdescr)
    return(invisible())
  })

##
## add quantity that is combination of other quantities
## - add value, error, covariamce, correlation
##
rc = StatComb$methods(
  quant.comb.add = function(add.name, add.comb, gamma=NULL, texdescr=NULL) {
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
    gamma.add.single(add.name, gamma)
    texdescr.add.single(add.name, texdescr)
    return(invisible())
  })

##
## add quantity defined as expression of existing quantities
## let the arg be resolved before entering this function
##
rc = StatComb$methods(
  quant.qexpr.add = function(add.name, add.expr, gamma=NULL, texdescr=NULL) {
    ##--- substitute parameters
    add.expr = esub.expr(add.expr, as.list(.param))
    add.deriv.expr = deriv(add.expr, all.vars(add.expr))
    add.val = eval(add.deriv.expr, as.list(.val))
    add.grad = attr(add.val, "gradient")
    quant.val.grad.add(add.name, add.val, add.grad)
    gamma.add.single(add.name, gamma)
    texdescr.add.single(add.name, texdescr)
    return(invisible())
  })

##
## add quantity defined as expression of existing quantities
##
rc = StatComb$methods(
  quant.expr.add = function(add.name, add.expr, gamma=NULL, texdescr=NULL) {
    add.expr = substitute(add.expr)
    quant.qexpr.add(add.name, add.expr)
    gamma.add.single(add.name, gamma)
    texdescr.add.single(add.name, texdescr)
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
  quant.fit.add = function(add.name, add.model.matrix, gamma=NULL, texdescr=NULL) {
    fit.comb = model.matrix.fit(add.model.matrix)
    quant.comb.add(add.name, fit.comb)
    gamma.add.single(add.name, gamma)
    texdescr.add.single(add.name, texdescr)
    return(invisible())
  })

##
## compute systematic contribution
##
rc = StatComb$methods(
  err.contrib = function(quant.name, ...) {
    err.contrib.list = list(...)
    if (length(err.contrib.list) == 0) {
      return(err(quant.name));
    } else if (length(err.contrib.list) == 1) {
      syst.name = err.contrib.list[[1]]
      return(.cov[quant.name, syst.name]/err(syst.name))
    } else {
      return(sqrt(sum(sapply(err.contrib.list, function(syst.name) {.cov[quant.name, syst.name]/err(syst.name)})^2)))
    }
  })

rc = StatComb$methods(
  err.contrib.rel = function(quant.name, ...) {
    err.contrib(quant.name, ...) / .val[quant.name]
  })

rc = StatComb$methods(
  err.contrib.perc = function(quant.name, ...) {
    err.contrib(quant.name, ...) / .val[quant.name] * 100
  })
