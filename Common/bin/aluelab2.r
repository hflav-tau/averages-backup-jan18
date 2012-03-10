#!/usr/bin/env Rscript

require(proto, quietly=TRUE)

## ////////////////////////////////////////////////////////////////////////////
## definitions

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

## ////////////////////////////////////////////////////////////////////////////
## definitions

##
## class to do computations on statistically correlated quantities
##
StatComb = proto()

##--- create object to store quantities values and covariance
StatComb$new = function(., val=numeric(0), cov=matrix(ncol=0, nrow=0)) {
  proto(., val=val, cov=cov)
}

StatComb$val = function(.) {
  .$val
}

StatComb$cov = function(.) {
  .$cov
}

StatComb$corr = function(.) {
  err = sqrt(diag(.$cov))
  .$cov / (err %o% err)
}

StatComb$err = function(.) {
  sqrt(diag(.$cov))
}

StatComb$val.err = function(., name=NULL) {
  if (is.null(name)) return(NULL)
  c(val=.$val[name], err=sqrt(.$cov[name, name]))
}

##
## add a vector of measurements and their correlation to the list
## to quant.val, quant.err, quant.corr, quant.cov
##
## aeb.meas.add = function(add.val, add.err, add.corr=NULL) {}
StatComb$meas.add = function(., add.val, add.err, add.corr=NULL) {
  ##
  ## if a correlation matrix is found, use it to determine which quantities where averaged
  ## and which quantities are non-averaged combined quantities
  ## retain only averaged quantities and their correlation
  ##
  if (!is.null(add.corr)) {
    quant.names.averaged = rownames(add.corr)
  } else {
    quant.names.averaged = names(add.val)
    add.corr = as.matrix(diag.m(rep(1, length(add.val))))
    rownames(add.corr) = names(add.val)
    colnames(add.corr) = names(add.val)
  }

  add.cov = add.corr * (add.err %o% add.err)
  
  ##--- assemble covariance matrix
  cov.right = matrix(0, dim(.$cov)[1], dim(add.cov)[2])
  colnames(cov.right) = colnames(add.cov)
  cov.top = cbind(.$cov, cov.right)
  cov.left = matrix(0, dim(add.cov)[1], dim(.$cov)[2])
  rownames(cov.left) = rownames(add.cov)
  cov.bottom = cbind(cov.left, add.cov)
  .$cov = rbind(cov.top, cov.bottom)
  
  ##--- assemble averaged quantities
  quant.averaged.sel = names(add.val) %in% quant.names.averaged
  .$val = c(.$val, add.val[quant.averaged.sel])
}

##
## add a single additional uncorrelated measurement
## to quant.val, quant.err, quant.corr, quant.cov
##
## aeb.meas.add.single = function(label, val, err) {}
StatComb$meas.add.single = function(., label, val, err) {
  val = as.numeric(val)
  err = as.numeric(err)
  names(val) = label
  names(err) = label
  .$meas.add(val, err)
}

##
## compute linear combination given measurements and their variance
##
StatComb$linear.comb.with.cov = function(., lc, val, cov) {
  lc.val = val %*% lc
  lc.var = drop(lc %*% cov %*% lc)
  lc.err = sign(lc.var)*sqrt(abs(lc.var))
  return(c(val=lc.val, err=lc.err))
}

##
## compute linear combinations using measurements stored
## in quant.val, quant.cov
##
## aeb.linear.comb.glob = function(lc) {}
StatComb$linear.comb = function(., label, val, err) {
  diff = setdiff(names(lc), names(.$val))
  if (length(diff) > 0) {
    stop("error: following quantities were not loaded: ", diff)
  }
  meas.lc = names(lc)
  return(.$linear.comb.with.cov(lc, .$val[meas.lc], .$cov[meas.lc, meas.lc]))
}

##
## add quantity, specifying its gradient to compute the covariance
## - add value, error, covariamce, correlation
##
## aeb.meas.val.grad.add = function(add.name, add.val, add.grad) {}
StatComb$meas.val.grad.add = function(., add.name, add.val, add.grad) {
  names(add.val) = add.name
  add.comb = drop(add.grad)
  add.comb.full = .$val * 0
  add.comb.full[names(add.comb)] = add.grad
  add.cov = add.comb.full %*% .$cov %*% add.comb.full
  names(add.cov) = add.name
  .$cov = rbind(
    cbind(.$cov, matrix(.$cov %*% add.comb.full, dimnames=list(NULL, add.name))),
    matrix(c(add.comb.full %*% .$cov, add.cov), 1, dim(.$cov)[2]+1, dimnames=list(add.name)))
  .$val = c(.$val, add.val)
  return(invisible())
}

##
## add quantity that is combination of other quantities
## - add value, error, covariamce, correlation
##
## aeb.meas.comb.add = function(add.name, add.comb) {}
StatComb$meas.comb.add = function(., add.name, add.comb) {
  add.comb = drop(add.comb)
  add.comb.full = .$val * 0
  add.comb.full[names(add.comb)] = add.comb
  add.val = add.comb.full %*% .$val
  names(add.val) = add.name
  add.err = sqrt(add.comb.full %*% .$cov %*% add.comb.full)
  names(add.err) = add.name
  .$val = c(.$val, add.val)
  .$cov = rbind(
    cbind(.$cov, matrix(.$cov %*% add.comb.full, dimnames=list(NULL, add.name))),
    matrix(c(add.comb.full %*% .$cov, add.err^2), 1, dim(.$cov)[2]+1, dimnames=list(add.name)))
  return(invisible())
}

##
## add quantity defined as expression of existing quantities
##
## aeb.meas.expr.add = function(add.name, add.expr) {}
StatComb$meas.expr.add = function(., add.name, add.expr) {
  add.deriv.expr = deriv(add.expr, all.vars(add.expr))
  add.val = eval(add.deriv.expr, as.list(.$val))
  add.grad = attr(add.val, "gradient")
  .$meas.val.grad.add(add.name, add.val, add.grad)
  return(invisible())
}

##
## given a model matrix that describes how a parameter relates to quantities,
## fit for it and return the fitted combination of quantities
##
## aeb.model.matrix.fit = function(model.matrix) {}
StatComb$model.matrix.fit = function(., model.matrix) {
  model.matrix.names = names(model.matrix)
  val = .$val[model.matrix.names]
  cov = .$cov[model.matrix.names, model.matrix.names]
  model.matrix = matrix(model.matrix, length(model.matrix), 1)
  invcov = solve(cov)
  fit.cov = solve(t(model.matrix) %*% invcov %*% model.matrix)
  fit.comb = fit.cov %*% t(model.matrix) %*% invcov
  return(fit.comb)
}

##
## fit theory parameter according to model matrix and
## add the result to the quantities
##
## aeb.meas.fit.add = function(add.name, add.model.matrix) {}
StatComb$meas.fit.add = function(., add.name, add.model.matrix) {
  fit.comb = .$model.matrix.fit(add.model.matrix)
  .$meas.comb.add(add.name, fit.comb)
  return(invisible())
}
