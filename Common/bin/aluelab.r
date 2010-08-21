## ////////////////////////////////////////////////////////////////////////////
##
## aluelab.r
##
## - get averages from alucomb.r log files
## - stores values, errors and correlation in global variables for further elaboration
##
## ////////////////////////////////////////////////////////////////////////////

library(methods)

source("../../../Common/bin/alu-utils.r")

## ////////////////////////////////////////////////////////////////////////////
## definitions

##-- convert current directory to the corresponding data (or log) directory
aeb.log.dir = file.path("../../../Data", sub("^.*/([^/]*/[^/]*/[^/]*)$", "\\1", getwd()))
aeb.log.dir.base = dirname(aeb.log.dir)

##-- if CL less than minimum, inflate variance with S-factors
aeb.cl.min = pnorm(3, mean=0, sd=1, lower.tail=FALSE)

## ////////////////////////////////////////////////////////////////////////////
## code

##-- global variables for averages and their variance
quant.val = numeric(0)
quant.err = numeric(0)
quant.corr = matrix(0,0,0)
quant.cov = matrix(0,0,0)

##
## add a vector of measurements and their correlation to the list
## to quant.val, quant.err, quant.corr, quant.cov
##
aeb.meas.add = function(add.val, add.err, add.corr=NULL) {
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

  ##-- assemble correlation matrix
  corr.right = matrix(0, dim(quant.corr)[1], dim(add.corr)[2])
  colnames(corr.right) = colnames(add.corr)
  corr.top = cbind(quant.corr, corr.right)
  corr.left = matrix(0, dim(add.corr)[1], dim(quant.corr)[2])
  rownames(corr.left) = rownames(add.corr)
  corr.bottom = cbind(corr.left, add.corr)
  quant.corr <<- rbind(corr.top, corr.bottom)

  ##-- assemble averaged quantities
  quant.averaged.sel = names(add.val) %in% quant.names.averaged
  quant.val <<- c(quant.val, add.val[quant.averaged.sel])
  quant.err <<- c(quant.err, add.err[quant.averaged.sel])

  ##-- update variance matrix
  quant.cov <<- quant.corr * (quant.err %o% quant.err)
}

##
## add a single additional uncorrelated measurement
## to quant.val, quant.err, quant.corr, quant.cov
##
aeb.meas.add.single = function(label, val, err) {
  names(val) = label
  names(err) = label
  aeb.meas.add(val, err)
}

##
## stores results in alucomb.r log files into global variables
## quant.val, quant.err, quant.corr, quant.cov
##
aeb.collect.data = function(items) {
  for (item in items) {
    if (regexpr("[.]rdata$", item) != -1) {
      file = item
    } else {
      file = file.path(aeb.log.dir.base, item, "average_alucomb.rdata")
    }
    if (!file_test("-f", file)) {
      cat("error: cannot find alucomb log file\n  ", file, "\n")
      next
    }
    ##-- get alucomb results
    rc = suppressWarnings(try(load(file), silent=TRUE))
    rc = attr(rc, "class")
    if (!is.null(rc) && rc == "try-error") {
      cat("error: cannot read alucomb .rdata file\n  ", file, "\n", sep="")
      ## cat("make -C", dir, " update_alucomb\n")
      ## warning("Cannot read alucomb log file\n  ", file.path(aeb.log.dir.base, item, "average_alucomb.log"))
      next
    }

    if (FALSE) {
    ##-- default correlation matrix (NULL if not read)
    corr.current = rc$corr
    conf.lev = pchisq(rc$chisq, df=rc$dof, lower.tail=FALSE)
    if (conf.lev < aeb.cl.min) {
      ##-- inflate errors if CL < CL_min
      rc$err = rc$err * rc$sfact
      ##-- correlation after S-factor error inflation
      corr.current = rc$corr2
    }
    }
    
    aeb.meas.add(quant.val, quant.err, quant.corr)
  }
}

##
## compute linear combination given measurements and their variance
##
aeb.linear.comb = function(lc, val, cov) {
  lc.val = val %*% lc
  lc.var = drop(lc %*% cov %*% lc)
  lc.err = sign(lc.var)*sqrt(abs(lc.var))
  return(c(val=lc.val, err=lc.err))
}

##
## compute linear combinations using measurements stored
## in quant.val, quant.cov
##
aeb.linear.comb.glob = function(lc) {
  if (length(quant.val) == 0) {
    stop("error: no data were loaded")
  }
  diff = setdiff(names(lc), names(quant.val))
  if (length(diff) > 0) {
    stop("error: following quantities were not loaded: ", diff)
  }
  meas.lc = names(lc)
  return(aeb.linear.comb(lc, quant.val[meas.lc], quant.cov[meas.lc, meas.lc]))
}

##
## add quantity, specifying its gradient to compute the covariance
## - add value, error, covariamce, correlation
##
aeb.meas.val.grad.add = function(add.name, add.val, add.grad) {
  names(add.val) = add.name
  add.comb = drop(add.grad)
  add.comb.full = quant.val * 0
  add.comb.full[names(add.comb)] = add.grad
  add.cov = add.comb.full %*% quant.cov %*% add.comb.full
  names(add.cov) = add.name
  quant.cov <<- rbind(
    cbind(quant.cov, matrix(quant.cov %*% add.comb.full, dimnames=list(NULL, add.name))),
    matrix(c(add.comb.full %*% quant.cov, add.cov), 1, dim(quant.cov)[2]+1, dimnames=list(add.name)))
  quant.val <<- c(quant.val, add.val)
  quant.err <<- c(quant.err, sqrt(add.cov))
  quant.corr <<- quant.cov / sqrt(diag(quant.cov)) %o% sqrt(diag(quant.cov))
  return(invisible())
}

##
## add quantity that is combination of other quantities
## - add value, error, covariamce, correlation
##
aeb.meas.comb.add = function(add.name, add.comb) {
  add.comb = drop(add.comb)
  add.comb.full = quant.val * 0
  add.comb.full[names(add.comb)] = add.comb
  add.val = add.comb.full %*% quant.val
  names(add.val) = add.name
  add.err = sqrt(add.comb.full %*% quant.cov %*% add.comb.full)
  names(add.err) = add.name
  quant.val <<- c(quant.val, add.val)
  quant.err <<- c(quant.err, add.err)
  quant.cov <<- rbind(
    cbind(quant.cov, matrix(quant.cov %*% add.comb.full, dimnames=list(NULL, add.name))),
    matrix(c(add.comb.full %*% quant.cov, add.err^2), 1, dim(quant.cov)[2]+1, dimnames=list(add.name)))
  quant.corr <<- quant.cov / sqrt(diag(quant.cov)) %o% sqrt(diag(quant.cov))
  return(invisible())
}

##
## add quantity defined as expression of existing quantities
##
aeb.meas.expr.add = function(add.name, add.expr) {
  rc = mapply(function(name, val) {assign(name, val, pos=1)}, names(quant.val), quant.val)
  add.deriv.expr = deriv(add.expr, all.vars(add.expr))
  add.val = eval(add.deriv.expr)
  add.grad = attr(add.val, "gradient")
  aeb.meas.val.grad.add(add.name, add.val, add.grad)
  return(invisible())
}

##
## given a model matrix that describes how a parameter relates to quantities,
## fit for it and return the fitted combination of quantities
##
aeb.model.matrix.fit = function(model.matrix) {
  model.matrix.names = names(model.matrix)
  val = quant.val[model.matrix.names]
  cov = quant.cov[model.matrix.names, model.matrix.names]
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
aeb.meas.fit.add = function(add.name, add.model.matrix) {
  fit.comb = aeb.model.matrix.fit(add.model.matrix)
  aeb.meas.comb.add(add.name, fit.comb)
  return(invisible())
}
