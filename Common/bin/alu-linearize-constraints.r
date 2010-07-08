#!/usr/bin/env Rscript

##
## alu-linear-constr.r
##
## - get fitted quantities from alucomb .rdata file
## - read constraint equations from constraints.txt
## - compute symbolic gradient vs. constraint variables
## - compute numerical linearized constraint
##
## c = f(x_i)
## c = f(x_i)|0 + Sum( ( x_i - x_i|0 ) * df/dx_i|0 )
## c - f(x_i)|0 +  Sum( x_i|0 * df/dx_i|0 ) = Sum( x_i* df/dx_i|0 )
##

load("average_alucomb.rdata")

lines  <- readLines("constraints.txt")
lines = lines[grep("^[*]\\s*\\S+\\s*=\\s*.*\\S\\s*$", lines, perl=TRUE)]
lines = gsub("^[*]\\s*", "", lines, perl=TRUE)
lines = gsub("\\s+$", "", lines, perl=TRUE)
lines.text = lines
lines = gsub("(\\S+)\\s*=\\s*(.*)", "\\1;\\2", lines, perl=TRUE)
constrs = strsplit(lines, ";")

env = new.env()
mapply(function(name, value) assign(name, value, env=env), names(quant), quant)

rc = lapply(constrs, function(constr) {
  constr.v = constr[[1]]
  constr.expr = parse(text=constr[[2]])
  constr.0 = as.numeric(eval(constr.expr, env=env))

  vars = all.vars(constr.expr, unique=TRUE)
  grad.expr = deriv(constr.expr, vars)
  grad.0 = as.numeric(eval(grad.expr, env=env))

  ##--- get symbolic expression for linearized constraint
  grad = eval(grad.expr, env=env)
  grad = attr(grad, "gradient")
  grad.text = paste(as.vector(grad), colnames(grad), sep="*", collapse=" + ")

  rc = paste(constr.v, " + ", grad.0 - constr.0, " = ", grad.text, sep="")
})

rc = mapply(function(constr, linear) {
  cat("constraint = ", constr, "\n", sep="")
  cat("linearized = ", linear, "\n", sep="")
}, lines.text, rc)
