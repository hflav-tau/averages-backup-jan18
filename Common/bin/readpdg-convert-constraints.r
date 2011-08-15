#!/usr/bin/env Rscript

## ////////////////////////////////////////////////////////////////////////////
##
## readpdg-convert-constraints.r
##
## - read .input file containing constraints comments produced
##   by readpdg.cc
## - convert to alucomb2.r style linear or non linear constraints
##
## usage:
## - readpdg-convert-constraints.r <file.input>
##
## ////////////////////////////////////////////////////////////////////////////

library(methods)
library(stringr)

args <- commandArgs(TRUE)

##--- option -s, use S-factors
flag.s.factors =  FALSE
if(any(args == "-h")) {
  ##--- help
  cat("constr-convert.r <input file>\n")
  args = args[args != "-h"]
  stop()
}

if (length(args) <= 0) {
  stop("please give an input file as argument")
}
file.name = args[1]

if (!file.exists(file.name)) {
  stop("cannot find file ", file.name, "\n")
}

lines = readLines(file.name)
cat("read file", file.name, "\n")

for (line in lines) {
  ##
  ## look for patterns like
  ## * Gamma150by66 = ( 1.0*Gamma800 + 1.0*Gamma151) / ( 0.2274*Gamma128 + 0.0153*Gamma152 )
  ##
  rc = str_match(line, "^[*]\\s+([[:alnum:]]+\\s+=\\s+.*\\S+)\\s*$")[2]
  if (!is.na(rc)) {
    rc = str_replace_all(rc, "\\s+", " ")
    rc = str_replace_all(rc, "\\(\\s+", "(")
    rc = str_replace_all(rc, "[.]?0+([*])", "\\1")
    rc = str_replace_all(rc, "(\\D)1[*]", "\\1")
    rc = str_replace_all(rc, "\\(([^*/+-]+)\\)", "\\1")
    parts = str_match(rc, "(.*\\S)\\s*=\\s*(\\S.*)")
    nlconstr.str =  paste("-", parts[2], " + ", parts[3], sep="")
    nlconstr.expr = parse(text=nlconstr.str)
    linear.flag = TRUE
    coeffs = list()
    vars = all.vars(nlconstr.expr)
    for (var in vars) {
      nlconstr.deriv = D(nlconstr.expr, var)
      if (length(all.vars(nlconstr.deriv)) > 0) {
        linear.flag = FALSE
      } else {
        coeffs[[var]] = eval(nlconstr.deriv)
      }
    }
    if (linear.flag) {
      cat("COMBOFQUANT", vars[1], paste(coeffs[vars[-1]], vars[-1]), "\n")
    } else {
      cat("NLCONSTRAINT", paste(parts[2],"c",sep="."), "0", paste("\"", nlconstr.str, "\"\n", sep=""))
    }
  }
}
