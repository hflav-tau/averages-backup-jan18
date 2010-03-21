#!/usr/bin/env Rscript

## ////////////////////////////////////////////////////////////////////////////
##
## plot measurements in pdg_average.input
##

## ////////////////////////////////////////////////////////////////////////////
## definitions

##--- trim leading and trailing spaces in strings
trim <- function(char) {
  return(sub("\\s+$", "", sub("^\\s+", "", char)))
}

##--- trim trailing spaces in strings
trim.trailing <-   function(char) {
  return(sub("\\s+$", "", char))
}

## ////////////////////////////////////////////////////////////////////////////
## code

fh <- file("pdg_average.input")
lines  <- readLines(fh)
close(fh)

measurement = ""

lines <- trim.trailing(lines)
meas.label = gsub("^[*c#;]\\s*", "", lines[1])

measurements = NULL
for (line in lines[-1]) {
  if (regexpr("^\\s*$",line) != -1 ||
      regexpr("^[*cC#;]",line) != -1) {
    next
  }
  temp = gsub("\\","\\b",line,fixed=TRUE)
  temp = gsub(";","\\s",temp,fixed=TRUE)
  fields = strsplit(gsub("\\s*(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(.*)","\\1;\\2;\\3;\\4", temp, perl=TRUE),";")
  fields = fields[[1]]
  fields = gsub("\\s",";",fields,fixed=TRUE)
  fields = gsub("\\b","\\",fields,fixed=TRUE)
  num.fields = fields[1:3]
  warning <- FALSE
  withCallingHandlers({meas <- as.numeric(num.fields)}, warning = function(w) {warning <<- TRUE; invokeRestart("muffleWarning")})
  if (warning) {
    cat("Cannot interpret the following line as a measurement\n")
    cat("  '",line,"'\n")
    print(num.fields)
  } else {
    measurements <- rbind(measurements,data.frame(val=meas[1], stat=meas[2], syst=meas[3], bibitem=fields[4])) 
  }
}
measurements

