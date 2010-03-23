#!/usr/bin/env Rscript

## ////////////////////////////////////////////////////////////////////////////
##
## plot measurements in pdg_average.input
##

library(plotrix)

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
    measurements =
      rbind(measurements,data.frame(val=meas[1],
                                    stat=meas[2],
                                    syst=meas[3],
                                    bibitem=fields[4]))
  }
}

meas = measurements$val
error = sqrt(measurements$stat^2 + measurements$syst^2)

x.min = min(meas-error)
x.max = max(meas+error)
x.range = x.max - x.min
x.min = x.min-0.2*x.range
x.max = x.max+0.2*x.range
ratio.labels = 1/2
x.max.plot = x.min + (x.max-x.min)/(1-ratio.labels)

y.num = length(meas)
y.tot = y.num+2
y.min = (y.tot - y.num +1)/2 + 0.5
y.max = (y.tot + y.num -1)/2 + 0.5

postscript("average.eps")
line.height = par("cin")[2]
dev.height = line.height*(par("mar")[1] + par("mar")[3] + y.num+2)
rc = dev.off()

postscript("average.eps", onefile=FALSE, width=8, height=dev.height,
           paper="special", horizontal=FALSE, title=meas.label)

plotCI(meas, y.min:y.max, uiw=error, err="x",
  xlim=c(x.min, x.max.plot), ylim=c(1,y.tot),
  xlab="", ylab="", yaxt = "n", pch=16, lwd=1.5, cex=0.7,
  main=meas.label)

text(x.max, y.min:y.max, pos=4, labels=measurements$bibitem, cex=0.9)
rc = dev.off()
