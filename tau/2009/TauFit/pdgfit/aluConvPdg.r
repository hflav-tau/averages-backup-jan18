#!/usr/bin/env Rscript

## ////////////////////////////////////////////////////////////////////////////
## definitions

##
## return list of lines from file
##
get.file.lines <- function(fname) {
  lines  <- readLines(fname)
  return(lines)
}

##
## test if pattern matches string irrespective of letters case
##
match.nocase = function(pattern, str) {
  return(regexpr(pattern, str, ignore.case=TRUE) != -1)
}
 
## ////////////////////////////////////////////////////////////////////////////
## code

file = "s035-fit-no-babar-belle-no-scaling-with-dummy-mode.fit"
args <- commandArgs(TRUE)
if (length(args) > 0) file = args[1]

lines = get.file.lines(file)
cat("read file", file, "\n", file=stderr())

par.beg = grep("^\\s*RESULTS FOR PARAMETERS", lines, perl=TRUE)[3] + 3
par.end = grep("^\\s*$", tail(lines, -par.beg), perl=TRUE)[1] + par.beg -1

params = lapply(strsplit(lines[par.beg:par.end], "\\s+", perl=TRUE),
  function(x) {c(x[c(2,5,6)], gsub(".*tau- --", "tau- --", paste(tail(x,-8), collapse=" ")))})

corr.beg = grep("^\\s*FITTED PARTIAL DECAY MODE BRANCHING", tail(lines, -par.end), perl=TRUE)[1] + par.end + 2
params.num = length(params)

corr = matrix(0, nrow=params.num, ncol=params.num)
col.beg = 1
line.beg = corr.beg
line.end = corr.beg+params.num
params.done = 0
params.next = params.num
repeat {
  par.per.line = length(unlist(gregexpr("P\\s*\\d+", lines[line.beg], perl=TRUE)))
  tmp = lapply(strsplit(gsub("P\\s*(\\d+)", "P\\1", lines[(line.beg+1):line.end], perl=TRUE),
    "\\s+", perl=TRUE), function(x) type.convert(tail(x,-2)))
  for (i in 1:params.next) {
    for (j in 1:par.per.line) {
      row = params.done + i
      col = params.done + j
      if (row == col) {
        corr[row, col] = 1
        corr[col, row] = 1
      } else {
        corr[row, col] = tmp[[i]][j] / 100
        corr[col, row] = tmp[[i]][j] / 100
      }
    }
  }
  params.done = params.done + par.per.line
  if (params.done >= params.num) break
  params.next = params.num - params.done
  line.beg = line.end+2
  line.end = line.beg + params.next
}

##      1 tau- --> mu- nubar(mu) nu(tau)
##      2 tau- --> eta K- nu(tau)
##      3 tau- --> h- 4pi0 nu(tau) (ex.K0,eta)
##      4 tau- --> h- omega pi0 nu(tau)
##      5 tau- --> K- 2pi0 nu(tau) (ex.K0)
##      6 tau- --> K- 3pi0 nu(tau) (ex.K0, eta)
##      7 tau- --> pi- Kbar0 nu(tau)
##      8 tau- --> pi- Kbar0 pi0 nu(tau)
##      9 tau- --> K- K0 pi0 nu(tau)
##     10 tau- --> pi- nu(tau)
##     11 tau- --> pi- pi0 nu(tau)
##     12 tau- --> K- pi0 nu(tau)
##     13 tau- --> dummy mode used  by the fit
##     14 tau- --> e- nubar(e) nu(tau)
##     15 tau- --> pi- 2pi0 nu(tau) (ex.K0)
##     16 tau- --> pi- 3pi0 nu(tau) (ex.K0)
##     17 tau- --> h- h- h+ 3pi0 nu(tau)
##     18 tau- --> pi- K(S)0 K(S)0 nu(tau)
##     19 tau- --> h- h- h+ 2pi0 nu(tau) (ex.K0,omega,eta)
##     20 tau- --> pi- K(S)0 K(L)0 nu(tau)
##     21 tau- --> K- K+ pi- pi0 nu(tau)
##     22 tau- --> pi- pi+ pi- nu(tau) (ex.K0,omega)
##     23 tau- --> K- pi+ pi- nu(tau) (ex.K0)
##     24 tau- --> pi- pi+ pi- pi0 nu(tau) (ex.K0,omega)
##     25 tau- --> K- pi+ pi- pi0 nu(tau) (ex.K0,eta)
##     26 tau- --> 3h- 2h+ nu(tau) (ex.K0)
##     27 tau- --> 3h- 2h+ pi0 nu(tau) (ex.K0)
##     28 tau- --> K- K+ pi- nu(tau)
##     29 tau- --> eta pi- pi0 nu(tau)
##     30 tau- --> K- K0 nu(tau)
##     31 tau- --> K- nu(tau)
##     32 tau- --> h- omega nu(tau)

if (FALSE) {
cat("* /////////////////////////////////////////////////////////////////////////////
* Gamma3   tau- --> mu- nubar_mu nu_tau
* Gamma5   tau- --> e- nubar_e nu_tau
* Gamma9   tau- --> pi- nu_tau
* Gamma10  tau- --> K- nu_tau
* Gamma14  tau- --> pi- pi0 nu_tau
* Gamma16  tau- --> K- pi0 nu_tau
* Gamma20  tau- --> pi- 2pi0 nu_tau (ex.K0)
* Gamma23  tau- --> K- 2pi0 nu_tau (ex.K0)
* Gamma27  tau- --> pi- 3pi0 nu_tau (ex.K0)
* Gamma28  tau- --> K- 3pi0 nu_tau (ex.K0,eta)
* Gamma30  tau- --> h- 4pi0 nu_tau (ex.K0,eta)
* Gamma35  tau- --> pi- Kbar0 nu_tau
* Gamma37  tau- --> K- K0 nu_tau
* Gamma40  tau- --> pi- Kbar0 pi0 nu_tau
* Gamma42  tau- --> K- K0 pi0 nu_tau
* Gamma47  tau- --> pi- K0S K0S nu_tau
* Gamma48  tau- --> pi- K0S K0L nu_tau
* Gamma62  tau- --> pi- pi+ pi- nu_tau (ex.K0,omega)
* Gamma70  tau- --> pi- pi+ pi- pi0 nu_tau  (ex.K0,omega)
* Gamma77  tau- --> h- h- h+ 2pi0 nu_tau  (ex.K0,omega,eta)
* Gamma78  tau- --> h- h- h+ 3pi0 nu_tau
* Gamma85  tau- --> K- pi+ pi- nu_tau (ex.K0)
* Gamma89  tau- --> K- pi+ pi- pi0 nu_tau (ex.K0,eta)
* Gamma93  tau- --> K- K+ pi- nu_tau
* Gamma94  tau- --> K- K+ pi- pi0 nu_tau
* Gamma103 tau- --> 3h- 2h+ nu_tau (ex.K0)
* Gamma104 tau- --> 3h- 2h+ pi0 nu_tau (ex.K0)
* Gamma126 tau- --> eta pi- pi0 nu_tau
* Gamma128 tau- --> eta K- nu_tau
* Gamma150 tau- --> h- omega nu_tau
* Gamma152 tau- --> h- omega pi0 nu_tau
* /////////////////////////////////////////////////////////////////////////////
")
}

##-- convert to "info file" numbering
conv.gamma = c(3, 128, 30, 152, 23, 28, 35, 40, 42, 9, 14, 16, 999, 5, 20,
  27, 78, 47, 77, 48, 94, 62, 85, 70, 89, 103, 104, 93, 126, 37, 10, 150)

##-- adjust conversion for no-dummy file
if (params.num == 31) {
  conv.gamma = c(conv.gamma[c(1:12, 14:32)])
}

conv.gamma.r = lapply(1:length(conv.gamma), function(x) x)
names(conv.gamma.r) = as.character(conv.gamma)
order.gamma = unlist(conv.gamma.r[as.character(conv.gamma[order(conv.gamma)])])

cat("* /////////////////////////////////////////////////////////////////////////////\n")
for (param in params[order.gamma]) {
  num = as.numeric(param[1])
  num.gamma = conv.gamma[num]
  cat(sprintf("* Gamma%-3d %s\n", num.gamma, param[4]))
}
cat("* /////////////////////////////////////////////////////////////////////////////\n")

for (param in params[order.gamma]) {
  cat("\n*", param[4], "\n")
  num = as.numeric(param[1])
  num.gamma = conv.gamma[num]
  cat("BEGIN       PDG Gamma", num.gamma, " published PDG09\n", sep="")
  cat("MEASUREMENT m_Gamma", num.gamma, " statistical systematic\n", sep="")
  cat("DATA        m_Gamma", num.gamma, " statistical systematic\n", sep="")
  cat("            ", param[2], " ", param[3], " ", 0, "\n", sep="")
  for (i in order.gamma) {
    if (i != num && corr[num, i] != 0) {
      cat("STAT_CORR_WITH PDG Gamma", conv.gamma[i], " published ", corr[num, i], "\n", sep="")
    }
  }
  cat("END", "\n")
}
