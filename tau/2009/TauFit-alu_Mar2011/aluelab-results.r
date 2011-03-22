#!/usr/bin/env Rscript

## ////////////////////////////////////////////////////////////////////////////
##
## aluelab-results.r
##
## - elaborate alucomb.r fit results
##
## ////////////////////////////////////////////////////////////////////////////

source("../../../Common/bin/aluelab.r")

args <- commandArgs(TRUE)

##--- option -s, use S-factors
flag.s.factors =  FALSE
if(any(args == "-h")) {
  ##--- help
  cat("aluelab-results.r [flags] [<.rdata file>]\n")
  cat("  no option: traditional Vus but without universality improved Be for B_had\n")
  cat("  -s use error scale factors\n")
  cat("  -u use unitarity constraint to compute Be_univ, B_had_VA, B_had_s\n")
  cat("  -vadirect determine B_VA directly rather than from 1-Be-Bmu-Bstrange\n")
  cat("  -lepuniv use R_had = (1 - (1+(f_mu/f_e))Be_univ)/Be_univ\n")
  args = args[args != "-s"]
  stop()
}
if(any(args == "-s")) {
  ##--- use S-factors
  flag.s.factors =  TRUE
  args = args[args != "-s"]
}
flag.unitarity =  FALSE
if(any(args == "-u")) {
  ##--- use unitarity constraint to compute Be_univ, B_had_VA, B_had_s
  flag.unitarity =  TRUE
  args = args[args != "-u"]
}
flag.vadirect =  FALSE
if(any(args == "-vadirect")) {
  ##--- determine B_VA from direct measurements rather than 1-Be-Bmu-Bstrange
  flag.vadirect =  TRUE
  args = args[args != "-vadirect"]
}
flag.lepuniv =  FALSE
if(any(args == "-lepuniv")) {
  ##--- use R_had = (1 - (1+(f_mu/f_e))Be_univ)/Be_univ
  flag.lepuniv =  TRUE
  args = args[args != "-lepuniv"]
}

if (length(args) > 0) {
  file.name = args[1]
} else {
  file.name = "average_alucomb.rdata"
}

##--- get alucomb results and data
load(file.name)

if (flag.s.factors) {
  quant.err = quant.sf.err
  quant.corr = quant.sf.corr
  quant.cov = quant.sf.cov
}

quant.names = names(quant.val)

if (FALSE) {
if (!any("Gamma801" == names(quant.val))) {
   rc = aeb.meas.expr.add("Gamma801", quote(1.699387*Gamma96))
}
if (!any("Gamma89" == names(quant.val))) {
  rc = aeb.meas.expr.add("Gamma89", quote(Gamma803 + 0.892*Gamma151))
}
}

##
## add measurements to compute universality improved Be
##

##--- from PDG 2009
aeb.meas.add.single("m_e", 0.510998910, 0.000000013)
aeb.meas.add.single("m_mu", 105.658367, 0.000004)
aeb.meas.add.single("tau_tau", 290.6e-15, 1.0e-15)

##--- from HFAG 2009
aeb.meas.add.single("m_tau", 1776.7673082, 0.1507259)

##--- from PDG 2010
aeb.meas.add.single("m_W", 80.399*1e3, 0.023*1e3)
aeb.meas.add.single("tau_mu", 2.197034e-6, 0.000021e-6)

##
## compute phase space factors for Bmu/Be universality
##

##--- phase space factor, function of lepton masses
phspf = quote(1 -8*x + 8*x^3 - x^4 - 12*x^2*log(x))
##--- phase space factors for e/tau, mu/tau, e/mu
rc = aeb.meas.expr.add("phspf_mebymtau",  eval(bquote(substitute(.(phspf), list(x=quote(m_e^2/m_tau^2))))))
rc = aeb.meas.expr.add("phspf_mmubymtau", eval(bquote(substitute(.(phspf), list(x=quote(m_mu^2/m_tau^2))))))
rc = aeb.meas.expr.add("phspf_mebymmu", eval(bquote(substitute(.(phspf), list(x=quote(m_e^2/m_mu^2))))))
rc = aeb.meas.expr.add("Bmu_by_Be_th", quote(phspf_mmubymtau/phspf_mebymtau))

aeb.meas.expr.add("Be", quote(Gamma5))
aeb.meas.expr.add("Bmu", quote(Gamma3))

##--- Be from Bmu
aeb.meas.expr.add("Be_from_Bmu", quote(phspf_mebymtau/phspf_mmubymtau *Bmu))

##
## rad. corrections from to get Be from tau lifetime
## values from 10.1103/RevModPhys.78.1043 p.1047, arXiv:hep-ph/0507078v2 p.7, could be recomputed
## - delta^L_gamma = 1 + alpha(mL)/2pi * (25/4 - pi^2)
## - delta^L_W = 1 + 3/5* m_L^2/M_W^2
##
delta.tau.gamma = 1 - 43.2e-4
delta.mu.gamma = 1 - 42.4e-4
delta.tau.W = 1 + 2.9e-4
delta.mu.W = 1 + 1.0e-6

##
## Be from tau lifetime
## Be = tau_tau / tau_mu (m_tau/m_mu)^5 f(m^2_e/m^2_tau)/f(m^2_e/m^2_mu) (delta^tau_gamma delta^tau_W)/(delta^mu_gamma delta^mu_W)
##
aeb.meas.expr.add("Be_from_taulife",
                  bquote(tau_tau/tau_mu * (m_tau/m_mu)^5 * phspf_mebymtau/phspf_mebymmu
                         *.(delta.tau.gamma) *.(delta.tau.W) /.(delta.mu.gamma) /.(delta.mu.W)))
##
## Bmu from tau lifetime
## Bmu = tau_tau / tau_mu (m_tau/m_mu)^5 f(m^2_mu/m^2_tau)/f(m^2_e/m^2_mu) (delta^tau_gamma delta^tau_W)/(delta^mu_gamma delta^mu_W)
##
aeb.meas.expr.add("Bmu_from_taulife",
                  bquote(tau_tau/tau_mu * (m_tau/m_mu)^5 * phspf_mmubymtau/phspf_mebymmu
                         *.(delta.tau.gamma) *.(delta.tau.W) /.(delta.mu.gamma) /.(delta.mu.W)))

##
## universality improved Be = B(tau -> e nu nubar (gamma))
## see arXiv:hep-ph/0507078v2 p.7, doi:10.1103/RevModPhys.78.1043 p.1047
##
## minimum chisq fit for Be_univ using, Be, Be from Bmu, Be from tau lifetime
##
## Bmu/Be = f(m_mu^2/m_tau^2) / f(m_e^2/m_tau^2)
## Be = tau_tau / tau_mu (m_tau/m_mu)^5 f(m^2_e/m^2_tau)/f(m^2_e/m^2_mu) (delta^tau_gamma delta^tau_W)/(delta^mu_gamma delta^mu_W)
##
aeb.meas.fit.add("Be_univ", c(Be=1, Be_from_Bmu=1, Be_from_taulife=1))

display.names = c(
  "Be", "Bmu", "Be_from_Bmu", "Be_from_taulife", "Be_univ",
  "Gamma9", "Gamma10", "Gamma66"
  )
show(rbind(cbind(val=quant.val[display.names], err=quant.err[display.names])))

if (FALSE) {
  ##--- print couplings ratio
  aeb.meas.expr.add("B_mu/B_e", quote(Gamma3/Gamma5))
  aeb.meas.expr.add("g_mu/g_e", quote(sqrt(Gamma3/Gamma5/Bmu_by_Be_th)))
  display.names = c(
    "Gamma5", 
    "Gamma3",
    "B_mu/B_e",
    "g_mu/g_e"
    )
  show(rbind(cbind(val=quant.val[display.names], err=quant.err[display.names])))
}
