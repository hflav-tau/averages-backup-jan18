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
  file.name = "average.rdata"
}

##--- get alucomb results and data
load(file.name)

if (flag.s.factors) {
  quant.err = quant.sf.err
  quant.corr = quant.sf.corr
  quant.cov = quant.sf.cov
}

quant.names = names(quant.val)

##
## recover some BRs as function of others
##
if (!any("Gamma801" == names(quant.val))) {
   rc = aeb.meas.expr.add("Gamma801", quote(1.699387*Gamma96))
}
if (!any("Gamma89" == names(quant.val))) {
  rc = aeb.meas.expr.add("Gamma89", quote(Gamma803 + 0.892*Gamma151))
}

##
## Gamma151 = K omega nu
## PDG 2009 ( 4.1 ± 0.9 ) × 10 ~ 4
##+++ do not include, included in K 3pi nu
##+++ aeb.meas.add.single("Gamma151", 4.1e-4, 0.9e-4)

##
## B[tau- -> (K3pi)- nu (ex. K0, omega,eta)]  : 1e-2 : 0.074 ± 0.030    
## Aleph estimate, not used
##
##+++ aeb.meas.add.single("Km3PiNu", 0.074e-2, 0.030e-2)

##
## B[tau- -> (K4pi)- nu]                      : 1e-2 : 0.011 ± 0.007
## Aleph estimate, not used
##
##+++ aeb.meas.add.single("Km4PiNu", 0.011e-2, 0.007e-2)

##
## PDG 2009 definition of Gamma110 = B(tau -> Xs nu)
##
Gamma110_pdg09.comb = c(
  Gamma10=1,  Gamma16=1, 
  Gamma23=1,  Gamma28=1,  Gamma35=1, 
  Gamma40=1,  Gamma85=1,  Gamma89=1,
  Gamma128=1)
aeb.meas.comb.add("Gamma110_pdg09", Gamma110_pdg09.comb)

##
## Gamma110 = B(tau -> Xs) -- local definition
##
## Gamma110 =
##   Gamma10+Gamma16+Gamma23+Gamma28+Gamma35+Gamma40+Gamma44+
##   Gamma53+Gamma85+Gamma89+Gamma128+Gamma130+Gamma132+
##   Gamma96+Gamma96*B(phi -> KS0 KL0)/B(phi -> K+ K-)
##
Gamma110.comb = c(
  Gamma10=1,  Gamma16=1, 
  Gamma23=1,  Gamma28=1,  Gamma35=1, 
  Gamma40=1,  Gamma44=1,  Gamma53=1,
  Gamma85=1,  Gamma89=1,  Gamma128=1,
  Gamma130=1, Gamma132=1,
  Gamma96=1.699387)

##
## Oct 2010, Gamma110 updated as follows
##
## Gamma110 = Gamma10  + Gamma16   + Gamma23   + Gamma28  + Gamma35  + Gamma40 + Gamma85 + Gamma89 + Gamma128
##          + Gamma151 + Gamma130  + Gamma132  + Gamma44  + Gamma53  + Gamma801
## 
## - replaced Gamma96 with Gamma801 = Gamma96 * [1 + B(phi -> KS0 KL0)/B(phi -> K+ K-)]
##   because of technical issues with readpdg.cc
## - added Gamma151 (tau -> K omega nu)
##
Gamma110.comb = c(
  Gamma10=1,  Gamma16=1,  Gamma23=1,  
  Gamma28=1,  Gamma35=1,  Gamma40=1,  
  Gamma44=1,  Gamma53=1,  Gamma85=1,  
  Gamma89=1,  Gamma128=1, Gamma130=1, 
  Gamma132=1, Gamma151=1, Gamma801=1)

if (!any("Gamma151" == quant.names)) {
  ##--- remove Gamma151 (tau -> K omega nu) if not present in fit for backward compatibility
  Gamma110.comb = Gamma110.comb[names(Gamma110.comb) != "Gamma151"]
}

##
## end Oct2010 update
## define Gamma110 using the alucomb.r Gamma110 constraint
## this is more robust against changes in the fitting 
##

Gamma110.comb = combination$constr.lin.comb[["Gamma110.coq"]]
Gamma110.comb = Gamma110.comb[names(Gamma110.comb) != "Gamma110"]
Gamma110.names = names(Gamma110.comb[Gamma110.comb != 0])

##--- use Gamma110 defined through COMBOFQUANT rather than here
## aeb.meas.comb.add("Gamma110", Gamma110.comb)

##--- list of all tau BRs that are not leptonic and not strange, i.e. not-strange-hadronic
B.tau.VA.names = setdiff(names(combination$constr.lin.comb[["GammaAll"]]), names(combination$constr.lin.comb[["Gamma110.coq"]]))
B.tau.VA.names = setdiff(B.tau.VA.names, c("Gamma5", "Gamma3", "Gamma998"))
aeb.meas.expr.add("B_tau_VA", parse(text=paste(B.tau.VA.names, collapse="+")))
aeb.meas.expr.add("B_tau_VA_unitarity", quote(1-Gamma5-Gamma3-Gamma110))
aeb.meas.expr.add("B_tau_s_unitarity", quote(1-Gamma5-Gamma3-B_tau_VA))

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
## Be, from unitarity = 1 - Bmu - B_VA - B_s
## Bmu, from unitarity = 1 - Be - B_VA - B_s
##
aeb.meas.expr.add("Be_unitarity", quote(1 - Gamma3 - B_tau_VA - Gamma110))
aeb.meas.expr.add("Bmu_unitarity", quote(1 - Gamma5 - B_tau_VA - Gamma110))

##
## fit best values for Be, Bmu, B_tau_VA, B_tau_s
## using both the direct measurements and the result from the unitarity constraint
##
if (any("Gamma998" == names(combination$constr.lin.comb[["GammaAll"]]))) {
  ##--- here we got results from unconstrained fit
  if (flag.unitarity) {
    ##--- use unitarity constraint to compute Be, Bmu, B_tau_VA/s
    aeb.meas.fit.add("Be_fit", c(Gamma5=1, Be_unitarity=1))
    aeb.meas.fit.add("Bmu_fit", c(Gamma3=1, Bmu_unitarity=1))
    aeb.meas.fit.add("B_tau_VA_fit", c(B_tau_VA=1, B_tau_VA_unitarity=1))
    aeb.meas.fit.add("B_tau_s_fit", c(Gamma110=1, B_tau_s_unitarity=1))
  } else {
    ##--- direct measurements for leptonic BRs and strange hadronic BR
    aeb.meas.expr.add("Be_fit", quote(Gamma5))
    aeb.meas.expr.add("Bmu_fit", quote(Gamma3))
    aeb.meas.expr.add("B_tau_s_fit", quote(Gamma110))
    if (flag.vadirect) {
      ##--- direct measurement for non-strange hadronic BR
      aeb.meas.expr.add("B_tau_VA_fit", quote(B_tau_VA))
    } else {
      ##
      ## non-strange hadronic BR by subtraction, however note that
      ## here we do not use the improved Be from universality
      ## (a Be fit using Be, Bmu and tau lifetime) as in
      ## M.Davier et al, RevModPhys.78.1043, arXiv:hep-ph/0507078
      ##
      aeb.meas.expr.add("B_tau_VA_fit", quote(1-Gamma5-Gamma3-Gamma110))
    }
  }
} else {
  ##--- here we got results from unitarity constrained fit
  rc = aeb.meas.expr.add("Be_fit", quote(Gamma5))
  rc = aeb.meas.expr.add("Bmu_fit", quote(Gamma3))
  rc = aeb.meas.expr.add("B_tau_VA_fit", quote(B_tau_VA))
  rc = aeb.meas.expr.add("B_tau_s_fit", quote(Gamma110))
}

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

##--- Be from Bmu
aeb.meas.expr.add("Be_from_Bmu", quote(phspf_mebymtau/phspf_mmubymtau *Bmu_fit))

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
aeb.meas.fit.add("Be_univ", c(Be_fit=1, Be_from_Bmu=1, Be_from_taulife=1))

##
## Vud
##
## arXiv:0710.3181v1 [nucl-th], 10.1103/PhysRevC.77.025501
## I.S.Towner, J.C.Hardy, An improved calculation of the isospin-symmetry-breaking corrections to superallowed Fermi beta decay
##
Vud.val = 0.97425
Vud.err = 0.00022
aeb.meas.add.single("Vud", Vud.val, Vud.err)

##
## SU3 breaking correction, straight from papers
##

##--- POS(KAON)08, A.Pich, Theoretical progress on the Vus determination from tau decays
deltaR.su3break.val = 0.216
deltaR.su3break.err = 0.016
##--- E. Gamiz et al., Nucl.Phys.Proc.Suppl.169:85-89,2007, arXiv:hep-ph/0612154v1
deltaR.su3break.val = 0.240
deltaR.su3break.err = 0.032
## aeb.meas.add.single("deltaR_su3break", deltaR.su3break.val, deltaR.su3break.err)

##
## SU3 breaking correction, recompute from data following
## E. Gamiz et al., Nucl.Phys.Proc.Suppl.169:85-89,2007, arXiv:hep-ph/0612154v1
##

##--- s quark mass, PhysRevD.74.074009
aeb.meas.add.single("m_s", 94, 6)

##--- E.Gamiz, M.Jamin, A.Pich, J.Prades, F.Schwab, |V_us| and m_s from hadronic tau decays
aeb.meas.add.single("deltaR_su3break_pheno", 0.1544, 0.0037)
aeb.meas.add.single("deltaR_su3break_msd2", 9.3, 3.4)
aeb.meas.add.single("deltaR_su3break_remain", 0.0034, 0.0028)
aeb.meas.expr.add("deltaR_su3break", quote(deltaR_su3break_pheno + deltaR_su3break_msd2*(m_s/1000)^2 + deltaR_su3break_remain))

if (any("Gamma998" == names(combination$constr.lin.comb[["GammaAll"]]))) {
  ##
  ## if using constrained fit with dummy mode (Gamma998), i.e. unconstrained fit
  ##
  if (flag.lepuniv) {
    ##--- determine R_tau using leptonic BRs and universality
    aeb.meas.expr.add("R_tau", quote(1/Be_univ -1 -phspf_mmubymtau/phspf_mebymtau))
    aeb.meas.expr.add("R_tau_s", quote(B_tau_s_fit/Be_univ))
    aeb.meas.expr.add("R_tau_VA", quote(R_tau - R_tau_s))
  } else {
    ##
    ## use "fit" values computed in this script (possibly incorporating unitarity constraint)
    ##
    ##--- add R_tau_VA = R - R_tau_s
    aeb.meas.expr.add("R_tau_VA", quote(B_tau_VA_fit/Be_univ))
    ##--- add R_tau_s = B(tau -> Xs nu) / Be_univ
    aeb.meas.expr.add("R_tau_s", quote(B_tau_s_fit/Be_univ))
    ##--- add R_tau as function of quantities
    aeb.meas.expr.add("R_tau", quote(R_tau_VA+R_tau_s))
  }
} else {
  ##
  ## if using constrained fit without dummy mode (Gamma998), i.e. constrained fit
  ##
  ##--- add R_tau as function of quantities
  if (flag.lepuniv) {
    ##--- determine R_tau using leptonic BRs and universality
    aeb.meas.expr.add("R_tau", quote(1/Be_univ -1 -phspf_mmubymtau/phspf_mebymtau))
  } else {
    ##--- use values computed in the alucomb.r fit, incorporating unitarity constraints
    aeb.meas.expr.add("R_tau", quote((B_tau_VA+Gamma110)/Be_univ))
  }
  ##--- add R_tau_s = B(tau -> Xs nu) / Be_univ
  aeb.meas.expr.add("R_tau_s", quote(Gamma110/Be_univ))
  ##--- add R_tau_VA = R_tau - R_tau_s
  aeb.meas.expr.add("R_tau_VA", quote(R_tau - R_tau_s))
}

##--- add Vus
aeb.meas.expr.add("Vus", quote(sqrt(R_tau_s/(R_tau_VA/Vud^2 - deltaR_su3break))))

##--- theory error
Vus.err.th = abs(quant.cov["Vus", "deltaR_su3break"])/quant.err["deltaR_su3break"]
Vus.err.exp = sqrt(quant.err["Vus"]^2 - Vus.err.th^2)

##-- using Hardy-Towner 2009
Vus.unitarity.val = 0.2255
Vus.unitarity.err = 0.0010
nsigma = (quant.val["Vus"] - Vus.unitarity.val)/quadrature(c(quant.err["Vus"], Vus.unitarity.err))

display.names = c(Gamma110.names,
  "Gamma5", "Be_unitarity", "Be_fit",
  "Gamma3", "Bmu_unitarity", "Bmu_fit",
  "Bmu_by_Be_th", "Be_from_Bmu", "Be_from_taulife", "Be_univ", "Bmu_from_taulife",
  "B_tau_VA", "B_tau_VA_unitarity", "B_tau_VA_fit",
  "Gamma110", "B_tau_s_unitarity", "B_tau_s_fit", "Gamma110_pdg09",
  "R_tau", "R_tau_s", "R_tau_VA", "deltaR_su3break", "Vus")
show(rbind(cbind(val=quant.val[display.names], err=quant.err[display.names]),
           Vus_err_perc = c(quant.err["Vus"]/quant.val["Vus"]*100, 0),
           Vus_err_exp = c(Vus.err.exp, 0),
           Vus_err_th = c(Vus.err.th, 0),
           Vus_err_th_perc = c(Vus.err.th /quant.val["Vus"]*100, 0),
           nsigma = c(val=nsigma, err=0)
           ))

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
