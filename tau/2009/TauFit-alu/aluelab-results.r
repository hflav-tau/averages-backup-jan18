#!/usr/bin/env Rscript

## ////////////////////////////////////////////////////////////////////////////
##
## aluelab-results.r
##
## - elaborate alucomb.r fit results
##
## ////////////////////////////////////////////////////////////////////////////

##
## Sudan's paper
## B(tau -> KSKL pi nu) = 0.00087 +- 0.000085 +- 0.000030
##+++ include if approved in time
##

source("../../../Common/bin/aluelab.r")

args <- commandArgs(TRUE)
if (length(args) > 0) {
  file.name = args[1]
} else {
  file.name = "average_alucomb.rdata"
}

##--- get alucomb results and data
load(file.name)

quant.val = quant
rm(quant)
## quant.err = quant3.err
## quant.corr = quant3.corr
## quant.cov = quant3.cov

meas.val = meas
rm(meas)

##
## Gamma43 = pi K0 >=1pi0 nu
## PDG 2009 -- 0.324 ± 0.074 ± 0.066 ABBIENDI 00C OPAL
##
##+++ aeb.meas.add.single("Gamma43", 0.324e-2, quadrature(0.074e-2, 0.0066e-2))

##
## Gamma44 = pi K0 pi0 pi0 nu
## PDG 2009 -- 0.26 ± 0.24 BARATE 99R ALEPH
##+++ aeb.meas.add.single("Gamma44", 2.6e-4, 2.4e-4)

##
## Gamma53 =  K0 h+ h- h- nu
## PDG 2009 ( ( 2.3 ± 2.0 ) × 10~4
## PDG 2009 -- (2.3 ± 1.9 ± 0.7) × 10-4	BARATE	 98E ALEPH
##+++ aeb.meas.add.single("Gamma53", 2.3e-4, 2.0e-4)

##
## Gamma144 = K phi nu
## KmPhiNu HFAG 2009 3.704330e-05 +- 3.259888e-06 (S=1.3)
##
##+++ aeb.meas.add.single("Gamma144", 3.704330e-05, 3.259888e-06)

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

##--- measurement names matching a specific PDG Gamma
meas.match = function(gamma) {
  return( regexpr(paste("[.]", gamma, "[.]", sep=""), names(meas.val)) != -1 )
}

##--- sum in quadrature
quadrature = function(x) {
  return(sqrt(sum(x^2)))
}

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
## compute Gamma110 = B(tau -> Xs)
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
Gamma110.names = names(Gamma110.comb[Gamma110.comb != 0])
## aeb.meas.comb.add("Gamma110", Gamma110.comb)

##
## B(tau -> Xs nu)
##
## can use Gamma110 via COMBOFQUANT in alucomb.r
## same definition as above
##
B_tau_s.val = quant.val["Gamma110"]
B_tau_s.err = quant.err["Gamma110"]

##
## universality Be = B(tau -> e nu nubar)
##
## minimum chisq fit constraining Bmu/Be from Standard Model, using Be, Bmu
## tau lifetime should also be included, see arXiv:hep-ph/0507078v2 p.89
##

##--- from PDG 2009
aeb.meas.add.single("m_e", 0.510998910, 0.000000013)
aeb.meas.add.single("m_mu", 105.658367, 0.000004)
aeb.meas.add.single("tau_tau", 290.6, 1.0)

##--- HFAG 2009
aeb.meas.add.single("m_tau", 1776.7673082, 0.1507259)

##
## compute phase space factors for Bmu/Be universality
##
aeb.meas.expr.add("xe", quote(m_e^2/m_tau^2))
aeb.meas.expr.add("xmu", quote(m_mu^2/m_tau^2))
rc = aeb.meas.expr.add("fx_e",
  quote(1 -8*xe + 8*xe^3 - xe^4 - 12*xe^2*log(xe)));
rc = aeb.meas.expr.add("fx_mu",
  quote(1 -8*xmu + 8*xmu^3 - xmu^4 - 12*xmu^2*log(xmu)));
aeb.meas.expr.add("fx_mu_by_e", quote(fx_mu/fx_e))

##
## Bmu/Be = f(m_mu^2/m_tau^2) / f(m_e^2/m_tau^2) = 0.972565 +- 0.000009 -- PDG 2004
## http://pi.physik.uni-bonn.de/~brock/teaching/vtp_ss06/doc/davier_0507078.pdf
## updated to HFAG 2009 by Swagato
##
## B_tau_mu_by_e_th.val = 0.972565
## B_tau_mu_by_e_th.err = 0.000009
## B_tau_mu_by_e_th.val = 0.972558
## B_tau_mu_by_e_th.err = 4.50333e-06
## aeb.meas.add.single("B_tau_mu_by_e_th", B_tau_mu_by_e_th.val, B_tau_mu_by_e_th.err)

##
## model matrix
## multiplied by the vector of theory parameters returns the vector of measurement types
##
## for universality Be = B(tau -> e nu nubar) = Be_univ
## -          1 * Be_univ = Be
## - fx_mu_by_e * Be_univ = Bmu
##
fx_mu_by_e.val = quant.val["fx_mu_by_e"]
names(fx_mu_by_e.val) = NULL
aeb.meas.fit.add("Gamma5univ", c(Gamma5=1, Gamma3=fx_mu_by_e.val))

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
## SU3 breaking correction
##
##--- POS(KAON)08, A.Pich, Theoretical progress on the Vus determination from tau decays
deltaR.su3break.val = 0.216
deltaR.su3break.err = 0.016
##--- E. Gamiz, M. Jamin, A. Pich, J. Prades, F. Schwab, V_us and m_s from hadronic tau decays
deltaR.su3break.val = 0.218
deltaR.su3break.err = 0.026

##--- PhysRevD.74.074009
aeb.meas.add.single("m_s", 94./1000., 6./1000.)
##--- E.Gamiz, M.Jamin, A.Pich, J.Prades, F.Schwab, |V_us| and m_s from hadronic tau decays
deltaR.su3break.val = 0.240
deltaR.su3break.err = 0.032
deltaR.su3break.val = 0.1544 + 9.3*0.094^2 + 0.0034 
deltaR.su3break.err = 0.032
aeb.meas.add.single("deltaR_su3break_pheno", 0.1544, 0.0037)
aeb.meas.add.single("deltaR_su3break_msd2", 9.3, 3.4)
aeb.meas.add.single("deltaR_su3break_remain", 0.0034, 0.0028)
aeb.meas.expr.add("deltaR_su3break", quote(deltaR_su3break_pheno + deltaR_su3break_msd2*m_s^2 + deltaR_su3break_remain))
## aeb.meas.add.single("deltaR_su3break", deltaR.su3break.val, deltaR.su3break.err)

##--- add R_tau as function of quantities
aeb.meas.expr.add("R_tau", quote(1/Gamma5univ -1 -fx_mu_by_e))
##--- add R_tau_s = B(tau -> Xs nu) / Be_univ
aeb.meas.expr.add("R_tau_s", quote(Gamma110/Gamma5univ))
##--- add R_tay_VA = R_tau - R_tau_s
aeb.meas.expr.add("R_tau_VA", quote(R_tau - R_tau_s))
##--- add Vus
aeb.meas.expr.add("Vus", quote(sqrt(R_tau_s/(R_tau_VA/Vud^2 - deltaR_su3break))))

##--- theory error
Vus.err.th = abs(quant.cov["Vus", "deltaR_su3break"])/quant.err["deltaR_su3break"]
Vus.err.exp = sqrt(quant.err["Vus"]^2 - Vus.err.th^2)

nsigma = (quant.val["Vus"] - 0.2255)/quadrature(c(quant.err["Vus"], 0.0010))

display.names = c(Gamma110.names,
  "Gamma110", "Gamma110_pdg09", "fx_mu_by_e", "Gamma5univ",
  "R_tau", "R_tau_s", "R_tau_VA", "deltaR_su3break", "Vus")
show(rbind(cbind(val=quant.val[display.names], err=quant.err[display.names]),
           Vus_err_exp = c(Vus.err.exp, 0),
           Vus_err_th = c(Vus.err.th, 0),
           Vus_err_perc = c(quant.err["Vus"]/quant.val["Vus"]*100, 0),
           Vus_th_err_perc = c(Vus.err.th /quant.val["Vus"]*100, 0),
           nsigma = c(val=nsigma, err=0)
           ))

