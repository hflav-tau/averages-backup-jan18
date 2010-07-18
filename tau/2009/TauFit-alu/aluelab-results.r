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

require(numDeriv)

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
quant.err = quant3.err
quant.corr = quant3.corr
quant.cov = quant3.cov

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
## minimum chisq fit constraining Bmu/Be from Standard Model,
## using Be, Bmu and their covariance
##

##
## Bmu/Be = f(m_mu^2/m_tau^2) / f(m_e^2/m_tau^2) = 0.972565 +- 0.000009 -- PDG 2004
## http://pi.physik.uni-bonn.de/~brock/teaching/vtp_ss06/doc/davier_0507078.pdf
##
B_tau_mu_by_e_th.val = 0.972565
B_tau_mu_by_e_th.err = 0.000009
aeb.meas.add.single("B_tau_mu_by_e_th", B_tau_mu_by_e_th.val, B_tau_mu_by_e_th.err)

##
## model matrix
## multiplied by the vector of theory parameters returns the vector of measurement types
##
## for universality Be = B(tau -> e nu nubar) = Be_univ
## -                    1 * Be_univ = Be
## - B_tau_mu_by_e_th.val * Be_univ = Bmu
##
aeb.meas.fit.add("Gamma5univ", c(Gamma5=1, Gamma3=B_tau_mu_by_e_th.val))
B_tau_e_univ.val = quant.val["Gamma5univ"]
B_tau_e_univ.err = quant.err["Gamma5univ"]

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
## POS(KAON)08, A.Pich, Theoretical progress on the Vus determination from tau decays
##
deltaR.su3break.val = 0.240
deltaR.su3break.err = 0.032
deltaR.su3break.val = 0.216
deltaR.su3break.err = 0.016
aeb.meas.add.single("deltaR_su3break", deltaR.su3break.val, deltaR.su3break.err)

##--- add R_tau as function of quantities
aeb.meas.expr.add("R_tau", quote(1/Gamma5univ -1 -B_tau_mu_by_e_th))
##--- add R_tau_s = B(tau -> Xs nu) / Be_univ
aeb.meas.expr.add("R_tau_s", quote(Gamma110/Gamma5univ))
##--- add R_tay_VA = R_tau - R_tau_s
aeb.meas.expr.add("R_tau_VA", quote(R_tau - R_tau_s))
##--- add Vus
aeb.meas.expr.add("Vus", quote(sqrt(R_tau_s/(R_tau_VA/Vud^2 - deltaR_su3break))))

display.names = c(Gamma110.names,
  "Gamma110", "Gamma110_pdg09", "B_tau_mu_by_e_th", "Gamma5univ",
  "R_tau", "R_tau_s", "R_tau_VA", "Vus")
show(rbind(cbind(val=quant.val[display.names], err=quant.err[display.names])
           ))
