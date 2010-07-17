#!/usr/bin/env Rscript

## ////////////////////////////////////////////////////////////////////////////
##
## aluelab-results.r
##
## - elaborate of alucomb.r fit results
##
## ////////////////////////////////////////////////////////////////////////////

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

## S035B33  Gamma40       G(pi- Kbar0 pi0 nu(tau))/G(total)
## S035C1   Gamma48       G(pi- K(S)0 K(L)0 nu(tau))/G(total)

##
## Sudan's paper
## B(tau -> KSKL pi nu) = 0.00087 +- 0.000085 +- 0.000030
##+++ include if approved in time
##

##--- measurement names matching a specific PDG Gamma
meas.match = function(gamma) {
  return( regexpr(paste("[.]", gamma, "[.]", sep=""), names(meas.val)) != -1 )
}

##--- sum in quadrature
quadrature = function(x) {
  return(sqrt(sum(x^2)))
}

tau.to.s.pdg08 = aeb.linear.comb.glob(c(
  Gamma10=1,  Gamma16=1, 
  Gamma23=1,  Gamma28=1,  Gamma35=1, 
  Gamma40=1,  Gamma85=1,  Gamma89=1,
  Gamma128=1))

##--- Gamma89  = K pi+pi- pi0 nu

##
## tau -> s Swagato
## - Gamma110 = Gamma10 + Gamma16 + Gamma35 + Gamma23 + Gamma40 + Gamma85 + Gamma128 + Gamma129 + Gamma144 + estimates
## - replace estimates with PDG values
##
## Gamma110 = Gamma10+Gamma16+Gamma23+Gamma28+Gamma35+Gamma40+Gamma44+Gamma53+Gamma85+Gamma89+Gamma128+Gamma130+Gamma132+Gamma144
##
sel.swb = c(
  Gamma10=1,  Gamma16=1, 
  Gamma23=1,  Gamma28=1,  Gamma35=1, 
  Gamma40=1,  Gamma44=1,  Gamma53=1,
  Gamma85=1,  Gamma89=1,
  Gamma128=1,
  Gamma130=1, Gamma132=1,
  Gamma96=1.699387)
tau.to.s.swb = aeb.linear.comb.glob(sel.swb)
sel.swb.names = names(sel.swb[sel.swb != 0])

##
## Bmu/Be = f(m_mu^2/m_tau^2) / f(m_e^2/m_tau^2) = 0.972565 +- 0.000009 -- PDG 2004
## http://pi.physik.uni-bonn.de/~brock/teaching/vtp_ss06/doc/davier_0507078.pdf
##
B.tau.mu.by.e.th = 0.972565
B.tau.mu.by.e.th.err = 0.000009

sel = c(Gamma5=1, Gamma3=B.tau.mu.by.e.th)
sel.names = names(sel)

mm.val = quant.val[sel.names]
mm.cov = quant.cov[sel.names, sel.names]
delta = matrix(sel, 2, 1)

mm.invcov = solve(t(delta) %*% solve(mm.cov) %*% delta)
vv.val = mm.invcov %*% t(delta) %*% solve(mm.cov) %*% mm.val
vv.cov = mm.invcov
vv.err = sqrt(diag(vv.cov))

##--- Vud
Vud.val = 0.97425
Vud.err = 0.00022

Rtau = 1/vv.val -1 -B.tau.mu.by.e.th
Rtau.s = tau.to.s.swb["val"] / vv.val
Rtau.ns = Rtau - Rtau.s

##--- 0.240 +- 0.032
DeltaR.ope = 0.240
DeltaR.ope.err = 0.032

Vus.val = sqrt(Rtau.s/(Rtau.ns/Vud.val^2 - DeltaR.ope))
Vus.err = 0

show(rbind(cbind(val=quant.val[sel.swb.names], err=quant.err[sel.swb.names]),
           Gamma110 = tau.to.s.swb,
           B_tau_e_univ = data.frame(val=vv.val, err=vv.err),
           R_tau_h = data.frame(val=Rtau, err=0),
           R_tau_S = data.frame(val=Rtau.s, err=0),
           R_tau_VA = data.frame(val=Rtau.ns, err=0),
           Vus = data.frame(val=Vus.val, err=Vus.err)
           ))
