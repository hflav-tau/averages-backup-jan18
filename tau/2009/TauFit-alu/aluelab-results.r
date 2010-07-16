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

##
## Gamma44 = pi K0 pi0 pi0 nu
## PDG 2009 ( 2.6 ± 2.4 ) × 10~3 
##
aeb.meas.add.single("Gamma44", 2.6e-4, 2.4e-4)

##
## Gamma53 =  K0 h+ h~ h~ nu
## PDG 2009 ( ( 2.3 ± 2.0 ) × 10~4
##
aeb.meas.add.single("Gamma53", 2.3e-4, 2.0e-4)

##
## Gamma129 = K*(892) eta nu
## PDG 2009 ( 1.38 ± 0.15 ) × 10~4 
##
aeb.meas.add.single("Gamma129", 1.38e-4, 0.15e-4)

##
## Gamma130 = eta K pi0 nu
## PDG 2009 ( 4.8 ± 1.2 ) × 10 ~ 5
##
aeb.meas.add.single("Gamma130", 4.8e-5, 1.2e-5)

##
## Gamma132 = eta K0 pi nu
## PDG 2009 ( 9.3 ± 1.5 ) × 10 ~ 5
##
aeb.meas.add.single("Gamma132", 9.3e-5, 1.5e-5)

##
## Gamma144 = K phi nu
## KmPhiNu HFAG 2009 3.704330e-05 +- 3.259888e-06 (S=1.3)
##
aeb.meas.add.single("Gamma144", 3.704330e-05, 3.259888e-06)

##
## Gamma151 = K omega nu
## PDG 2009 ( 4.1 ± 0.9 ) × 10 ~ 4
##
aeb.meas.add.single("Gamma151", 4.1e-4, 0.9e-4)

##
## B[tau- -> (K3pi)- nu (ex. K0, omega,eta)]  : 1e-2 : 0.074 ± 0.030    
##
aeb.meas.add.single("Km3PiNu", 0.074e-2, 0.030e-2)

##
## B[tau- -> (K4pi)- nu]                      : 1e-2 : 0.011 ± 0.007
##
aeb.meas.add.single("Km4PiNu", 0.011e-2, 0.007e-2)

## S035B33  Gamma40       G(pi- Kbar0 pi0 nu(tau))/G(total)
## S035C1   Gamma48       G(pi- K(S)0 K(L)0 nu(tau))/G(total)

##
## Sudan's paper
## B(tau -> KSKL pi nu) = 0.00087 +- 0.000085 +- 0.000030
## - va aggiunta alle misure se approvata in tempo
##

##--- measurement names matching a specific PDG Gamma
meas.match = function(gamma) {
  return( regexpr(paste("[.]", "Gamma40", "[.]", sep=""), names(meas)) != -1 )
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
sel.swb = c(
  Gamma10=1,  Gamma16=1, 
  Gamma23=1,  Gamma28=1,  Gamma35=1, 
  Gamma40=1,  Gamma44=1,  Gamma53=1,
  Gamma85=1,  Gamma89=1,
  Gamma128=1, Gamma129=0,
  Gamma130=1, Gamma132=1,
  Gamma144=1, Gamma151=1,
  Km3PiNu=0, Km4PiNu=0)
tau.to.s.swb = aeb.linear.comb.glob(sel.swb)
sel.swb.names = names(sel.swb[sel.swb != 0])

##
## Bmu/Be = f(m_mu^2/m_tau^2) / f(m_e^2/m_tau^2) = 0.972565 +- 0.000009 -- PDG 2004
## http://pi.physik.uni-bonn.de/~brock/teaching/vtp_ss06/doc/davier_0507078.pdf
##
sel = c(Gamma5=1, Gamma3=1/0.972565)
sel.names = names(sel)

mm.val = quant.val[sel.names]
mm.cov = quant.cov[sel.names, sel.names]
delta = matrix(c(1, 0.972565), 2, 1)

mm.invcov = solve(t(delta) %*% solve(mm.cov) %*% delta)
vv.val = mm.invcov %*% t(delta) %*% solve(mm.cov) %*% mm.val
vv.cov = mm.invcov
vv.err = sqrt(diag(vv.cov))

##--- Vud
Vud.val = 0.97425
Vud.err = 0.00022

Rtau = (1 - quant.val["Gamma3"])/quant.val["Gamma5"] -1
Rtau.s = tau.to.s.swb["val"] / quant.val["Gamma5"]
Rtau.ns = Rtau-Rtau.s

##--- 0.240 +- 0.032
DeltaR.ope = 0.240
DeltaR.ope.err = 0.032

Vus.val = sqrt(Rtau.s/(Rtau.ns/Vud.val^2 - DeltaR.ope))

show(rbind(cbind(val=quant.val[sel.swb.names], err=quant.err[sel.swb.names]),
           Gamma110=tau.to.s.swb,
           B.tau.e.univ=data.frame(val=vv.val, err=vv.err),
           B.tau.s=data.frame(val=Vus.val, err=0)))
