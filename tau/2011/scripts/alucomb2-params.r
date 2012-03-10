#!/usr/bin/env Rscript

##
## alucomb2-params.r
##
## produce file parameters.input with constants used for
## - updating systematic terms
## - defining constraint equations
##

require(proto, quietly=TRUE)
source("../../../Common/bin/alucomb2-getpdg.r")
source("../../../Common/bin/aluelab2.r")

## ////////////////////////////////////////////////////////////////////////////
## functions

alucomb2.fmt.param = function(label, val, err.p, err.m=NULL) {
  val.str = sprintf("%-12.5g", val)
  if (!is.null(err.m) && err.p != -err.m) {
    exc.str = paste("+", sprintf("%-11.5g", err.p), " ", sprintf("%.5g", err.m), sep="")
  } else {
    exc.str = paste("+-", sprintf("%.5g", err.p), sep="")
  }
  paste(format(label, width=16), val.str, exc.str)
}

alucomb2.print.param = function(label, val, err.p, err.m=NULL) {
  cat("  ", alucomb2.fmt.param(label, val, err.p, err.m), "\n", sep="")
}

##
## class for updatinf params
##
Alucomb2Param = proto()

Alucomb2Param$new = function(., val, cov) {
  proto(.)
}

Alucomb2Param$get.pdg.node = function(., node) {
  pdginfo = alucomb2.getpdg(node)
  pdg.src = "PDG11 fit"
  rc = pdginfo["PDG FIT",]
  if (is.null(rc[1])) {
    rc = pdginfo[["PDG AVERAGE"]]
    pdg.src = "PDG11 average"
  }
  if (is.null(rc)) return(c(NA, NA))
  quant = as.numeric(rc[1:2])
  names(quant) = colnames(rc[1:2])
  attr(quant, "pdg.source") = pdg.src
  return(quant)
}

Alucomb2Param$get.pdg.tau.gamma = function(., gamma) {
  .$get.pdg.node(.$node.from.tau.gamma[gamma])
}

Alucomb2Param$print.pdg.node = function(., comment, label, node, factor=1) {
  rc = .$get.pdg.node(node)
  pdg.src = attr(rc, "pdg.source")
  rc = rc * factor
  cat("  # --- ", comment, ", ", pdg.src, "\n", sep="")
  alucomb2.print.param(label, rc[1], rc[2])
}

Alucomb2Param$print.pdg.tau.gamma = function(., comment, label, gamma, factor = 1) {
  .$print.pdg.node(comment, label, .$node.from.tau.gamma[gamma], factor)
}

## ////////////////////////////////////////////////////////////////////////////
## code

alucomb2.params = function() {
  ##--- get tau gamma to PDG node conversion
  quant.cards = readLines("../Common-Mar2011/quantities.input")
  quant.cards = quant.cards[grep("^QUANTITY", quant.cards)]
  rc = str_match(quant.cards, ignore.case("^QUANTITY\\s+(\\S+)\\s+node\\s+(\\S+)"))
  Alucomb2Param$node.from.tau.gamma = character(0)
  Alucomb2Param$node.from.tau.gamma[rc[,2]] = rc[,3]
  rm(rc)

  cat("#\n")
  cat("# parameters used for possibly updating systematic terms\n")
  cat("# use PDG 2011 when possible\n")
  cat("#\n")
  cat("PARAMETERS\n")
  
  cat("  # --- e+e- -> tau+tau- cross section from Swagato et al. paper\n")
  alucomb2.print.param("sigmataupmy4s", 0.919, 0.003)
  
  Alucomb2Param$print.pdg.tau.gamma("tau -> pi- pi0 nu, Gamma14", "PimPizNu", "Gamma14")
  Alucomb2Param$print.pdg.tau.gamma("tau -> K- K0S, 1/2 * Gamma37", "KmKzsNu", "Gamma37", 1/2)
  
  ## Alucomb2Param$print.pdg.tau.gamma("tau -> pi- pi0 K0S, 1/2 * Gamma40", "PimPizKzsNu", "Gamma40", 1/2)
  ##--- use HFAG-tau mid-2011 fit, which uses a BaBar preliminary result
  htag.tau.mar11.PimPizKzNu = c(val=0.00344957, err=0.000147272)
  cat("  # --- tau -> pi- pi0 K0S, 1/2 * Gamma40, HFAG-tau mid-2011 fit using BaBar prelim. meas.\n")
  alucomb2.print.param("PimPizKzsNu", 1/2*htag.tau.mar11.PimPizKzNu[1], 1/2*htag.tau.mar11.PimPizKzNu[2])
  
  Alucomb2Param$print.pdg.tau.gamma("tau -> K- pi0 K0S, 1/2 * Gamma42", "KmPizKzsNu", "Gamma42", 1/2)
  Alucomb2Param$print.pdg.tau.gamma("tau -> pi- K0S K0S, Gamma47", "PimKzsKzsNu", "Gamma47")
  Alucomb2Param$print.pdg.tau.gamma("tau -> pi-pi-pi+ (ex.K0), Gamma60", "PimPimPipNu", "Gamma60")
  Alucomb2Param$print.pdg.tau.gamma("tau -> pi-pi-pi+pi0 (ex.K0), Gamma69", "PimPimPipPizNu", "Gamma69")
  Alucomb2Param$print.pdg.tau.gamma("tau -> K-pi-pi+, Gamma85", "PimKmPipNu", "Gamma85")
  Alucomb2Param$print.pdg.tau.gamma("tau -> K-K+pi, Gamma93", "PimKmKpNu", "Gamma93")

  ##
  ## see http://hfag.googlecode.com/svn/trunk/tau/2009/Docs/kaons-isospin-notes.pdf
  ## Gamma46 = G(pi- K0 Kbar0 nu(tau)) / G(total)
  ## Gamma47 = G(pi- K(S)0 K(S)0 nu(tau)) / G(total)
  ## Gamma48 = G(pi- K(S)0 K(L)0 nu(tau)) / G(total)
  ## Gamma804 = G(pi- K(L)0 K(L)0 nu(tau)) / G(total)
  ## Gamma93 = G(pi- K- K+ nu(tau)) / G(total)
  ## one has:
  ## - GammaK0K0 = GammaKPKM
  ## - GammaKSKS = GammaKLKL
  ## - GammaK0K0 = GammaKSKL + GammaKSKS + GammaKLKL
  ## but one does not have:
  ## - GammaKSKS = GammaKSKL
  ## to improve GammaKSKL one can use the above equalities
  ## GammaKSKS_fit = average(KSKS, KLKL) (however, there no PDG value for KLKL)
  ## GammaKSKL_fit = average(KSKL, KPKM - 2*GammaKSKS_fit)
  ## 

  ##--- most straightforward definition
  ## Alucomb2Param$print.pdg.tau.gamma("tau -> pi- K0S K0S, Gamma47", "PimKzsKzlNu", "Gamma47")

  tau.piksks = Alucomb2Param$get.pdg.tau.gamma("Gamma47")
  tau.pikskl = Alucomb2Param$get.pdg.tau.gamma("Gamma48")
  tau.pikpkm = Alucomb2Param$get.pdg.tau.gamma("Gamma93")
  
  quant = StatComb$new()
  
  quant$meas.add.single("BSS", tau.piksks[1], tau.piksks[2])
  quant$meas.add.single("BSL", tau.pikskl[1], tau.pikskl[2])
  quant$meas.add.single("BPM", tau.pikpkm[1], tau.pikpkm[2])
  
  quant$meas.expr.add("BSL_other", quote(BPM - 2*BSS))
  quant$meas.fit.add("BSL_fit", c(BSL=1, BSL_other=1))

  rc = quant$val.err("BSL_fit")
  cat("  #\n")
  cat("  # tau -> pi KS KL nu\n")
  cat("  # best fit of PDG11 fit results (neglecting correlations)\n")
  cat("  # weighted average of BRs to KSKL and K+K- - 2KSKS\n")
  cat("  # KSKL = Gamma48, K+K- = Gamma93, KSKS = Gamma47\n")
  cat("  #\n")
  alucomb2.print.param("PimKzsKzlNu", rc[1], rc[2])
  rm(quant)

  cat("\n")
  cat("#\n")
  cat("# the following parameters are used in the constraint equations\n")
  cat("# (the errors are actually not used)\n")
  cat("#\n")
  cat("PARAMETERS\n")

  BR_eta_neutral = Alucomb2Param$get.pdg.node("S014R21")
  BR_eta_charged = c(1 - BR_eta_neutral[1], BR_eta_neutral[2])

  Alucomb2Param$print.pdg.node("eta -> gamma gamma, node S014R34", "BR_eta_2gam", "S014R34")

  ##Alucomb2Param$print.pdg.node("eta -> neutral, node S014R21", "BR_eta_neutral", "S014R21")
  cat("  # eta -> neutral, node S014R21, ", attr(BR_eta_neutral, "pdg.source"), "\n", sep="")
  alucomb2.print.param("BR_eta_neutral", BR_eta_neutral[1], BR_eta_neutral[2])
  
  Alucomb2Param$print.pdg.node("eta -> 3pi0, node S014R52", "BR_eta_3piz", "S014R52")
  Alucomb2Param$print.pdg.node("eta -> pi+pi-pi0, node S014R53", "BR_eta_pimpippiz", "S014R53")

  cat("  # eta -> charged, node S014R21, ", attr(BR_eta_charged, "pdg.source"), "\n", sep="")
  alucomb2.print.param("BR_eta_charged", BR_eta_charged[1], BR_eta_charged[2])
  
  Alucomb2Param$print.pdg.node("K0S -> 2piz, node S012R2", "BR_KS_2piz", "S012R2")
  Alucomb2Param$print.pdg.node("K0S -> pi+pi-, node S012R1", "BR_KS_pimpip", "S012R1")

  ## Alucomb2Param$print.pdg.node("omega -> pi+pi-pi0, node M001R21", "BR_om_pimpippiz", "M001R21")
  cat("  # --- omega -> pi+pi-pi0, from PDG11 fit\n")
  alucomb2.print.param("BR_om_pimpippiz", 89.2/100, 0.7/100)
  
  ## Alucomb2Param$print.pdg.node("omega -> pi+pi-, node M001R15", "BR_om_pimpip", "M001R15")
  cat("  # --- omega -> pi+pi-, from PDG11 fit\n")
  alucomb2.print.param("BR_om_pimpip", 1.53/100, +0.11/100, -0.13/100)  

  ## Alucomb2Param$print.pdg.node("omega -> pi0gamma, node M001R28", "BR_om_pizgamma", "M001R28")
  cat("  # --- omega -> pi0gamma, from PDG11 fit\n")
  alucomb2.print.param("BR_om_pizgamma", 8.28/100, 0.28/100)  

  Alucomb2Param$print.pdg.node("phi -> K+K-, node M004R1", "BR_phi_KmKp", "M004R1")
  Alucomb2Param$print.pdg.node("phi -> K0SK0L, node M004R2", "BR_phi_KSKL", "M004R2")

  cat("  #\n")
  cat("  BRA_Kz_KS_KET		0.5	    +-0		# |<Kz | KS>|^2\n")
  cat("  BRA_Kz_KL_KET		0.5         +-0         # |<Kz | KS>|^2\n")
  cat("  BRA_Kzbar_KS_KET	0.5	    +-0		# |<Kzbar | KS>|^2\n")
  cat("  BRA_Kzbar_KL_KET	0.5	    +-0		# |<Kzbar | KS>|^2\n")
  cat("  # --- |<Kz Kzbar | KS KS>|^2\n")
  cat("  BRA_KzKzbar_KLKL_KET_by_BRA_KzKzbar_KSKS_KET  1 +-0\n")
}

args = commandArgs(TRUE)
alucomb2.params() 
