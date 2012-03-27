#!/usr/bin/env Rscript

require(stringr, quietly=TRUE)
source("../../../Common/bin/aluelab2.r")

## ////////////////////////////////////////////////////////////////////////////
## functions

##
## get list of quantity names in the definition of a quantity
## using the relevant constraint equation
## used since Oct 2010
## since March 2012, use aluelab.get.str.expr, needs alucomb2 format
##
aluelab.get.quant.names = function(quant.name, combination) {
  var.comb = combination$constr.all.comb[[paste(quant.name, ".c", sep="")]]
  if (is.null(var.comb) || is.na(var.comb)) {
    var.comb = combination$constr.all.comb[[paste(quant.name, ".coq", sep="")]]
  }
  if (is.null(var.comb)) {
    var.comb = combination$constr.lin.comb[[paste(quant.name, ".coq", sep="")]]
  }
  if (is.null(var.comb)) {
    var.comb = combination$constr.all.comb[[quant.name]]
  }
  if (is.null(var.comb)) {
    var.comb = combination$constr.lin.comb[[quant.name]]
  }
  var.comb = var.comb[names(var.comb) != quant.name]
  var.names = names(var.comb[var.comb != 0])
  return(var.names)
}

##
## get string expression for a quantity, from its constraint equations
## needs alucomb2 format
##
aluelab.get.str.expr = function(quant.name, combination) {
  str.expr = combination$constr.all.str.expr[[paste(quant.name, ".c", sep="")]]
  str.expr = str_match(str.expr, "-([[:alnum:]]+)\\s+[+]\\s+(.*\\S)\\s*$")[,3]
  ##--- remove outer braces
  str.expr = gsub("(\\(((?:[^()]++|(?1))+)*+\\))", "\\2", str.expr, perl=TRUE)
  return(str.expr)
}

## ////////////////////////////////////////////////////////////////////////////
## code

aluelab.results = function(args) {
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
    quant.cov = quant.sf.cov
  }

  quant = StatComb$new(quant.val, quant.cov)
  quant.names = names(quant.val)

  comb.params = lapply(combination$params, function(x) unname(x["value"]))
  quant$params.add(comb.params)
  
  ##
  ## recover some BRs as function of others
  ##
  if (!any("Gamma89" == quant.names)) {
    rc = quant$meas.expr.add("Gamma89", quote(Gamma803 + BR_om_pimpippiz*Gamma151))
  }
    
  ##
  ## PDG 2009 definition of Gamma110 = B(tau -> Xs nu)
  ##
  Gamma110_pdg09.comb = c(
    Gamma10=1,  Gamma16=1, 
    Gamma23=1,  Gamma28=1,  Gamma35=1, 
    Gamma40=1,  Gamma85=1,  Gamma89=1,
    Gamma128=1)
  quant$meas.comb.add("Gamma110_pdg09", Gamma110_pdg09.comb)
  
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

  ##
  ## end Oct2010 update
  ## define Gamma110 using the alucomb.r Gamma110 constraint
  ## this is more robust against changes in the fitting 
  ##
  Gamma110.names = aluelab.get.quant.names("Gamma110", combination)
  
  ##
  ## March 2012, get Gamma110 from its constraint equation
  ##
  Gamma110.str.expr = aluelab.get.str.expr("Gamma110", combination)
  Gamma110.comb = quant$str.to.comb(Gamma110.str.expr)
  Gamma110.names = names(Gamma110.comb)
  
  Gamma110.str.expr = paste(Gamma110.comb, Gamma110.names, sep="*", collapse=" + ")
  cat("Gamma110 = ", Gamma110.str.expr, "\n", sep="")
  
  if (!any("Gamma151" == quant.names)) {
    ##--- remove Gamma151 (tau -> K omega nu) if not present in fit for backward compatibility
    Gamma110.comb = Gamma110.comb[Gamma110.names != "Gamma151"]
    cat("warning, Gamma151 removed from Gamma110 definition\n")
  }
  
  ##--- list of all tau BRs that are not leptonic and not strange, i.e. not-strange-hadronic
  B.tau.VA.names = setdiff(aluelab.get.quant.names("GammaAll", combination), Gamma110.names)
  B.tau.VA.names = setdiff(B.tau.VA.names, c("Gamma5", "Gamma3", "Gamma998"))
  quant$meas.expr.add("B_tau_VA", parse(text=paste(B.tau.VA.names, collapse="+")))
  quant$meas.expr.add("B_tau_VA_unitarity", quote(1-Gamma5-Gamma3-Gamma110))
  quant$meas.expr.add("B_tau_s_unitarity", quote(1-Gamma5-Gamma3-B_tau_VA))
  
  ##
  ## add measurements to compute universality improved Be
  ##
  
  ##--- from PDG 2009, 2011
  quant$meas.add.single("m_e", 0.510998910, 0.000000013)
  quant$meas.add.single("m_mu", 105.658367, 0.000004)
  quant$meas.add.single("tau_tau", 290.6e-15, 1.0e-15)
  ##--- m_tau HFAG 2009
  ## quant$meas.add.single("m_tau", 1776.7673082, 0.1507259)
  ##--- m_tau PDG 2011
  quant$meas.add.single("m_tau", 1776.82, 0.16)
 
  quant$meas.add.single("m_pi", 139.57018, 0.00035)
  quant$meas.add.single("tau_pi", 2.6033e-8, 0.0005e-8)
  quant$meas.add.single("m_K", 493.677, 0.016)
  quant$meas.add.single("tau_K", 1.2380e-8, 0.0021e-8)

  ##--- from PDG 2010, 2011
  quant$meas.add.single("m_W", 80.399*1e3, 0.023*1e3)
  quant$meas.add.single("tau_mu", 2.197034e-6, 0.000021e-6)
  
  ##
  ## Be, from unitarity = 1 - Bmu - B_VA - B_s
  ## Bmu, from unitarity = 1 - Be - B_VA - B_s
  ##
  quant$meas.expr.add("Be_unitarity", quote(1 - Gamma3 - B_tau_VA - Gamma110))
  quant$meas.expr.add("Bmu_unitarity", quote(1 - Gamma5 - B_tau_VA - Gamma110))
  
  ##--- understand if there was no unitarity constraint
  no.unit.constr.flag =
    ("Gamma998" %in% aluelab.get.quant.names("GammaAll", combination) ||
     "Gamma998" %in% aluelab.get.quant.names("Unitarity", combination))

  if (no.unit.constr.flag) {
    ##
    ## unconstrained fit
    ## fit best values for Be, Bmu, B_tau_VA, B_tau_s
    ## using both the direct measurements and the result from the unitarity constraint
    ##
    if (flag.unitarity) {
      ##--- use unitarity constraint to compute Be, Bmu, B_tau_VA/s
      quant$meas.fit.add("Be_fit", c(Gamma5=1, Be_unitarity=1))
      quant$meas.fit.add("Bmu_fit", c(Gamma3=1, Bmu_unitarity=1))
      quant$meas.fit.add("B_tau_VA_fit", c(B_tau_VA=1, B_tau_VA_unitarity=1))
      quant$meas.fit.add("B_tau_s_fit", c(Gamma110=1, B_tau_s_unitarity=1))
    } else {
      ##--- direct measurements for leptonic BRs and strange hadronic BR
      quant$meas.expr.add("Be_fit", quote(Gamma5))
      quant$meas.expr.add("Bmu_fit", quote(Gamma3))
      quant$meas.expr.add("B_tau_s_fit", quote(Gamma110))
      if (flag.vadirect) {
        ##--- direct measurement for non-strange hadronic BR
        quant$meas.expr.add("B_tau_VA_fit", quote(B_tau_VA))
      } else {
        ##
        ## non-strange hadronic BR by subtraction, however note that
        ## here we do not use the improved Be from universality
        ## (a Be fit using Be, Bmu and tau lifetime) as in
        ## M.Davier et al, RevModPhys.78.1043, arXiv:hep-ph/0507078
        ##
        quant$meas.expr.add("B_tau_VA_fit", quote(1-Gamma5-Gamma3-Gamma110))
      }
    }
  } else {
    ##
    ## unitarity constrained fit
    ## the fitted values have the unitarity constraint already
    ##
    rc = quant$meas.expr.add("Be_fit", quote(Gamma5))
    rc = quant$meas.expr.add("Bmu_fit", quote(Gamma3))
    rc = quant$meas.expr.add("B_tau_VA_fit", quote(B_tau_VA))
    rc = quant$meas.expr.add("B_tau_s_fit", quote(Gamma110))
  }
  
  ##
  ## compute phase space factors for Bmu/Be universality
  ##
  ##--- phase space factor, function of lepton masses
  phspf = quote(1 -8*x + 8*x^3 - x^4 - 12*x^2*log(x))
  ##--- phase space factors for e/tau, mu/tau, e/mu
  rc = quant$meas.expr.add("phspf_mebymtau",  esub.expr(phspf, list(x=quote(m_e^2/m_tau^2))))
  rc = quant$meas.expr.add("phspf_mmubymtau", esub.expr(phspf, list(x=quote(m_mu^2/m_tau^2))))
  rc = quant$meas.expr.add("phspf_mebymmu", esub.expr(phspf, list(x=quote(m_e^2/m_mu^2))))
  rc = quant$meas.expr.add("Bmu_by_Be_th", quote(phspf_mmubymtau/phspf_mebymtau))
  
  ##--- Be from Bmu
  quant$meas.expr.add("Be_from_Bmu", quote(phspf_mebymtau/phspf_mmubymtau *Bmu_fit))
  
  ##
  ## rad. corrections from to get Be from tau lifetime
  ## values from 10.1103/RevModPhys.78.1043 p.1047, arXiv:hep-ph/0507078v2 p.7, could be recomputed
  ## - delta^L_gamma = 1 + alpha(mL)/2pi * (25/4 - pi^2)
  ## - delta^L_W = 1 + 3/5* m_L^2/M_W^2
  ##
  quant$params.add(c(delta_mu_gamma=(1-42.4e-4), delta_tau_gamma=(1-43.2e-4)))
  quant$meas.expr.add("delta_mu_W", quote(1 + 3/5*m_mu^2/m_W^2))
  quant$meas.expr.add("delta_tau_W", quote(1 + 3/5*m_tau^2/m_W^2))
  
  ##
  ## Be from tau lifetime
  ## Be= tau_tau / tau_mu (m_tau/m_mu)^5 f(m^2_e/m^2_tau)/f(m^2_e/m^2_mu) (delta^tau_gamma delta^tau_W)/(delta^mu_gamma delta^mu_W)
  ##
  quant$meas.expr.add("Be_from_taulife",
                      quote(tau_tau/tau_mu * (m_tau/m_mu)^5 * phspf_mebymtau/phspf_mebymmu
                            * (delta_tau_gamma*delta_tau_W) / (delta_mu_gamma*delta_mu_W)))
  ##
  ## Bmu from tau lifetime
  ## Bmu= tau_tau/tau_mu (m_tau/m_mu)^5 f(m^2_mu/m^2_tau)/f(m^2_e/m^2_mu) (delta^tau_gamma delta^tau_W)/(delta^mu_gamma delta^mu_W)
  ##
  quant$meas.expr.add("Bmu_from_taulife",
                      quote(tau_tau/tau_mu * (m_tau/m_mu)^5 * phspf_mmubymtau/phspf_mebymmu
                             * (delta_tau_gamma*delta_tau_W) / (delta_mu_gamma*delta_mu_W)))
  
  ##
  ## universality improved Be = B(tau -> e nu nubar (gamma))
  ## see arXiv:hep-ph/0507078v2 p.7, doi:10.1103/RevModPhys.78.1043 p.1047
  ##
  ## minimum chisq fit for Be_univ using, Be, Be from Bmu, Be from tau lifetime
  ##
  ## Bmu/Be = f(m_mu^2/m_tau^2) / f(m_e^2/m_tau^2)
  ## Be= tau_tau / tau_mu (m_tau/m_mu)^5 f(m^2_e/m^2_tau)/f(m^2_e/m^2_mu) (delta^tau_gamma delta^tau_W)/(delta^mu_gamma delta^mu_W)
  ##
  quant$meas.fit.add("Be_univ", c(Gamma5=1, Be_from_Bmu=1, Be_from_taulife=1))
  
  ##
  ## Vud
  ##
  ## arXiv:0710.3181v1 [nucl-th], 10.1103/PhysRevC.77.025501
  ## I.S.Towner, J.C.Hardy, An improved calculation of the isospin-symmetry-breaking corrections to superallowed Fermi beta decay
  ## also PDG 2010 review
  ##
  Vud.val = 0.97425
  Vud.err = 0.00022
  quant$meas.add.single("Vud", Vud.val, Vud.err)
  
  ##
  ## SU3 breaking correction, straight from papers
  ##
  
  ##--- POS(KAON)08, A.Pich, Theoretical progress on the Vus determination from tau decays
  deltaR.su3break.val = 0.216
  deltaR.su3break.err = 0.016
  ##--- E. Gamiz et al., Nucl.Phys.Proc.Suppl.169:85-89,2007, arXiv:hep-ph/0612154v1
  deltaR.su3break.val = 0.240
  deltaR.su3break.err = 0.032
  ## quant$meas.add.single("deltaR_su3break", deltaR.su3break.val, deltaR.su3break.err)
  
  ##
  ## SU3 breaking correction, recompute from data following
  ## E. Gamiz et al., Nucl.Phys.Proc.Suppl.169:85-89,2007, arXiv:hep-ph/0612154v1
  ##
  
  ##--- s quark mass, PhysRevD.74.074009
  quant$meas.add.single("m_s", 94, 6)
  ##--- PDG 2011
  ## quant$meas.add.single("m_s", 100, sqrt((20.^2 + 30.^2)/2.))
  
  ##--- E.Gamiz, M.Jamin, A.Pich, J.Prades, F.Schwab, |V_us| and m_s from hadronic tau decays
  quant$meas.add.single("deltaR_su3break_pheno", 0.1544, 0.0037)
  quant$meas.add.single("deltaR_su3break_msd2", 9.3, 3.4)
  quant$meas.add.single("deltaR_su3break_remain", 0.0034, 0.0028)
  quant$meas.expr.add("deltaR_su3break", quote(deltaR_su3break_pheno + deltaR_su3break_msd2*(m_s/1000)^2 + deltaR_su3break_remain))
  
  if (no.unit.constr.flag) {
    ##
    ## if using constrained fit with dummy mode (Gamma998), i.e. unconstrained fit
    ##
    if (flag.lepuniv) {
      ##--- determine R_tau using leptonic BRs and universality
      quant$meas.expr.add("R_tau", quote(1/Be_univ -1 -phspf_mmubymtau/phspf_mebymtau))
      quant$meas.expr.add("R_tau_s", quote(B_tau_s_fit/Be_univ))
      quant$meas.expr.add("R_tau_VA", quote(R_tau - R_tau_s))
    } else {
      ##
      ## use "fit" values computed in this script (possibly incorporating unitarity constraint)
      ##
      ##--- add R_tau_VA = R - R_tau_s
      quant$meas.expr.add("R_tau_VA", quote(B_tau_VA_fit/Be_univ))
      ##--- add R_tau_s = B(tau -> Xs nu) / Be_univ
      quant$meas.expr.add("R_tau_s", quote(B_tau_s_fit/Be_univ))
      ##--- add R_tau as function of quantities
      quant$meas.expr.add("R_tau", quote(R_tau_VA+R_tau_s))
    }
  } else {
    ##
    ## if using constrained fit without dummy mode (Gamma998), i.e. constrained fit
    ##
    ##--- add R_tau as function of quantities
    if (flag.lepuniv) {
      ##--- determine R_tau using leptonic BRs and universality
      quant$meas.expr.add("R_tau", quote(1/Be_univ -1 -phspf_mmubymtau/phspf_mebymtau))
    } else {
      ##--- use values computed in the alucomb.r fit, incorporating unitarity constraints
      quant$meas.expr.add("R_tau", quote((B_tau_VA+Gamma110)/Be_univ))
    }
    ##--- add R_tau_s = B(tau -> Xs nu) / Be_univ
    quant$meas.expr.add("R_tau_s", quote(Gamma110/Be_univ))
    ##--- add R_tau_VA = R_tau - R_tau_s
    quant$meas.expr.add("R_tau_VA", quote(R_tau - R_tau_s))
  }
  
  ##--- add Vus
  quant$meas.expr.add("Vus", quote(sqrt(R_tau_s/(R_tau_VA/Vud^2 - deltaR_su3break))))
  
  ##--- theory error
  Vus.err.th = abs(quant$cov()["Vus", "deltaR_su3break"])/quant$err()["deltaR_su3break"]
  Vus.err.exp = sqrt(quant$err()["Vus"]^2 - Vus.err.th^2)

  ##--- using Hardy-Towner 2009, PDG 2010
  Vus.unitarity.val = 0.2255
  Vus.unitarity.err = 0.0010
  quant$meas.add.single("Vus_uni", Vus.unitarity.val, Vus.unitarity.err)

  quant$meas.expr.add("Vus_mism", quote(Vus - Vus_uni))
  nsigma = quant$val()["Vus_mism"] / quant$err()["Vus_mism"]

  display.names = c(Gamma110.names,
    "Gamma5", "Be_unitarity", "Be_fit",
    "Gamma3", "Bmu_unitarity", "Bmu_fit",
    "Bmu_by_Be_th", "Be_from_Bmu", "Be_from_taulife", "Be_univ", "Bmu_from_taulife",
    "B_tau_VA", "B_tau_VA_unitarity", "B_tau_VA_fit",
    "Gamma110", "B_tau_s_unitarity", "B_tau_s_fit", "Gamma110_pdg09",
    "R_tau", "R_tau_s", "R_tau_VA", "deltaR_su3break", "Vus")

  print(rbind(cbind(val=quant$val()[display.names], err=quant$err()[display.names]),
             Vus_err_perc = c(quant$err()["Vus"]/quant$val()["Vus"]*100, 0),
             Vus_err_exp = c(Vus.err.exp, 0),
             Vus_err_th = c(Vus.err.th, 0),
             Vus_err_th_perc = c(Vus.err.th /quant$val()["Vus"]*100, 0),
             nsigma = c(val=nsigma, err=0)
             ))

  ##
  ## Lattice averages from http://arxiv.org/abs/0910.2928 and
  ## http://krone.physik.unizh.ch/~lunghi/webpage/LatAves/page7/page7.html
  ##
  quant$meas.add.single("f_K_by_f_pi", 1.1936, 0.0053)
  quant$meas.add.single("f_K", 156.1, 1.1)

  ##
  ## gtau/gmu using tau -> hnu / h -> mu nu
  ##
  quant$meas.add.single("pitoENu", 1.230e-4, 0.004e-4)
  quant$meas.add.single("pitoMuNu", 99.98770e-2, 0.00004e-2)
  quant$meas.add.single("KtoENu", 1.584e-5, 0.020e-5)
  quant$meas.add.single("KtoMuNu", 63.55e-2, 0.11e-2)
  ##--- from Marciano:1993sh,Decker:1994ea,Decker:1994dd
  quant$meas.add.single("delta_pi", 0.16e-2, 0.14e-2)
  quant$meas.add.single("delta_K", 0.90e-2, 0.22e-2)

  ##--- gtau/gmu using tau -> pi nu / pi -> mu nu
  quant$meas.expr.add("gtaubygmu_pi",
                      quote(sqrt(Gamma9/pitoMuNu *(2*m_pi*m_mu^2*tau_pi) /((1+delta_pi)*m_tau^3*tau_tau)*
                                 ((1-m_mu^2/m_pi^2)/(1-m_pi^2/m_tau^2))^2)))

  ##--- gtau/gmu using tau -> K nu / K -> mu nu
  quant$meas.expr.add("gtaubygmu_K",
                      quote(sqrt(Gamma10/KtoMuNu *(2*m_K*m_mu^2*tau_K) /((1+delta_K)*m_tau^3*tau_tau)*
                                 ((1-m_mu^2/m_K^2)/(1-m_K^2/m_tau^2))^2)))

  ##--- gtau/gmu using tau lifetime
  quant$meas.expr.add("gtaubygmu_tau", quote(sqrt(Gamma5/Be_from_taulife)))

  ##--- gtau/gmu average
  quant$meas.fit.add("gtaubygmu_fit", c(gtaubygmu_tau=1, gtaubygmu_pi=1, gtaubygmu_K=1))
  
  ##--- gtau/ge using tau lifetime
  quant$meas.expr.add("gtaubyge_tau", quote(sqrt(Gamma3/Bmu_from_taulife)))
  
  ##--- gmu / ge from tau -> mu / tau -> e
  quant$meas.expr.add("gmubyge_tau", quote(sqrt(Gamma3/Gamma5 * phspf_mebymtau/phspf_mmubymtau)))

  display.names = c(
    "gtaubygmu_tau", "gtaubygmu_pi", "gtaubygmu_K", "gtaubygmu_fit",
    "gtaubyge_tau", "gmubyge_tau"
    )
  print(rbind(cbind(val=quant$val()[display.names], err=quant$err()[display.names])))

  ##
  ## define LeTeX commands with elaborated values
  ##
  toTex = TrStr$new("0123456789_", "zothfvsneiU")
  specialFormat = c(
    tau_tau="%.1f"
    )
  specialFactor = c(
    tau_tau=1e15
    )
  rc = mapply(function(name, val, err) {
    fmt = ifelse(is.na(specialFormat[name]), "%.4f", specialFormat[name])
    val = ifelse(is.na(specialFactor[name]), val, val*specialFactor[name])
    err = ifelse(is.na(specialFactor[name]), err, err*specialFactor[name])
    name.orig = name
    name = gsub("(\\d+)", "N\\1N", name)
    fmt = paste("}{\\ensuremath{", fmt, " \\pm ", fmt, "}}", sep="")
    rc = paste("\\newcommand{\\", toTex$tr(name), sprintf(fmt, val, err), "% ", name.orig, sep="") 
  }, names(quant$val()), quant$val(), quant$err())
  fname.short = gsub("average[^-]*-*", "", file.name)
  fname.short = gsub("[.]rdata", "", fname.short)
  fname = "../report/tau-elab"
  if (fname.short != "") fname = paste(fname, "-", fname.short, sep="")
  if (flag.unitarity) fname = paste(fname, "-unitar", sep="")
  if (flag.vadirect) fname = paste(fname, "-vadirect", sep="")
  if (flag.lepuniv) fname = paste(fname, "-lepuniv", sep="")
  fname = paste(fname, ".tex", sep="")
  cat(rc, sep="\n", file=fname)
  cat("produced file '", fname, "'\n", sep="")
}

args = commandArgs(TRUE)
aluelab.results(args) 
