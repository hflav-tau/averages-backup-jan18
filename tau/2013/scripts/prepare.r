#!/usr/bin/env Rscript

require(stringr, quietly=TRUE)
source("../../../Common/bin/aluelab3.r")

## ////////////////////////////////////////////////////////////////////////////
## functions

## ////////////////////////////////////////////////////////////////////////////
## code

quant = StatComb$new()

##--- constants used by Belle for arXiv:1402.5213
rc = quant$meas.add.single("Be", 0.1783, 0.0004)
rc = quant$meas.add.single("Bmu", 0.1741, 0.0004)

##
## Belle KS paper results
##
## some BR depend on e-mu normalization
## and have syst. contr. from Be and Bmu
## the BR expression has the following factor
##
quant$meas.expr.add("e_mu_norm", (Be * Bmu) / (Be + Bmu))

##--- rel. error due to normalization
e_mu_norm.rel.syst.perc =  quant$err()["e_mu_norm"] / quant$val()["e_mu_norm"] * 100

e_mu_norm.rel.syst.Be.perc = quant$syst.contrib.perc("e_mu_norm", "Be")
e_mu_norm.rel.syst.Bmu.perc = quant$syst.contrib.perc("e_mu_norm", "Bmu")

##
## Belle KS paper results
##
## correlated systematic contribution due to e-mu normalization
##

##--- from the paper
e_mu_norm.re.syst.total.perc = 0.5

##--- part of e-mu normalization syst. contribution not due to leptonic BRs
e_mu_norm.re.syst.nonBl.perc = sqrt(e_mu_norm.re.syst.total.perc^2 - e_mu_norm.rel.syst.perc^2)

##
## print
##

cat("#\n")
cat("# Belle paper on KS BRs\n")
cat("#\n")

cat("\n")
cat(sprintf("# total systematic contr. for e-mu normalization: %.4g%%\n", e_mu_norm.rel.syst.perc))
cat("SYSTEMATICS\n")
cat(sprintf(" %8.4g%%         EmNuNumb         # e-mu normalization, Be contr.\n",  e_mu_norm.rel.syst.Be.perc))
cat(sprintf(" %8.4g%%         MumNuNumb        # e-mu normalization, Bmu contr.\n",  e_mu_norm.rel.syst.Bmu.perc))
cat("\n")

cat("SYSTPAPER\n")
cat(sprintf(" %8.4g%%         e_mu_norm        # e-mu normalization, non Bl contr.\n",  e_mu_norm.re.syst.nonBl.perc))
cat("\n")

cat("PARAMETERS\n")
cat(sprintf("  EmNuNumb        %8.4g          +-%.4g\n",  quant$val("Be"), quant$err("Be")))
cat(sprintf("  MumNuNumb       %8.4g          +-%.4g\n",  quant$val("Bmu"), quant$err("Bmu")))

##
## BaBar KS KS paper
##
## a common systematic is quoted, which includes the tau cross-section and BaBar Luminosity
## we separate between tau cross-section and everything else
## we do not account for BaBar luminosity uncertainty correlations
##

syst.sigmataupmy4s = -0.3/100

babar.ksks.syst.common = 0.034
babar.kskspi0.syst.common = 0.03

babar.ksks.syst.common.noxs = sqrt(babar.ksks.syst.common^2 - syst.sigmataupmy4s^2)
babar.kskspi0.syst.common.noxs = sqrt(babar.kskspi0.syst.common^2 - syst.sigmataupmy4s^2)

cat("\n")

cat("#\n")
cat("# BaBar paper on KS KS BRs\n")
cat("#\n")

cat("\n")
cat("# split of common systematics for BaBar tau -> pi KS KS nu\n")
cat("SYSTEMATICS\n")
cat(sprintf(" %8.4g%%         sigmataupmy4s\n",  syst.sigmataupmy4s*100))

cat("\n")
cat("SYSTPAPER\n")
cat(sprintf(" %8.4g&         common_noxs      # common syst except for tau x-section\n",  babar.ksks.syst.common.noxs))

cat("\n")
cat("PARAMETERS\n")
cat(sprintf("  sigmataupmy4s  %8.4g           +-%.4g\n",  0.919, 0.003))

cat("\n")
cat("# split of common systematics for BaBar tau -> pi KS KS pi0 nu\n")
cat("SYSTEMATICS\n")
cat(sprintf(" %8.4g%%         sigmataupmy4s\n",  syst.sigmataupmy4s*100))

cat("\n")
cat("SYSTPAPER\n")
cat(sprintf(" %8.4g&         common_noxs      # common syst except for tau x-section\n",  babar.kskspi0.syst.common.noxs))

cat("\n")
cat("PARAMETERS\n")
cat(sprintf("  sigmataupmy4s  %8.4g           +-%.4g\n",  0.919, 0.003))
