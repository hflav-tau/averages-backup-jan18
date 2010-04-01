#!/usr/bin/env Rscript

##
## print in standard output the Combos cards to combine the BaBar results on
## tau -> 3h nu BRs
##
## 1) PimPimPipNu
## 2) PimKmPipNu
## 3) PimKmKpNu
## 4) KmKmKpNu
##

name.1 = "PimPimPipNu"
name.2 = "PimKmPipNu"
name.3 = "PimKmKpNu"
name.4 = "KmKmKpNu"

BR.1 = 8.8337e-2;
BR.2 = 0.27257e-2; 
BR.3 = 0.13461e-2;
BR.4 = 1.5777e-5;

stat.1 = 0.0074e-2;
stat.2 = 0.0018e-2;
stat.3 = 0.001e-2;
stat.4 = 0.13e-5;

syst.1 = 0.126724e-2;
syst.2 = 0.00924406e-2;
syst.3 = 0.00364131e-2;
syst.4 = 0.12308e-5;

##
## total correlation coefficients
##
rhotot.12 = 0.543535;
rhotot.13 = 0.390346;
rhotot.14 = 0.031469;
rhotot.23 = 0.177495;
rhotot.24 = 0.0931907;
rhotot.34 = 0.0870484;

##-- print result
pr.result = function(br, stat, syst) {
  ##-- set stat. error to total error
  cat("           ", sprintf("%-11.5e   %-11.5e %-11.5e\n", br, sqrt(sum(c(stat,syst)^2)), 0))
  ##-- stat/syst errors as usual
  ## cat("           ", sprintf("%-11.5e   %-11.5e %-11.5e\n", br, stat, syst))
}

##
## begin BR 1
##
cat("* B(tau- -> pi- pi- pi+ nu) [ex K0]\n")
cat("BEGIN    BaBar PimPimPipNu published PRL100:011801,2008\n")
cat("\n")

##
## stat. error = sum of paper stat. and syst. errors
## syst. error = 0 since we neglect negligible common systematics
##
cat("MEASUREMENT m_PimPimPipNu statistical systematic\n")
cat("DATA        m_PimPimPipNu statistical systematic\n")
pr.result(BR.1, stat.1, syst.1)
cat("\n")

cat("STAT_CORR_WITH BaBar", name.2, "published", rhotot.12, "\n")
cat("STAT_CORR_WITH BaBar", name.3, "published", rhotot.13, "\n")
cat("STAT_CORR_WITH BaBar", name.4, "published", rhotot.14, "\n")
cat("\n")
cat("END\n")

cat("\n")
##
## begin BR 2
##
cat("* B(tau- -> pi- K- pi+ nu) [ex. K0]\n")
cat("BEGIN    BaBar PimKmPipNu published PRL100:011801,2008\n")
cat("\n")

##
## stat. error = sum of paper stat. and syst. errors
## syst. error = 0 since we neglect negligible common systematics
##
cat("MEASUREMENT m_PimKmPipNu statistical systematic\n")
cat("DATA        m_PimKmPipNu statistical systematic\n")
pr.result(BR.2, stat.2, syst.2)
cat("\n")

cat("STAT_CORR_WITH BaBar", name.1, "published", rhotot.12, "\n")
cat("STAT_CORR_WITH BaBar", name.3, "published", rhotot.23, "\n")
cat("STAT_CORR_WITH BaBar", name.4, "published", rhotot.24, "\n")
cat("\n")
cat("END\n")

cat("\n")
##
## begin BR 3
##
cat("* B(tau- -> pi- K- K+ nu)\n")
cat("BEGIN    BaBar PimKmKpNu published PRL100:011801,2008\n")
cat("\n")

##
## stat. error = sum of paper stat. and syst. errors
## syst. error = 0 since we neglect negligible common systematics
##
cat("MEASUREMENT m_PimKmKpNu statistical systematic\n")
cat("DATA        m_PimKmKpNu statistical systematic\n")
pr.result(BR.3, stat.3, syst.3)
cat("\n")

cat("STAT_CORR_WITH BaBar", name.1, "published", rhotot.13, "\n")
cat("STAT_CORR_WITH BaBar", name.2, "published", rhotot.23, "\n")
cat("STAT_CORR_WITH BaBar", name.4, "published", rhotot.34, "\n")
cat("\n")
cat("END\n")

cat("\n")
##
## begin BR 4
##
cat("* B(tau- -> K- K- K+ nu)\n")
cat("BEGIN    BaBar KmKmKpNu published PRL100:011801,2008\n")
cat("\n")

##
## stat. error = sum of paper stat. and syst. errors
## syst. error = 0 since we neglect negligible common systematics
##
cat("MEASUREMENT m_KmKmKpNu statistical systematic\n")
cat("DATA        m_KmKmKpNu statistical systematic\n")
pr.result(BR.4, stat.4, syst.4)
cat("\n")

cat("STAT_CORR_WITH BaBar", name.1, "published", rhotot.14, "\n")
cat("STAT_CORR_WITH BaBar", name.2, "published", rhotot.24, "\n")
cat("STAT_CORR_WITH BaBar", name.3, "published", rhotot.34, "\n")
cat("\n")
cat("END\n")
