#!/usr/bin/env Rscript

##
## print in standard output the Combos cards to combine the Belle results on
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

BR.1 = 8.42e-2
BR.2 = 3.30e-3
BR.3 = 1.55e-3
BR.4 = 3.29e-5

stat.1 = 0.01e-2
stat.2 = 0.01e-3
stat.3 = 0.01e-3
stat.4 = 0.17e-5

syst.1 = 0.255e-2
syst.2 = 0.165e-3
syst.3 = 0.055e-3
syst.4 = 0.195e-5

sigma.1 = sqrt(stat.1**2 + syst.1**2)
sigma.2 = sqrt(stat.2**2 + syst.2**2)
sigma.3 = sqrt(stat.3**2 + syst.3**2)
sigma.4 = sqrt(stat.4**2 + syst.4**2)

## The total covariance matrix:
covtot.11 = 6.698e-06
covtot.12 = 7.551e-08
covtot.13 = 7.167e-09
covtot.14 =-3.587e-10
covtot.21 = 7.551e-08
covtot.22 = 2.78e-08
covtot.23 = 7.49e-10
covtot.24 = 1.515e-11
covtot.31 = 7.167e-09
covtot.32 = 7.49e-10
covtot.33 = 3.132e-09
covtot.34 =-1.206e-12
covtot.41 =-3.587e-10
covtot.42 = 1.515e-11
covtot.43 =-1.206e-12
covtot.44 = 6.72e-12

## The total correlation coefficients between channels are:
rhotot.12 = covtot.12/sqrt(covtot.11*covtot.22)
rhotot.13 = covtot.13/sqrt(covtot.11*covtot.33)
rhotot.14 = covtot.14/sqrt(covtot.11*covtot.44)
rhotot.23 = covtot.23/sqrt(covtot.22*covtot.33)
rhotot.24 = covtot.24/sqrt(covtot.22*covtot.44)
rhotot.34 = covtot.34/sqrt(covtot.33*covtot.44)

## Update with higher precision numbers
sigma.1 = sqrt(covtot.11)
sigma.2 = sqrt(covtot.22)
sigma.3 = sqrt(covtot.33)
sigma.4 = sqrt(covtot.44)

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
cat("BEGIN    Belle PimPimPipNu published arXiv:1001.0083 (2010)\n")
cat("\n")

##
## stat. error = sum of paper stat. and syst. errors
## syst. error = 0 since we neglect negligible common systematics
##
cat("MEASUREMENT m_PimPimPipNu statistical systematic\n")
cat("DATA        m_PimPimPipNu statistical systematic\n")
pr.result(BR.1, stat.1, syst.1)
cat("\n")

cat("STAT_CORR_WITH Belle", name.2, "published", rhotot.12, "\n")
cat("STAT_CORR_WITH Belle", name.3, "published", rhotot.13, "\n")
cat("STAT_CORR_WITH Belle", name.4, "published", rhotot.14, "\n")
cat("\n")
cat("END\n")

cat("\n")
##
## begin BR 2
##
cat("* B(tau- -> pi- K- pi+ nu) [ex. K0]\n")
cat("BEGIN    Belle PimKmPipNu published arXiv:1001.0083 (2010)\n")
cat("\n")

##
## stat. error = sum of paper stat. and syst. errors
## syst. error = 0 since we neglect negligible common systematics
##
cat("MEASUREMENT m_PimKmPipNu statistical systematic\n")
cat("DATA        m_PimKmPipNu statistical systematic\n")
pr.result(BR.2, stat.2, syst.2)
cat("\n")

cat("STAT_CORR_WITH Belle", name.1, "published", rhotot.12, "\n")
cat("STAT_CORR_WITH Belle", name.3, "published", rhotot.23, "\n")
cat("STAT_CORR_WITH Belle", name.4, "published", rhotot.24, "\n")
cat("\n")
cat("END\n")

cat("\n")
##
## begin BR 3
##
cat("* B(tau- -> pi- K- K+ nu)\n")
cat("BEGIN    Belle PimKmKpNu published arXiv:1001.0083 (2010)\n")
cat("\n")

##
## stat. error = sum of paper stat. and syst. errors
## syst. error = 0 since we neglect negligible common systematics
##
cat("MEASUREMENT m_PimKmKpNu statistical systematic\n")
cat("DATA        m_PimKmKpNu statistical systematic\n")
pr.result(BR.3, stat.3, syst.3)
cat("\n")

cat("STAT_CORR_WITH Belle", name.1, "published", rhotot.13, "\n")
cat("STAT_CORR_WITH Belle", name.2, "published", rhotot.23, "\n")
cat("STAT_CORR_WITH Belle", name.4, "published", rhotot.34, "\n")
cat("\n")
cat("END\n")

cat("\n")
##
## begin BR 4
##
cat("* B(tau- -> K- K- K+ nu)\n")
cat("BEGIN    Belle KmKmKpNu published arXiv:1001.0083 (2010)\n")
cat("\n")

##
## stat. error = sum of paper stat. and syst. errors
## syst. error = 0 since we neglect negligible common systematics
##
cat("MEASUREMENT m_KmKmKpNu statistical systematic\n")
cat("DATA        m_KmKmKpNu statistical systematic\n")
pr.result(BR.4, stat.4, syst.4)
cat("\n")

cat("STAT_CORR_WITH Belle", name.1, "published", rhotot.14, "\n")
cat("STAT_CORR_WITH Belle", name.2, "published", rhotot.24, "\n")
cat("STAT_CORR_WITH Belle", name.3, "published", rhotot.34, "\n")
cat("\n")
cat("END\n")
