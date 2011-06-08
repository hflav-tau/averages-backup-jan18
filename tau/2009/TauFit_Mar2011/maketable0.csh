#!/bin/csh
#
# prints numbers from this *.results AND *.log files
#
setenv INPUT $1
if ( $INPUT == 0 ) then
  setenv LINPUT unconstrained
else if ( $INPUT == 1 ) then
  setenv LINPUT unconstrained_aleph_hcorr
else if ( $INPUT == 2 ) then
  setenv LINPUT constrained
else if ( $INPUT == 3 ) then
  setenv LINPUT constrained_aleph_hcorr
endif

setenv LOGFILE readpdg_${LINPUT}.log
setenv RESFILE readpdg_${LINPUT}.results
#
rm -rf maketable0.temp
mkdir  maketable0.temp
#
touch  maketable0.temp/table.tail
#
echo "\\begin{table}[\!hbtp] " >> maketable0.temp/table.tail
echo "\\begin{center} " >> maketable0.temp/table.tail
if ( $INPUT == 0 ) then
echo "\\caption{Results from unitarity unconstrained fit. PDG modification of ALEPH correlation matrix is used as input.}" >> maketable0.temp/table.tail
echo "\\label{tab:TauGlobalFit_unconstrained}" >> maketable0.temp/table.tail
else if ( $INPUT == 1 ) then
echo "\\caption{Results from unitarity unconstrained fit.}" >> maketable0.temp/table.tail
echo "\\label{tab:TauGlobalFit_unconstrained_aleph_hcorr}" >> maketable0.temp/table.tail
else if ( $INPUT == 2 ) then
echo "\\caption{Results from unitarity constrained fit. PDG modification of ALEPH correlation matrix is used as input.}" >> maketable0.temp/table.tail
echo "\\label{tab:TauGlobalFit_constrained}" >> maketable0.temp/table.tail
else if ( $INPUT == 3 ) then
echo "\\caption{Results from unitarity constrained fit.}" >> maketable0.temp/table.tail
echo "\\label{tab:TauGlobalFit_constrained_aleph_hcorr}" >> maketable0.temp/table.tail
endif
echo "{\\scalebox{.9}{\\begin{tabular}{|l@{ = }c|r@{ = }r|} \hline"  >> maketable0.temp/table.tail
echo '\\multicolumn{1}{|l}{$\\tau^-$ decay} & \multicolumn{1}{c|}{Branching Ratios (\\%)} & \\multicolumn{2}{c|}{Correlations and Weights (\\%)} \\\\ \\hline ' >> maketable0.temp/table.tail
#
setenv Be_Val      `grep -F "B(tau- -> e- nub nu) = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Be_Err      `grep -F "B(tau- -> e- nub nu) = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $3}'`
setenv Corr_Be_Bmu `grep -F "Corr between e- nubar(e) nu(tau) and mu- nubar(mu) nu(tau) = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Corr_Be_Bpi `grep -F "Corr between e- nubar(e) nu(tau) and pi- nu(tau) = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Corr_Be_BK  `grep -F "Corr between e- nubar(e) nu(tau) and K- nu(tau) = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Corr_Be_Bs  `grep -F "Corr between e- nubar(e) nu(tau) and strange nu(tau) = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Corr_Be_Bt  `grep -F "Corr between e- nubar(e) nu(tau) and sum_of_all_decays = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
#
setenv Bmu_Val      `grep -F "B(tau- -> mu- nub nu) = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Bmu_Err      `grep -F "B(tau- -> mu- nub nu) = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $3}'`
setenv Corr_Bmu_Bpi `grep -F "Corr between mu- nubar(mu) nu(tau) and pi- nu(tau) = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Corr_Bmu_BK  `grep -F "Corr between mu- nubar(mu) nu(tau) and K- nu(tau) = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Corr_Bmu_Bs  `grep -F "Corr between mu- nubar(mu) nu(tau) and strange nu(tau) = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Corr_Bmu_Bt  `grep -F "Corr between mu- nubar(mu) nu(tau) and sum_of_all_decays = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
#
setenv Bpi_Val      `grep -F "B(tau- -> pi- nu) = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Bpi_Err      `grep -F "B(tau- -> pi- nu) = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $3}'`
setenv Corr_Bpi_BK  `grep -F "Corr between pi- nu(tau) and K- nu(tau) = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Corr_Bpi_Bs  `grep -F "Corr between pi- nu(tau) and strange nu(tau) = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Corr_Bpi_Bt  `grep -F "Corr between pi- nu(tau) and sum_of_all_decays = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
#
setenv BK_Val      `grep -F "B(tau- -> K- nu) = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv BK_Err      `grep -F "B(tau- -> K- nu) = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $3}'`
setenv Corr_BK_Bs  `grep -F "Corr between K- nu(tau) and strange nu(tau) = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Corr_BK_Bt  `grep -F "Corr between K- nu(tau) and sum_of_all_decays = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
#
setenv Bs_Val      `grep -F "B(tau -> strange) = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Bs_Err      `grep -F "B(tau -> strange) = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $3}'`
setenv Corr_Bs_Bt  `grep -F "Corr between strange nu(tau) and sum_of_all_decays" $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
#
setenv Bt_Val      `grep -F "B(tau -> total) = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Bt_Err      `grep -F "B(tau -> total) = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $3}'`
#
setenv Be_from_Bmu_Val      `grep -F "B(tau- -> e- nub nu)_Bmu  = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Be_from_Bmu_Err      `grep -F "B(tau- -> e- nub nu)_Bmu  = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $3}'`
setenv Corr_Be_from_Bmu_Be  `grep -F "Corr between Be_from_Bmu and Be             = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Corr_Be_from_Bmu_tt  `grep -F "Corr between Be_from_Bmu and Be_from_tautau = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Corr_Be_from_Bmu_Bpi `grep -F "Corr between Be_from_Bpi and Be_from_Bmu    = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Corr_Be_from_Bmu_BK  `grep -F "Corr between Be_from_BK  and Be_from_Bmu    = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
#
setenv Be_from_tt_Val       `grep -F "B(tau- -> e- nub nu)_tautau  = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Be_from_tt_Err       `grep -F "B(tau- -> e- nub nu)_tautau  = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $3}'`
setenv Corr_Be_from_tt_Be   "0.00"
setenv Corr_Be_from_tt_Bpi  `grep -F "Corr between Be_from_Bpi and Be_from_tautau = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Corr_Be_from_tt_BK   `grep -F "Corr between Be_from_BK  and Be_from_tautau = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
#
setenv Be_from_Bpi_Val       `grep -F "B(tau- -> e- nub nu)_Bpi  = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Be_from_Bpi_Err       `grep -F "B(tau- -> e- nub nu)_Bpi  = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $3}'`
setenv Corr_Be_from_Bpi_Be   `grep -F "Corr between Be_from_Bpi and Be             = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Corr_Be_from_Bpi_BK   `grep -F "Corr between Be_from_BK  and Be_from_Bpi    = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
#
setenv Be_from_BK_Val       `grep -F "B(tau- -> e- nub nu)_BK   = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Be_from_BK_Err       `grep -F "B(tau- -> e- nub nu)_BK   = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $3}'`
setenv Corr_Be_from_BK_Be   `grep -F "Corr between Be_from_BK  and Be             = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
#
setenv Be_univ3_Val         `grep -F "<B(tau- -> e- nub nu)_univ3> = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Be_univ3_Err         `grep -F "<B(tau- -> e- nub nu)_univ3> = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $3}'`
setenv Be_univ3_wt_Be       `grep -F "<B(tau- -> e- nub nu)_univ3> = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $3}' | awk '{print 100*$1}'`
setenv Be_univ3_wt_Bmu      `grep -F "<B(tau- -> e- nub nu)_univ3> = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $3}' | awk '{print 100*$3}'`
setenv Be_univ3_wt_tt       `grep -F "<B(tau- -> e- nub nu)_univ3> = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $3}' | awk '{print 100*$5}'`
#
setenv Be_univ4_Val         `grep -F "<B(tau- -> e- nub nu)_univ4> = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Be_univ4_Err         `grep -F "<B(tau- -> e- nub nu)_univ4> = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $3}'`
setenv Be_univ4_wt_Be       `grep -F "<B(tau- -> e- nub nu)_univ4> = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $3}' | awk '{print 100*$1}'`
setenv Be_univ4_wt_Bmu      `grep -F "<B(tau- -> e- nub nu)_univ4> = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $3}' | awk '{print 100*$3}'`
setenv Be_univ4_wt_tt       `grep -F "<B(tau- -> e- nub nu)_univ4> = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $3}' | awk '{print 100*$5}'`
setenv Be_univ4_wt_Bpi      `grep -F "<B(tau- -> e- nub nu)_univ4> = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $3}' | awk '{print 100*$7}'`
#
setenv Be_univ5_Val         `grep -F "<B(tau- -> e- nub nu)_univ5> = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Be_univ5_Err         `grep -F "<B(tau- -> e- nub nu)_univ5> = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $3}'`
setenv Be_univ5_wt_Be       `grep -F "<B(tau- -> e- nub nu)_univ5> = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $3}' | awk '{print 100*$1}'`
setenv Be_univ5_wt_Bmu      `grep -F "<B(tau- -> e- nub nu)_univ5> = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $3}' | awk '{print 100*$3}'`
setenv Be_univ5_wt_tt       `grep -F "<B(tau- -> e- nub nu)_univ5> = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $3}' | awk '{print 100*$5}'`
setenv Be_univ5_wt_Bpi      `grep -F "<B(tau- -> e- nub nu)_univ5> = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $3}' | awk '{print 100*$7}'`
setenv Be_univ5_wt_BK       `grep -F "<B(tau- -> e- nub nu)_univ5> = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $3}' | awk '{print 100*$9}'`
#
setenv Bhad_Val `grep -F "B(tau -> hadrons) = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Bhad_Err `grep -F "B(tau -> hadrons) = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $3}'`
setenv Rhad_Val `grep -F "R(tau -> hadrons) = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Rhad_Err `grep -F "R(tau -> hadrons) = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $3}'`
setenv Rs_Val  `grep -F "R(tau -> strange) = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Rs_Err  `grep -F "R(tau -> strange) = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $3}'`
setenv Rns_Val `grep -F "R(tau -> non-str) = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Rns_Err `grep -F "R(tau -> non-str) = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $3}'`
setenv Vus_Val  `grep -F "|Vus|_strange = " $RESFILE | grep -v Difference | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Vus_Err  `grep -F "|Vus|_strange = " $RESFILE | grep -v Difference | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $3}'`
#
setenv BmuBe_Val  `grep -F "B(tau- -> mu- nub nu)/B(tau- -> e- nub nu) = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv BmuBe_Err  `grep -F "B(tau- -> mu- nub nu)/B(tau- -> e- nub nu) = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $3}'`
setenv fmufe_Val  `grep -F "f(m_mu^2/m_tau^2)/f(m_e^2/m_tau^2) = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv fmufe_Err  `grep -F "f(m_mu^2/m_tau^2)/f(m_e^2/m_tau^2) = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $3}'`
setenv gmuge_Val  `grep -F "gmu/ge = " $RESFILE | grep -v Corr | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv gmuge_Err  `grep -F "gmu/ge = " $RESFILE | grep -v Corr | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $3}'`
#
setenv Bpiuniv_Val `grep -F "B(tau- -> pi- nu)_univ = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Bpiuniv_Err `grep -F "B(tau- -> pi- nu)_univ = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $3}'`
setenv gtaugmu_pi_Val  `grep -F "(gtau/gmu)_pi = " $RESFILE | grep -v Corr | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv gtaugmu_pi_Err  `grep -F "(gtau/gmu)_pi = " $RESFILE | grep -v Corr | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $3}'`
setenv Vus_KPi_Val `grep -F "|Vus|_TauToK/Pi = " $RESFILE | grep -v Diff | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Vus_KPi_Err `grep -F "|Vus|_TauToK/Pi = " $RESFILE | grep -v Diff | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $3}'`
#
setenv BKuniv_Val     `grep -F "B(tau- -> K- nu)_univ = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv BKuniv_Err     `grep -F "B(tau- -> K- nu)_univ = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $3}'`
setenv gtaugmu_K_Val  `grep -F "(gtau/gmu)_K = " $RESFILE | grep -v Corr | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv gtaugmu_K_Err  `grep -F "(gtau/gmu)_K = " $RESFILE | grep -v Corr | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $3}'`
setenv BKBpi_Val      `grep -F "G(K- nu(tau)) / G(pi- nu(tau)) = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv BKBpi_Err      `grep -F "G(K- nu(tau)) / G(pi- nu(tau)) = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $3}'`
setenv Vus_K_Val      `grep -F "|Vus|_TauToKmNu = " $RESFILE | grep -v Diff | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Vus_K_Err      `grep -F "|Vus|_TauToKmNu = " $RESFILE | grep -v Diff | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $3}'`
#
setenv gtaugmu_piK_Val  `grep -F "<(gtau/gmu)_pik> = " $RESFILE | grep -v Corr | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv gtaugmu_piK_Err  `grep -F "<(gtau/gmu)_pik> = " $RESFILE | grep -v Corr | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $3}'`
setenv gtaugmu_piK_piWt `grep -F "<(gtau/gmu)_pik> = " $RESFILE | grep -v Corr | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $3}' | awk '{print 100.*$1}'`
setenv gtaugmu_piK_KWt  `grep -F "<(gtau/gmu)_pik> = " $RESFILE | grep -v Corr | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $3}' | awk '{print 100.*$3}'`
setenv gtaugmu_piK_corr `grep -F "Corr between (gtau/gmu)_pi and (gtau/gmu)_K = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
#
setenv gtaugmu_piKe_Val  `grep -F "<(gtau/gmu)_pike> = " $RESFILE | grep -v Corr | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv gtaugmu_piKe_Err  `grep -F "<(gtau/gmu)_pike> = " $RESFILE | grep -v Corr | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $3}'`
setenv gtaugmu_piKe_piWt `grep -F "<(gtau/gmu)_pike> = " $RESFILE | grep -v Corr | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $3}' | awk '{print 100.*$1}'`
setenv gtaugmu_piKe_KWt  `grep -F "<(gtau/gmu)_pike> = " $RESFILE | grep -v Corr | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $3}' | awk '{print 100.*$3}'`
setenv gtaugmu_piKe_eWt  `grep -F "<(gtau/gmu)_pike> = " $RESFILE | grep -v Corr | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $3}' | awk '{print 100.*$5}'`
setenv gtaugmu_piKe_corrpi `grep -F "Corr between (gtau/gmu)_pi and (gtau/gmu)_e = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv gtaugmu_piKe_corrK  `grep -F "Corr between (gtau/gmu)_K  and (gtau/gmu)_e = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
#
setenv gtaugmu_e_Val   `grep -F "(gtau/gmu)_e = " $RESFILE | grep -v Corr | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv gtaugmu_e_Err   `grep -F "(gtau/gmu)_e = " $RESFILE | grep -v Corr | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $3}'`
setenv gtauge_mu_Val   `grep -F "(gtau/ge)_mu = " $RESFILE | grep -v Corr | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv gtauge_mu_Err   `grep -F "(gtau/ge)_mu = " $RESFILE | grep -v Corr | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $3}'`
setenv Bmu_from_tt_Val `grep -F "B(tau- -> mu- nub nu)_tautau  = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Bmu_from_tt_Err `grep -F "B(tau- -> mu- nub nu)_tautau  = " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $3}'`
#
echo '$\\mathrm{B_e}$                       & ' ${Be_Val} ' $\pm$ ' ${Be_Err} ' & [ Corr. with $\\mathrm{B_\\mu}$, $\\mathrm{B_\\pi}$, $\\mathrm{B_K}$, $\\mathrm{B_s}$, $\\mathrm{B_{tot}}$ ] & [' ${Corr_Be_Bmu} ',' ${Corr_Be_Bpi} ',' ${Corr_Be_BK} ',' ${Corr_Be_Bs} ',' ${Corr_Be_Bt} '] \\\\' >> maketable0.temp/table.tail
#
echo '$\\mathrm{B_\\mu}$                   & ' ${Bmu_Val} ' $\pm$ ' ${Bmu_Err}' & [ Corr. with    $\\mathrm{B_\\pi}$, $\\mathrm{B_K}$, $\\mathrm{B_s}$, $\\mathrm{B_{tot}}$ ] & [' ${Corr_Bmu_Bpi} ',' ${Corr_Bmu_BK} ',' ${Corr_Bmu_Bs} ',' ${Corr_Bmu_Bt} '] \\\\' >> maketable0.temp/table.tail
#
echo '$\\mathrm{B_\\pi}$                   & ' ${Bpi_Val} ' $\pm$ ' ${Bpi_Err}' & [ Corr. with       $\\mathrm{B_K}$, $\\mathrm{B_{strange}}$, $\\mathrm{B_{total}}$ ] & [' ${Corr_Bpi_BK} ',' ${Corr_Bpi_Bs} ',' ${Corr_Bpi_Bt} '] \\\\' >> maketable0.temp/table.tail
#
echo '$\\mathrm{B_K}$                    & ' ${BK_Val} ' $\pm$ ' ${BK_Err}' & [ Corr. with       $\\mathrm{B_{strange}}$, $\\mathrm{B_{total}}$ ] & [' ${Corr_BK_Bs} ',' ${Corr_BK_Bt} '] \\\\' >> maketable0.temp/table.tail
#
echo '$\\mathrm{B_{strange}}$              & ' ${Bs_Val} ' $\pm$ ' ${Bs_Err}' & [ Corr. with      $\\mathrm{B_{total}}$ ] & [' ${Corr_Bs_Bt} '] \\\\ \\cline{3-4}' >> maketable0.temp/table.tail
#
echo '$\\mathrm{B_{total}}$               & ' ${Bt_Val} ' $\pm$ ' ${Bt_Err}' & \multicolumn{2}{c|}{[Note: this is used to calculate $\\mathrm{B_{non-strange}}$ = $\\mathrm{B_{total}}$ - $\\mathrm{B_{strange}}$ - $\\mathrm{B_{leptonic}}$]}\\\\ \\hline' >> maketable0.temp/table.tail 
#
echo '$\\mathrm{B_e}$ (from $\\mathrm{B_\\mu}$)         & ' ${Be_from_Bmu_Val} ' $\pm$ ' ${Be_from_Bmu_Err}' & [ Corr. with $\\mathrm{B_e}$ from ($\\mathrm{B_e}$, $\\tau_\\tau$, $\\mathrm{B_\\pi}$, $\\mathrm{B_K}$)] & [' ${Corr_Be_from_Bmu_Be} ',' ${Corr_Be_from_Bmu_tt} ',' ${Corr_Be_from_Bmu_Bpi} ',' ${Corr_Be_from_Bmu_BK} '] \\\\' >> maketable0.temp/table.tail
#
echo '$\\mathrm{B_e}$ (from $\\tau_\\tau$)      & ' ${Be_from_tt_Val} ' $\pm$ ' ${Be_from_tt_Err}' & [ Corr. with $\\mathrm{B_e}$ from ($\\mathrm{B_e}$, $\\mathrm{B_\\pi}$, $\\mathrm{B_K}$)] & [' ${Corr_Be_from_tt_Be} ',' ${Corr_Be_from_tt_Bpi} ',' ${Corr_Be_from_tt_BK} '] \\\\' >> maketable0.temp/table.tail
#
echo '$\\mathrm{B_e}$ (from $\\mathrm{B_\\pi}$)      & ' ${Be_from_Bpi_Val} ' $\pm$ ' ${Be_from_Bpi_Err}' & [ Corr. with $\\mathrm{B_e}$ from ($\\mathrm{B_e}$, $\\mathrm{B_K}$)] & [' ${Corr_Be_from_Bpi_Be} ',' ${Corr_Be_from_Bpi_BK} '] \\\\' >> maketable0.temp/table.tail
#
echo '$\\mathrm{B_e}$ (from $\\mathrm{B_K}$)      & ' ${Be_from_BK_Val} ' $\pm$ ' ${Be_from_BK_Err}' & [ Corr. with $\\mathrm{B_e}$ from ($\\mathrm{B_e}$)] & [' ${Corr_Be_from_BK_Be} '] \\\\ \\hline ' >> maketable0.temp/table.tail
#
echo '$ \\langle \\mathrm{B_e} \\rangle_{\mathrm{(univ3)}} $  & ' ${Be_univ3_Val} ' $\pm$ ' ${Be_univ3_Err}' & [ Weights for $\\mathrm{B_e}$ from ($\\mathrm{B_e}$, $\\mathrm{B_\\mu}$, $\\tau_\\tau$)] & [' ${Be_univ3_wt_Be} ',' ${Be_univ3_wt_Bmu} ',' ${Be_univ3_wt_tt} '] \\\\' >> maketable0.temp/table.tail
#
echo '$ \\langle \\mathrm{B_e} \\rangle_{\mathrm{(univ4)}} $    & ' ${Be_univ4_Val} ' $\pm$ ' ${Be_univ4_Err}' & [ Weights for $\\mathrm{B_e}$ from ($\\mathrm{B_e}$, $\\mathrm{B_\\mu}$, $\\tau_\\tau$, $\\mathrm{B_\\pi}$)] & [' ${Be_univ4_wt_Be} ',' ${Be_univ4_wt_Bmu} ',' ${Be_univ4_wt_tt} ',' ${Be_univ4_wt_Bpi} '] \\\\' >> maketable0.temp/table.tail
#
echo '$ \\langle \\mathrm{B_e} \\rangle_{\mathrm{(univ5)}}$      & ' ${Be_univ5_Val} ' $\pm$ ' ${Be_univ5_Err}' & [ Weights for $\\mathrm{B_e}$ from ($\\mathrm{B_e}$, $\\mathrm{B_\\mu}$, $\\tau_\\tau$, $\\mathrm{B_\\pi}$, $\\mathrm{B_K}$)] & [' ${Be_univ5_wt_Be} ',' ${Be_univ5_wt_Bmu} ',' ${Be_univ5_wt_tt} ',' ${Be_univ5_wt_Bpi} ',' ${Be_univ5_wt_BK} '] \\\\ \\hline' >> maketable0.temp/table.tail
#
echo '$\mathrm{B_{had}}$ & \multicolumn{1}{c}{' ${Bhad_Val} ' $\pm$ ' ${Bhad_Err} '} &  \multicolumn{1}{c|}{$\mathrm{B_{had}}/ \\langle \\mathrm{B_e} \\rangle_{\mathrm{univ5}} $ = ' ${Rhad_Val} ' $\pm$ ' ${Rhad_Err} ' } & $|V_{us}|$ = ' ${Vus_Val} ' $\pm$ ' ${Vus_Err} '\\\\ \\cline{4-4}' >> maketable0.temp/table.tail
#
echo '$\mathrm{B_{strange}}/ \\langle \\mathrm{B_e} \\rangle_{\mathrm{univ5}}$ & \multicolumn{1}{c}{' ${Rs_Val} ' $\pm$ ' ${Rs_Err} '} &  \multicolumn{1}{c|}{$\mathrm{B_{non-strange}}/ \\langle \\mathrm{B_e} \\rangle_{\mathrm{univ5}}$ = ' ${Rns_Val} ' $\pm$ ' ${Rns_Err} ' } & $\\mathrm{B_K/B_\\pi}$ = ' ${BKBpi_Val} '$\pm$' ${BKBpi_Err} ' \\\\ \\cline{1-3}' >> maketable0.temp/table.tail
#
echo '$(\\mathrm{B_\\pi})_{univ}$ & \multicolumn{1}{c}{' ${Bpiuniv_Val} ' $\pm$ ' ${Bpiuniv_Err} '} &  \multicolumn{1}{c|}{$\mathrm{(g_\\tau/g_\\mu)_\\pi} $  = ' ${gtaugmu_pi_Val} ' $\pm$ ' ${gtaugmu_pi_Err} ' } & $|V_{us}|_{\mathrm{B_K}/{\mathrm{B_\\pi}}}$ = ' ${Vus_KPi_Val} '$\pm$' ${Vus_KPi_Err} ' \\\\ \\cline{1-3}' >> maketable0.temp/table.tail
#
echo '$(\\mathrm{B_K})_{univ}$  & \multicolumn{1}{c}{' ${BKuniv_Val} ' $\pm$ ' ${BKuniv_Err} '} &  \multicolumn{1}{c|}{ $\mathrm{(g_\\tau/g_\\mu)_K} $  = ' ${gtaugmu_K_Val} ' $\pm$ ' ${gtaugmu_K_Err} ' } & $|V_{us}|_{\mathrm{B_K}}$ = ' ${Vus_K_Val} '$\pm$' ${Vus_K_Err} ' \\\\ \\hline' >> maketable0.temp/table.tail
#
echo '$\\mathrm{B_\\mu}$ (from $\\tau_\\tau$)      & \multicolumn{1}{c}{' ${Bmu_from_tt_Val} ' $\pm$ ' ${Bmu_from_tt_Err} '} &  \multicolumn{1}{c|}{$\mathrm{(g_\\tau/g_e)_\\mu}$  = ' ${gtauge_mu_Val} ' $\pm$ ' ${gtauge_mu_Err} ' } & $\mathrm{(g_\\tau/g_\\mu)_e}$  = ' ${gtaugmu_e_Val} ' $\pm$ ' ${gtaugmu_e_Err} ' \\\\ \\hline' >> maketable0.temp/table.tail
#
echo '$\mathrm{\\langle(g_\\tau/g_\\mu)\\rangle_{\\pi,K}}$  & ' ${gtaugmu_piK_Val} ' $\pm$ ' ${gtaugmu_piK_Err} '&  \multicolumn{1}{c|}{[Weights for   $\mathrm{()_\\pi} $, $\mathrm{()_K} $ ] = [ ' ${gtaugmu_piK_piWt} ',' ${gtaugmu_piK_KWt} ']} &  Corr. [ $\mathrm{()_{\\pi K}} $] = [ ' ${gtaugmu_piK_corr} ' ]\\\\' >> maketable0.temp/table.tail
#
echo '$\mathrm{\\langle(g_\\tau/g_\\mu)\\rangle_{\\pi,K,e}}$  & ' ${gtaugmu_piKe_Val} ' $\pm$ ' ${gtaugmu_piKe_Err} ' & \multicolumn{1}{c|}{[Weights for   $\mathrm{()_\\pi} $, $\mathrm{()_K} $ , $\mathrm{()_e} $ ] = [ ' ${gtaugmu_piKe_piWt} ',' ${gtaugmu_piKe_KWt} ',' ${gtaugmu_piKe_eWt} '] }& Corr. [ $\mathrm{()_{\\pi e}} $, $\mathrm{()_{Ke}} $] = [ ' ${gtaugmu_piKe_corrpi} ',' ${gtaugmu_piKe_corrK} ' ]\\\\ \\hline' >> maketable0.temp/table.tail
#
echo '$\mathrm{(B_\\mu/B_e)}$ & \multicolumn{1}{c}{' ${BmuBe_Val} ' $\pm$ ' ${BmuBe_Err} '} & \multicolumn{1}{c}{$\mathrm{(f_\\mu/f_e)} $ = ' ${fmufe_Val} ' $\pm$ ' ${fmufe_Err} '} & $\mathrm{(g_\\mu/g_e)} $ = ' ${gmuge_Val} ' $\pm$ ' ${gmuge_Err} ' \\\\' >> maketable0.temp/table.tail
#
echo "\\hline\\end{tabular}}}" >> maketable0.temp/table.tail
echo "\\end{center}" >> maketable0.temp/table.tail
echo "\\end{table}" >> maketable0.temp/table.tail

# cat maketable0.head 
cat maketable0.temp/table.tail | tee readpdg_${LINPUT}_table0.tex

