#!/bin/csh
setenv INPUT $1
if ( $INPUT == 0 ) then
  setenv LINPUT unconstrained
else if ( $INPUT == 1 ) then
  setenv LINPUT unconstrained_aleph_hcorr
else if ( $INPUT == 2 ) then
  setenv LINPUT unconstrained_only_lepcorr
else if ( $INPUT == 3 ) then
  setenv LINPUT constrained
else if ( $INPUT == 4 ) then
  setenv LINPUT constrained_aleph_hcorr
endif

setenv LOGFILE readpdg_${LINPUT}.log
setenv RESFILE readpdg_${LINPUT}.results

rm -rf maketable.temp
mkdir  maketable.temp

grep -A42 "Final Results from fit with" NoBB/$LOGFILE >  maketable.temp/a1
grep -A42 "Final Results from fit with" $LOGFILE > maketable.temp/a2

grep -F "&" maketable.temp/a1 | sed -e 's/\\\\/ /' > maketable.temp/a1_1 
grep -F "&" maketable.temp/a2 | awk -F "&" '{print "& ",$2}' > maketable.temp/a2_1 

paste maketable.temp/a1_1 maketable.temp/a2_1 > maketable.temp/table.body

head -2 maketable.temp/table.body > maketable.temp/table_lep.body

head -25 maketable.temp/table.body | tail -23 > maketable.temp/table_ns.body

tail -15 maketable.temp/table.body > maketable.temp/table_s.body

setenv GammaStrVal `grep Gamma110 $LOGFILE | awk '{printf "%8.4f\n", $(NF-2)}'`
setenv GammaStrErr `grep Gamma110 $LOGFILE | awk '{printf "%8.4f\n", $(NF-0)}'`
#
setenv GammaStrValNoBB `grep Gamma110 NoBB/$LOGFILE | awk '{printf "%8.4f\n", $(NF-2)}'`
setenv GammaStrErrNoBB `grep Gamma110 NoBB/$LOGFILE | awk '{printf "%8.4f\n", $(NF-0)}'`
#
setenv GammaAllVal `grep GammaAll $LOGFILE | awk '{printf "%8.4f\n", 100.*$(NF-2)}'`
setenv GammaAllErr `grep GammaAll $LOGFILE | awk '{printf "%8.4f\n", 100.*$(NF-0)}'`
#
setenv GammaAllValNoBB `grep GammaAll NoBB/$LOGFILE | awk '{printf "%8.4f\n", 100.*$(NF-2)}'`
setenv GammaAllErrNoBB `grep GammaAll NoBB/$LOGFILE | awk '{printf "%8.4f\n", 100.*$(NF-0)}'`

echo "\\hline" > maketable.temp/table.tail
echo 'Sum of strange modes                                                   & '    ${GammaStrValNoBB}' $\pm$ '${GammaStrErrNoBB}' & '${GammaStrVal}' $\pm$ '${GammaStrErr}' \\\\ \\hline' >> maketable.temp/table.tail
echo 'Sum of all modes                                                       & '    ${GammaAllValNoBB}' $\pm$ '${GammaAllErrNoBB}' & '${GammaAllVal}' $\pm$ '${GammaAllErr}' \\\\ \\hline' >> maketable.temp/table.tail

setenv Be_Val_NoBB `grep -F "B(tau- -> e- nub nu) " NoBB/$RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Be_Err_NoBB `grep -F "B(tau- -> e- nub nu) " NoBB/$RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $3}'`

setenv Be_Val      `grep -F "B(tau- -> e- nub nu) " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Be_Err      `grep -F "B(tau- -> e- nub nu) " $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $3}'`

# echo 'Be                               & ' ${Be_Val_NoBB} ' $\pm$ ' ${Be_Err_NoBB} ' & ' ${Be_Val} ' $\pm$ ' ${Be_Err} ' \\\\' >> maketable.temp/table.tail

setenv Be_from_Bmu_Val_NoBB `grep -F "B(tau- -> e- nub nu)_Bmu" NoBB/$RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Be_from_Bmu_Err_NoBB `grep -F "B(tau- -> e- nub nu)_Bmu" NoBB/$RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $3}'`

setenv Be_from_Bmu_Val      `grep -F "B(tau- -> e- nub nu)_Bmu" $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Be_from_Bmu_Err      `grep -F "B(tau- -> e- nub nu)_Bmu" $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $3}'`


setenv Corr_Be_Bmu_NoBB `grep -F "Corr between Be_from_Bmu and Be " NoBB/$RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Corr_Be_Bmu      `grep -F "Corr between Be_from_Bmu and Be "      $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`

echo 'Be (from Bmu)  [ Corr. with Be ]                & ' ${Be_from_Bmu_Val_NoBB} ' $\pm$ ' ${Be_from_Bmu_Err_NoBB} '[' ${Corr_Be_Bmu_NoBB} '\%]' ' & ' ${Be_from_Bmu_Val} ' $\pm$ ' ${Be_from_Bmu_Err} '[' ${Corr_Be_Bmu} '\%]' '\\\\' >> maketable.temp/table.tail

# echo 'Corr (Be, Be (from Bmu))         & ' ${Corr_Be_Bmu_NoBB} '                & ' ${Corr_Be_Bmu} '                \\\\' >> maketable.temp/table.tail

setenv Be_from_tautau_Val_NoBB `grep -F "B(tau- -> e- nub nu)_tautau" NoBB/$RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Be_from_tautau_Err_NoBB `grep -F "B(tau- -> e- nub nu)_tautau" NoBB/$RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $3}'`

setenv Be_from_tautau_Val      `grep -F "B(tau- -> e- nub nu)_tautau" $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Be_from_tautau_Err      `grep -F "B(tau- -> e- nub nu)_tautau" $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $3}'`

setenv Corr_Be_Bett_NoBB `grep -F "Corr between Be_from_Bmu and Be_from_tautau" NoBB/$RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Corr_Be_Bett      `grep -F "Corr between Be_from_Bmu and Be_from_tautau"      $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`

echo 'Be (from  $\\tau_\\tau$)  [ Corr. with Be ]          & ' ${Be_from_tautau_Val_NoBB} ' $\pm$ ' ${Be_from_tautau_Err_NoBB} '[' ${Corr_Be_Bett_NoBB} '\%]' ' & ' ${Be_from_tautau_Val} ' $\pm$ ' ${Be_from_tautau_Err} '[' ${Corr_Be_Bett} '\%]' ' \\\\' >> maketable.temp/table.tail

# echo 'Corr (Be, Be (from $\\tau_\\tau$)) & ' ${Corr_Be_Bett_NoBB} '                & ' ${Corr_Be_Bett} '                \\\\' >> maketable.temp/table.tail

setenv Be_univ_Val_NoBB `grep -F "<B(tau- -> e- nub nu)_univ>" NoBB/$RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Be_univ_Err_NoBB `grep -F "<B(tau- -> e- nub nu)_univ>" NoBB/$RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $3}'`

setenv Be_univ_Val      `grep -F "<B(tau- -> e- nub nu)_univ>"      $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $1}'`
setenv Be_univ_Err      `grep -F "<B(tau- -> e- nub nu)_univ>"      $RESFILE | tail -1 | sed -e 's/(/ /g' | sed -e 's/)/ /g' | awk -F= '{print $2}' | awk '{print $3}'`

echo 'Be$_{univ}$                        & ' ${Be_univ_Val_NoBB} ' $\pm$ ' ${Be_univ_Err_NoBB} ' & ' ${Be_univ_Val} ' $\pm$ ' ${Be_univ_Err} ' \\\\' >> maketable.temp/table.tail

echo "\\hline\\end{tabular}}}" >> maketable.temp/table.tail
if ( $INPUT == 0 ) then
echo "\\caption{Results for branching fractions (in \%) from unconstrained fit to data from non-B-Factories and including those from B-Factories. PDG modification of ALEPH correlation matrix is used as input.}" >> maketable.temp/table.tail
echo "\\label{tab:TauGlobalFit_unconstrained}" >> maketable.temp/table.tail
else if ( $INPUT == 1 ) then
echo "\\caption{Results for branching fractions (in \%) from unconstrained fit to data from non-B-Factories and including those from B-Factories. }" >> maketable.temp/table.tail
echo "\\label{tab:TauGlobalFit_unconstrained_aleph_hcorr}" >> maketable.temp/table.tail
else if ( $INPUT == 2 ) then
echo "\\caption{Results for branching fractions (in \%) from unconstrained fit to data from non-B-Factories and including those from B-Factories. Correlations between the measured decay modes have been turned off, except those between the leptonic mdoes.}" >> maketable.temp/table.tail
echo "\\label{tab:TauGlobalFit_unconstrained_only_lepcorr}" >> maketable.temp/table.tail
else if ( $INPUT == 3 ) then
echo "\\caption{Results for branching fractions (in \%) from constrained fit to data from non-B-Factories and including those from B-Factories. PDG modification of ALEPH correlation matrix is used as input.}" >> maketable.temp/table.tail
echo "\\label{tab:TauGlobalFit_constrained}" >> maketable.temp/table.tail
else if ( $INPUT == 4 ) then
echo "\\caption{Results for branching fractions (in \%) from constrained fit to data from non-B-Factories and including those from B-Factories.}" >> maketable.temp/table.tail
echo "\\label{tab:TauGlobalFit_constrained_aleph_hcorr}" >> maketable.temp/table.tail
endif
echo "\\end{center}" >> maketable.temp/table.tail
echo "\\end{table}" >> maketable.temp/table.tail

cat maketable.head maketable.temp/table_lep.body maketable_ns.head maketable.temp/table_ns.body maketable_s.head maketable.temp/table_s.body maketable.temp/table.tail | tee readpdg_${LINPUT}_table.tex

