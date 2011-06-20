#!/bin/csh
#
# prints numbers from this *.log file and corresponding file in NoBB directory
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

rm -rf maketable1.temp
mkdir  maketable1.temp

grep -A42 "Final Results from fit with" NoBB/$LOGFILE >  maketable1.temp/a1
grep -A42 "Final Results from fit with" $LOGFILE > maketable1.temp/a2

grep -F "&" maketable1.temp/a1 | sed -e 's/\\\\/ /' > maketable1.temp/a1_1 
grep -F "&" maketable1.temp/a2 | awk -F "&" '{print "& ",$2}' > maketable1.temp/a2_1 

paste maketable1.temp/a1_1 maketable1.temp/a2_1 > maketable1.temp/table.body

head -2 maketable1.temp/table.body > maketable1.temp/table_lep.body

head -25 maketable1.temp/table.body | tail -23 > maketable1.temp/table_ns.body

tail -15 maketable1.temp/table.body > maketable1.temp/table_s.body

setenv GammaStrVal `grep Gamma110 $LOGFILE | awk '{printf "%8.4f\n", 100.*$(NF-2)}'`
setenv GammaStrErr `grep Gamma110 $LOGFILE | awk '{printf "%8.4f\n", 100.*$(NF-0)}'`
#
setenv GammaStrValNoBB `grep Gamma110 NoBB/$LOGFILE | awk '{printf "%8.4f\n", 100.*$(NF-2)}'`
setenv GammaStrErrNoBB `grep Gamma110 NoBB/$LOGFILE | awk '{printf "%8.4f\n", 100.*$(NF-0)}'`
#
setenv GammaAllVal `grep GammaAll $LOGFILE | awk '{printf "%8.4f\n", 100.*$(NF-2)}'`
setenv GammaAllErr `grep GammaAll $LOGFILE | awk '{printf "%8.4f\n", 100.*$(NF-0)}'`
#
setenv GammaAllValNoBB `grep GammaAll NoBB/$LOGFILE | awk '{printf "%8.4f\n", 100.*$(NF-2)}'`
setenv GammaAllErrNoBB `grep GammaAll NoBB/$LOGFILE | awk '{printf "%8.4f\n", 100.*$(NF-0)}'`

echo "\\hline" > maketable1.temp/table.tail
echo 'Sum of strange modes                                                   & '    ${GammaStrValNoBB}' $\pm$ '${GammaStrErrNoBB}' & '${GammaStrVal}' $\pm$ '${GammaStrErr}' \\\\ \\hline' >> maketable1.temp/table.tail
echo 'Sum of all modes                                                       & '    ${GammaAllValNoBB}' $\pm$ '${GammaAllErrNoBB}' & '${GammaAllVal}' $\pm$ '${GammaAllErr}' \\\\ \\hline' >> maketable1.temp/table.tail
#
setenv NMeas `grep chisquare_tot $LOGFILE | tail -2 | head -1 | awk -F= '{print $3}' | awk '{printf "%3d\n", $1}'`
setenv NBase `grep chisquare_tot $LOGFILE | tail -2 | head -1 | awk -F= '{print $4}' | awk '{printf "%3d\n", $1}'`
setenv Chi2  `grep chisquare_tot $LOGFILE | tail -2 | head -1 | awk -F= '{print $2}' | awk '{printf "%5.1f\n", $1}'`
setenv Prob  `grep chisquare_tot $LOGFILE | tail -2 | head -1  | awk -F= '{print $6}' | awk '{printf "%5.1f\n", $1*100}'`
#
setenv NMeas_NoBB `grep chisquare_tot NoBB/$LOGFILE | tail -2 | head -1 | awk -F= '{print $3}' | awk '{printf "%3d\n", $1}'`
setenv NBase_NoBB `grep chisquare_tot NoBB/$LOGFILE | tail -2 | head -1 | awk -F= '{print $4}' | awk '{printf "%3d\n", $1}'`
setenv Chi2_NoBB  `grep chisquare_tot NoBB/$LOGFILE | tail -2 | head -1 | awk -F= '{print $2}' | awk '{printf "%5.1f\n", $1}'`
setenv Prob_NoBB  `grep chisquare_tot NoBB/$LOGFILE | tail -2 | head -1 | awk -F= '{print $6}' | awk '{printf "%5.1f\n", $1*100}'`
#
#
echo "\\hline" >> maketable1.temp/table.tail
echo '\# of measurements                                                     & '    ${NMeas_NoBB}' & '${NMeas}' \\\\ \\hline' >> maketable1.temp/table.tail
echo '\# of base modes                                                       & '    ${NBase_NoBB}' & '${NBase}' \\\\ \\hline' >> maketable1.temp/table.tail
echo '$\\chi^2$                                                              & '    ${Chi2_NoBB}'  & '${Chi2}'  \\\\ \\hline' >> maketable1.temp/table.tail
echo 'CL (\\%)                                                               & '    ${Prob_NoBB}'  & '${Prob}'  \\\\ \\hline' >> maketable1.temp/table.tail
#
echo "\\hline\\end{tabular}}}" >> maketable1.temp/table.tail
if ( $INPUT == 0 ) then
echo "\\caption{Results for branching fractions (in \%) from unconstrained fit to data from non-B-Factories and including those from B-Factories. PDG modification of ALEPH correlation matrix is used as input.}" >> maketable1.temp/table.tail
echo "\\label{tab:TauGlobalFit_unconstrained_NoBB_BB}" >> maketable1.temp/table.tail
else if ( $INPUT == 1 ) then
echo "\\caption{Results for branching fractions (in \%) from unconstrained fit to data from non-B-Factories and including those from B-Factories. }" >> maketable1.temp/table.tail
echo "\\label{tab:TauGlobalFit_unconstrained_aleph_hcorr_NoBB_BB}" >> maketable1.temp/table.tail
else if ( $INPUT == 2 ) then
echo "\\caption{Results for branching fractions (in \%) from constrained fit to data from non-B-Factories and including those from B-Factories. PDG modification of ALEPH correlation matrix is used as input.}" >> maketable1.temp/table.tail
echo "\\label{tab:TauGlobalFit_constrained_NoBB_BB}" >> maketable1.temp/table.tail
else if ( $INPUT == 3 ) then
echo "\\caption{Results for branching fractions (in \%) from constrained fit to data from non-B-Factories and including those from B-Factories.}" >> maketable1.temp/table.tail
echo "\\label{tab:TauGlobalFit_constrained_aleph_hcorr_NoBB_BB}" >> maketable1.temp/table.tail
endif
echo "\\end{center}" >> maketable1.temp/table.tail
echo "\\end{table}" >> maketable1.temp/table.tail
echo "\\end{document}" >>  maketable1.temp/table.tail

cat maketable1.head maketable1.temp/table_lep.body maketable1_ns.head maketable1.temp/table_ns.body maketable1_s.head maketable1.temp/table_s.body maketable1.temp/table.tail | tee readpdg_${LINPUT}_table.tex
