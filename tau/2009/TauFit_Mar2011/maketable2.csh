#!/bin/csh
#
# prints numbers from this *.log file and corresponding file in BB_PiKUniv/ directory
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

rm -rf maketable2.temp
mkdir  maketable2.temp

grep -A42 "Final Results from fit with" BB_PiKUniv/$LOGFILE >  maketable2.temp/a1
grep -A42 "Final Results from fit with" $LOGFILE > maketable2.temp/a2

grep -F "&" maketable2.temp/a1 | sed -e 's/\\\\/ /' > maketable2.temp/a1_1 
grep -F "&" maketable2.temp/a2 | awk -F "&" '{print "& ",$2}' > maketable2.temp/a2_1 

paste maketable2.temp/a1_1 maketable2.temp/a2_1 > maketable2.temp/table.body

head -2 maketable2.temp/table.body > maketable2.temp/table_lep.body

head -25 maketable2.temp/table.body | tail -23 > maketable2.temp/table_ns.body

tail -15 maketable2.temp/table.body > maketable2.temp/table_s.body

setenv GammaStrVal `grep Gamma110 $LOGFILE | awk '{printf "%8.4f\n", 100.*$(NF-2)}'`
setenv GammaStrErr `grep Gamma110 $LOGFILE | awk '{printf "%8.4f\n", 100.*$(NF-0)}'`
#
setenv GammaStrValBBUniv `grep Gamma110 BB_PiKUniv/$LOGFILE | awk '{printf "%8.4f\n", 100.*$(NF-2)}'`
setenv GammaStrErrBBUniv `grep Gamma110 BB_PiKUniv/$LOGFILE | awk '{printf "%8.4f\n", 100.*$(NF-0)}'`
#
setenv GammaAllVal `grep GammaAll $LOGFILE | awk '{printf "%8.4f\n", 100.*$(NF-2)}'`
setenv GammaAllErr `grep GammaAll $LOGFILE | awk '{printf "%8.4f\n", 100.*$(NF-0)}'`
#
setenv GammaAllValBBUniv `grep GammaAll BB_PiKUniv/$LOGFILE | awk '{printf "%8.4f\n", 100.*$(NF-2)}'`
setenv GammaAllErrBBUniv `grep GammaAll BB_PiKUniv/$LOGFILE | awk '{printf "%8.4f\n", 100.*$(NF-0)}'`

echo "\\hline" > maketable2.temp/table.tail
echo 'Sum of strange modes                                                   & '    ${GammaStrValBBUniv}' $\pm$ '${GammaStrErrBBUniv}' & '${GammaStrVal}' $\pm$ '${GammaStrErr}' \\\\ \\hline' >> maketable2.temp/table.tail
echo 'Sum of all modes                                                       & '    ${GammaAllValBBUniv}' $\pm$ '${GammaAllErrBBUniv}' & '${GammaAllVal}' $\pm$ '${GammaAllErr}' \\\\ \\hline' >> maketable2.temp/table.tail
#
setenv NMeas `grep chisquare_tot $LOGFILE | tail -2 | head -1 | awk -F= '{print $3}' | awk '{printf "%3d\n", $1}'`
setenv NBase `grep chisquare_tot $LOGFILE | tail -2 | head -1 | awk -F= '{print $4}' | awk '{printf "%3d\n", $1}'`
setenv Chi2  `grep chisquare_tot $LOGFILE | tail -2 | head -1 | awk -F= '{print $2}' | awk '{printf "%5.1f\n", $1}'`
setenv Prob  `grep chisquare_tot $LOGFILE | tail -2 | head -1  | awk -F= '{print $6}' | awk '{printf "%5.1f\n", $1*100}'`
#
setenv NMeas_BBUniv `grep chisquare_tot BB_PiKUniv/$LOGFILE | tail -2 | head -1 | awk -F= '{print $3}' | awk '{printf "%3d\n", $1}'`
setenv NBase_BBUniv `grep chisquare_tot BB_PiKUniv/$LOGFILE | tail -2 | head -1 | awk -F= '{print $4}' | awk '{printf "%3d\n", $1}'`
setenv Chi2_BBUniv  `grep chisquare_tot BB_PiKUniv/$LOGFILE | tail -2 | head -1 | awk -F= '{print $2}' | awk '{printf "%5.1f\n", $1}'`
setenv Prob_BBUniv  `grep chisquare_tot BB_PiKUniv/$LOGFILE | tail -2 | head -1 | awk -F= '{print $6}' | awk '{printf "%5.1f\n", $1*100}'`
#
#
echo "\\hline" >> maketable2.temp/table.tail
echo '\# of measurements                                                     & '    ${NMeas_BBUniv}' & '${NMeas}' \\\\ \\hline' >> maketable2.temp/table.tail
echo '\# of base modes                                                       & '    ${NBase_BBUniv}' & '${NBase}' \\\\ \\hline' >> maketable2.temp/table.tail
echo '$\\chi^2$                                                              & '    ${Chi2_BBUniv}'  & '${Chi2}'  \\\\ \\hline' >> maketable2.temp/table.tail
echo 'CL (\\%)                                                               & '    ${Prob_BBUniv}'  & '${Prob}'  \\\\ \\hline' >> maketable2.temp/table.tail
#
echo "\\hline\\end{tabular}}}" >> maketable2.temp/table.tail
if ( $INPUT == 0 ) then
echo "\\caption{Results for branching fractions (in \%) from unconstrained fit to data from non-B-Factories and including those from B-Factories. PDG modification of ALEPH correlation matrix is used as input.}" >> maketable2.temp/table.tail
echo "\\label{tab:TauGlobalFit_unconstrained_BBUniv_BB}" >> maketable2.temp/table.tail
else if ( $INPUT == 1 ) then
echo "\\caption{Results for branching fractions (in \%) from unconstrained fit to data from non-B-Factories and including those from B-Factories. }" >> maketable2.temp/table.tail
echo "\\label{tab:TauGlobalFit_unconstrained_aleph_hcorr_BBUniv_BB}" >> maketable2.temp/table.tail
else if ( $INPUT == 2 ) then
echo "\\caption{Results for branching fractions (in \%) from constrained fit to data from non-B-Factories and including those from B-Factories. PDG modification of ALEPH correlation matrix is used as input.}" >> maketable2.temp/table.tail
echo "\\label{tab:TauGlobalFit_constrained_BBUniv_BB}" >> maketable2.temp/table.tail
else if ( $INPUT == 3 ) then
echo "\\caption{Results for branching fractions (in \%) from constrained fit to data from non-B-Factories and including those from B-Factories.}" >> maketable2.temp/table.tail
echo "\\label{tab:TauGlobalFit_constrained_aleph_hcorr_BBUniv_BB}" >> maketable2.temp/table.tail
endif
echo "\\end{center}" >> maketable2.temp/table.tail
echo "\\end{table}" >> maketable2.temp/table.tail
echo "\\end{document}" >>  maketable2.temp/table.tail

cat maketable2.head maketable2.temp/table_lep.body maketable2_ns.head maketable2.temp/table_ns.body maketable2_s.head maketable2.temp/table_s.body maketable2.temp/table.tail | tee readpdg_${LINPUT}_table2.tex
