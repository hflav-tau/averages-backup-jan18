#!/bin/csh -f
setenv LOGFILE $1
setenv TXTFILE `echo ${LOGFILE} | awk -F. '{print $1}' | awk -F"readpdg_" '{print $2}' | awk '{print "getresult_readpdg_average_rescaled_"$1".txt"}'`
rm -f ${TXTFILE}_1 ; grep -A38 -F "Comparison of Results from fit with [errors rescaled in a Ad-Hoc style for kkk only] w.r.t original fit:" $LOGFILE | tail -n +3 | awk '{printf "Gamma%-3s %8.6f %8.6f\n",$3,$8,$9}' > ${TXTFILE}_1
rm -f ${TXTFILE}_2 ; grep strange $LOGFILE | head -4 | tail -1 | awk -F= '{print $3}' | sed -e 's/)//' -e 's/(//' | awk '{printf "Gamma110 %8.6f %8.6f\n",$1/100,$3/100}'  > ${TXTFILE}_2
cat ${TXTFILE}_1 ${TXTFILE}_2 | tee ${TXTFILE}
rm -f ${TXTFILE}_1 ${TXTFILE}_2
