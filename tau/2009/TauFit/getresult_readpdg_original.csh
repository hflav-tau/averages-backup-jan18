#!/bin/csh -f
setenv LOGFILE $1
setenv TXTFILE `echo ${LOGFILE} | awk -F. '{print $1}' | awk -F"readpdg_" '{print $2}' | awk '{print "getresult_readpdg_average_"$1".txt"}'`
rm -f ${TXTFILE}_1 ; grep -A38 -F "Results from original fit:" $LOGFILE | tail -n +3 | awk '{printf "Gamma%-3s %8.6f %8.6f\n",$3,$7,$8}' > ${TXTFILE}_1
rm -f ${TXTFILE}_2 ; grep strange $LOGFILE | head -1 | awk -F= '{print $3}' | awk '{printf "Gamma110 %8.6f %8.6f\n",$1,$3}'  > ${TXTFILE}_2
cat ${TXTFILE}_1 ${TXTFILE}_2 | tee ${TXTFILE}
rm -f ${TXTFILE}_1 ${TXTFILE}_2
