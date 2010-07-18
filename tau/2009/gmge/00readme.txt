root -l -b -q phasespace.C | tee phasespace.log

../../../Common/epc -s0 phasespace_err.input | phasespace_err.log

../../../Common/epc -s0 gmge.input | tee gmge.log


