root -l -b -q predicted.C | tee predicted.log

../../../Common/epc -s0 tautopinu.input | tee tautopinu.log

../../../Common/epc -s0 gtaugm_pi.input | tee gtaugm_pi.log

