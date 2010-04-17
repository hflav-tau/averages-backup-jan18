root -l -b -q predicted.C | tee predicted.log

../../../Common/epc -s0 tautoknu.input | tee tautoknu.log

../../../Common/epc -s0 gtaugm_k.input | tee gtaugm_k.log

