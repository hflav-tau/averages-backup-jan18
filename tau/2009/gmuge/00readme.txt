root -l -b -q phasespace.C++ | tee phasespace.log

# historic, not used anymore [ treats error due to mtau in f(me/mtau) and f(mu/mtau) to be independent]
# ../../../Common/epc -s0 phasespace_err.input | tee phasespace_err.log

../../../Common/epc -s0 gmuge.input | tee gmuge.log


