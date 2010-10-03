#!/bin/csh
foreach txtfile ( average_constrained.txt \
                  average_constrained_aleph_hcorr.txt \
                  average_rescaled_constrained.txt \
                  average_rescaled_constrained_aleph_hcorr.txt \
                  average_unconstrained.txt \
                  average_unconstrained_aleph_hcorr.txt \
                  average_rescaled_unconstrained.txt \
                  average_rescaled_unconstrained_aleph_hcorr.txt )
setenv TXTFILE_1 `echo $txtfile | awk '{print "getresult_readpdg_"$1}'`
setenv TXTFILE_2 `echo $txtfile | awk '{print "getresult_alucomb_"$1}'`
setenv TXTFILE_3 `echo $txtfile | awk '{print "compare_readpdg_alucomb_"$1}'`

paste $TXTFILE_1 $TXTFILE_2  | awk '{printf "%s  %10.6e  %10.6e\n", $0, $2-$5, $3-$6}' | tee $TXTFILE_3

end
