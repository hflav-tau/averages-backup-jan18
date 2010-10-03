#!/bin/csh -f
foreach txtfile ( average_constrained.txt \
                  average_constrained_aleph_hcorr.txt \
                  average_rescaled_constrained.txt \
                  average_rescaled_constrained_aleph_hcorr.txt \
                  average_unconstrained.txt \
                  average_unconstrained_aleph_hcorr.txt \
                  average_rescaled_unconstrained.txt \
                  average_rescaled_unconstrained_aleph_hcorr.txt )
setenv infile getresult_readpdg_$txtfile
setenv listfile `echo $txtfile | awk -F. '{print "alucomb_"$1".list"}'`
setenv outfile getresult_alucomb_$txtfile

less $infile | awk '{printf "grep -F \"%-8s \" ../TauFit-alu_2/%s\n",$1,LIST}' LIST=${listfile} | sh | awk '{printf "%-8s %8.6f %8.6f\n",$1,$2,$3}' | tee $outfile

end
