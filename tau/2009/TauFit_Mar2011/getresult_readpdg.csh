#!/bin/csh -f
foreach logfile ( readpdg_constrained.log readpdg_constrained_aleph_hcorr.log readpdg_unconstrained.log readpdg_unconstrained_aleph_hcorr.log )
 ./getresult_readpdg_original.csh $logfile
 ./getresult_readpdg_rescaled.csh $logfile
end
