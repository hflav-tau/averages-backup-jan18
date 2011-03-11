#!/bin/csh -f
foreach logfile ( combos_average_constrained.log \
                  combos_average_constrained_aleph_hcorr.log \
                  combos_average_rescaled_constrained.log \
                  combos_average_rescaled_constrained_aleph_hcorr.log )
getresult_combos_constrained.csh $logfile
end

foreach logfile ( combos_average_unconstrained.log \
                  combos_average_unconstrained_aleph_hcorr.log \
                  combos_average_rescaled_unconstrained.log \
                  combos_average_rescaled_unconstrained_aleph_hcorr.log )
getresult_combos_unconstrained.csh $logfile
end
