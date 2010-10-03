cd ../../../../Common/ ; source ./config.csh ; cd -
./any.com readpdg

source run_readpdg

source run_combos_constrained
source run_combos_constrained_aleph_hcorr
source run_combos_unconstrained
source run_combos_unconstrained_aleph_hcorr

cp -p *.log log

cp -p combos*.input log/input_dir
cp -p alucomb*.input log/input_dir


