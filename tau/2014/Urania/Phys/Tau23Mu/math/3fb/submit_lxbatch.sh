cat > ./combinelimits_2012_tau23mu_V4_dev.job <<EOF
cd /afs/cern.ch/user/s/swaban/scratch1/hfagtau-averages_2018/tau/2014/Urania
source ./setup.sh
cd /afs/cern.ch/user/s/swaban/scratch1/hfagtau-averages_2018/tau/2014/Urania/Phys/Tau23Mu/math/3fb
python ./combinelimits_2012_tau23mu_V4_dev.py
EOF
chmod +x ./combinelimits_2012_tau23mu_V4_dev.job
rm -f combinelimits_2012_tau23mu_V4_dev.log
/afs/cern.ch/user/s/swaban/public/bsub_lxplus 1nd combinelimits_2012_tau23mu_V4_dev.log 1 ./combinelimits_2012_tau23mu_V4_dev.job

cat > run_batch_tau23mu_0_0.job <<EOF
cd /afs/cern.ch/user/s/swaban/scratch1/hfagtau-averages_2018/tau/2014/Urania
source ./setup.sh
cd /afs/cern.ch/user/s/swaban/scratch1/hfagtau-averages_2018/tau/2014/Urania/Phys/Tau23Mu/math/3fb
python ./run_batch_tau23mu.py 0 0
EOF
chmod +x run_batch_tau23mu_0_0.job
rm -f run_batch_tau23mu_0_0.log
/afs/cern.ch/user/s/swaban/public/bsub_lxplus 1nd run_batch_tau23mu_0_0.log 1 ./run_batch_tau23mu_0_0.job
