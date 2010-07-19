cd ../../../../Common/ ; source ./config.csh ; cd -
./any.com readpdg
./readpdg.exe | tee readpdg.log

cp -p alucomb_average_pdginput_no_babar_belle.input log/
cp -p alucomb_pdginput_measurements.input log/

cp -p combos_average_pdginput_no_babar_belle.input log/
cp -p combos_pdginput_measurements.input log/

cp -p combos_average_pdginput_no_babar_belle_scaled.input log/
cp -p combos_pdginput_measurements_scaled.input log/

../../../../combos/combos combos_average_pdginput_no_babar_belle.input | tee combos.log
../../../../combos/combos combos_average_pdginput_no_babar_belle_scaled.input | tee combos_scaled.log
./compare_combos_pdg.csh | tee compare_combos_pdg.log

cp -p combos.log log/

cp -p combos_scaled.log log/

cp -p compare_combos_pdg.log log/

cp -p compare_alucomb_combos_pdg.log log/



