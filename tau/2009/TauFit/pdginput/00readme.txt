./any.com readpdg
./readpdg.exe | tee readpdg.log
../../../../combos/combos combos_average_pdginput_no_babar_belle.input | tee combos.log
../../../../combos/combos combos_average_pdginput_no_babar_belle_scaled.input | tee combos_scaled.log
./compare_combos_pdg.csh | tee compare_combos_pdg.log
