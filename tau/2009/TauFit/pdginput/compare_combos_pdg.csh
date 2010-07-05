#!/bin/csh
rm -f a1_1 ; grep "Tot Err:" combos.log | tail -30 | awk '{print $1, $2, $(NF-3)}' > a1_1
rm -f a1_2 ; grep "All: Sum of Multi-Quantity Average" combos.log | tail -1 | awk -F"=" '{print $2}' | awk '{print "M_Gamma103",1-$1,$3}' > a1_2
rm -f a1 ; cat a1_1 a1_2 > a1
rm -f a2 ; grep basefitvalue readpdg.cc | grep push_back | head -31 | awk -F"basefitvalue" '{print $2}' | sed -e 's/=/ /g' | sed -e 's/\[ibase\]//g' | sed -e 's/; basefiterror//g'  | sed -e 's/; baserescalederror //g' | sed -e 's/; basescalefactor//g ' | awk -F";" '{print $1}' > a2

echo "DECAYMODE COMBOS_VAL COMBOS_ERR  PDG_VAL    ERROR   RESCALED  SFAC 1-COMBOS/PDG_VAL 1-COMBOS/PDG_ERR    PDG_VAL/ERR |COMBOS/PDG_VAL|"
paste a1 a2 |  awk '{printf "%s %12.6f %s %12.6f %s %12.6f %s \n", $0, 100.*(1-($2/$4)), "%", 100.*(1-($3/$5)), "%", 100.*($5/$4), "%"}' | awk '{if ($8>0) {print $0,$8}else{print $0,-$8}}' | sort -k14 -g
