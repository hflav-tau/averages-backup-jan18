#!/bin/csh
rm -f a1_1 ; grep "Tot Err:" combos.log | tail -30 | awk '{printf "%s%10.6f%10.6f\n", $1, $2, $(NF-3)}' | sed -e 's/M_GAMMA/Gamma/g' > a1_1
rm -f a1_2 ; grep "All: Sum of Multi-Quantity Average" combos.log | tail -1 | awk -F"=" '{print $2}' | awk '{printf "%s%10.6f%10.6f\n", "Gamma103",1-$1,$3}' > a1_2
rm -f a1 ; cat a1_1 a1_2 > a1

rm -f a1s_1 ; grep "Tot Err:" combos_scaled.log | tail -30 | awk '{print $1, $2, $(NF-3)}' | sed -e 's/M_GAMMA/Gamma/g' > a1s_1
rm -f a1s_2 ; grep "All: Sum of Multi-Quantity Average" combos_scaled.log | tail -1 | awk -F"=" '{print $2}' | awk '{print "Gamma103",1-$1,$3}' > a1s_2
rm -f a1s ; cat a1s_1 a1s_2 > a1s

rm -f a1ss ; paste a1 a1s | awk '{if ($NF > $3) {printf "%12s %9.6f %9.6f %9.6f %4.2f\n", $1,$2,$3,$NF,$NF/$3}else{printf "%12s %9.6f %9.6f %9.6f %4.2f\n", $1,$2,$3,$3,1}}' > a1ss

rm -f a2 ; grep basefitvalue readpdg.cc | grep push_back | head -31 | awk -F"basefitvalue" '{print $2}' | sed -e 's/=/ /g' | sed -e 's/\[ibase\]//g' | sed -e 's/; basefiterror//g'  | sed -e 's/; baserescalederror //g' | sed -e 's/; basescalefactor//g ' | awk -F";" '{print $1}' | awk '{printf "%9.6f %9.6f %9.6f %4.2f\n", $1,$2,$3,$4 }' > a2

rm -f a0; grep basefitvalue readpdg.cc | grep push_back | head -31 | awk -F"back" '{print $2,$3}' | awk -F"baseseed" '{print $1}' | sed -e 's/)/ /g' -e 's/(/ /g' | awk '{print $4}' > a0

echo "PARNUMBER  DECAYMODE COMBOS_VAL COMBOS_ERR  SCALED SFAC  PDG_VAL   PDG_ERR PDG_SCALED SFAC PDG-COMBOS_VAL PDG-COMBOS_ERR PDG-COMBOS_SCALED PDG-COMBOS_SFAC"
paste a0 a1ss a2 |  awk '{printf "%s %14.6f %14.6f %14.6f %14.2f\n", $0, ($7-$3), ($8-$4), ($9-$5), ($10-$6) }'

rm -f a1_1 a1_2 a1 
rm -f a1s_1 a1s_2 a1s
rm -f a1ss
rm -f a2 a0

