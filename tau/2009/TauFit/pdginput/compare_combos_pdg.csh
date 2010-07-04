#!/bin/csh
grep "Tot Err:" combos.log | tail -30 >! a1
grep basefitvalue readpdg.cc | grep push_back | head -30 | awk -F"basefitvalue" '{print $2}' | sed -e 's/=/ /g' | sed -e 's/\[ibase\]//g' | sed -e 's/; basefiterror//g'  | sed -e 's/; baserescalederror //g' | sed -e 's/; basescalefactor//g ' | awk -F";" '{print $1}' >! a2
paste a1 a2 |  awk '{print $1, $2,$9,$(NF-3),$(NF-2),$(NF-1),$NF,100.*(($2/$(NF-3))-1),"%",($9/$2)*100,"%"}' | awk '{if ($8>0) {print $0,$8}else{print $0,-$8}}' | sort -k12 -g
# echo "M_GAMMA104 0.0001749 0.0000256 0.000181 0.000026 0.000026 1.01 -3.37017 % 14.6369 %" | awk '{if ($8>0) {print $0,$8}else{print $0,-$8}}' | sort -k12 -g

# echo "M_GAMMA30 0.0014215 0.0003084 0.001024 0.000392 0.000398 1.01 38.8184 % 21.6954 %"   | awk '{if ($8>0) {print $0,$8}else{print $0,-$8}}' | sort -k12 -g

# echo "M_GAMMA103 0.0008194 0.0000516 0.000810 0.000053 0.000055 1.05  1.16049 % 6.29729 %" | awk '{if ($8>0) {print $0,$8}else{print $0,-$8}}' | sort -k12 -g

# echo "M_GAMMA103 0.0008102 0.0000526 0.000810 0.000053 0.000055 1.05  0.0246914 % 6.49222 %" | awk '{if ($8>0) {print $0,$8}else{print $0,-$8}}' | sort -k12 -g
  echo "M_GAMMA103 0.0008101 0.0000526 0.000810 0.000053 0.000055 1.05  0.0123457 % 6.49303 %" | awk '{if ($8>0) {print $0,$8}else{print $0,-$8}}' | sort -k12 -g
