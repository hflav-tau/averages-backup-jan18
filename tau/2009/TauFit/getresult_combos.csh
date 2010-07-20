#!/bin/csh -f
rm -f a1_1 ; grep "Tot Err:" average.log | tail -35 | awk '{print $1, $2, $(NF-3)}' | sed -e 's/M_GAMMA/Gamma/g' > a1_1
rm -f a1_2 ; grep "All: Sum of Multi-Quantity Average" average.log | tail -1 | awk -F"=" '{print $2}' | awk '{print "Gamma103",1-$1,$3}' > a1_2
rm -f a1_3 ; grep "P1: Sum of Multi-Quantity Average" average.log | tail -1 |  awk -F"=" '{print $2}' | awk '{print "Gamma110",$1,$3}' > a1_3
rm -f a1 ; cat a1_1 a1_2 a1_3 | tee a1
