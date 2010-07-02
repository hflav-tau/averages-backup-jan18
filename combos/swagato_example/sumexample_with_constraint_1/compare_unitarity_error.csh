#!/bin/csh

setenv LOGFILE compare_unitarity_error.summary
rm -f $LOGFILE
touch $LOGFILE

foreach i ( 0 1 2 3 4 5 6 7)
 rm -f unitarity.input 
 sed -e 's/SUBSTITUTE/1.0e-'$i'/g' unitarity_template.input > unitarity.input 
 rm -f average_unitarity_1e-$i.log
 ../../combos average.input > average_unitarity_1e-$i.log
 echo "  // Unitarity Error = 1.0e-$i"  >> $LOGFILE
 grep "Tot Err:" average_unitarity_1e-$i.log | tail -2 | awk -v I=$i 'BEGIN{printf "  float val_"I"[nx] = {"} {if (NR==2) {printf $2     }else{printf $2     ", "}} END{printf "};\n"}' >> $LOGFILE
 grep "Tot Err:" average_unitarity_1e-$i.log | tail -2 | awk -v I=$i 'BEGIN{printf "  float err_"I"[nx] = {"} {if (NR==2) {printf $(NF-3)}else{printf $(NF-3)", "}} END{printf "};\n"}' >> $LOGFILE
end

echo grep "All: Sum of Multi-Quantity Average" average_unitarity_1e-*.log  | awk '{if (NR%2==0) print $0}' >> $LOGFILE
     grep "All: Sum of Multi-Quantity Average" average_unitarity_1e-*.log  | awk '{if (NR%2==0) print $0}' >> $LOGFILE

cat $LOGFILE
