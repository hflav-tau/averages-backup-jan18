#!/bin/csh

#setenv LOGFILE compare_unitarity_error.log
#rm -f $LOGFILE
#touch $LOGFILE
#
#foreach i ( 0 1 2 3 4 5 6 7)
# rm -f unitarity.input 
# sed -e 's/SUBSTITUTE/1.0e-'$i'/g' unitarity_template.input > unitarity.input 
# make combos
# cp -p average.log log/average_unitarity_1e-$i.log
# echo "  // Unitarity Error = 1.0e-$i"  >> $LOGFILE
# grep "Tot Err:" log/average_unitarity_1e-$i.log | tail -34 | awk -v I=$i 'BEGIN{printf "  float val_"I"\[nx\] = \{"}{printf $2     ", "}END{printf "\};\n"}' >> $LOGFILE
# grep "Tot Err:" log/average_unitarity_1e-$i.log | tail -34 | awk -v I=$i 'BEGIN{printf "  float err_"I"\[nx\] = \{"}{printf $(NF-3)", "}END{printf "\};\n"}' >> $LOGFILE
#end

# rm -f unitarity.input ; ln -s unitarity_default.input unitarity.input

grep NDOFNEW log/average_unitarity_1e-*.log | awk '{if (NR%4==3||NR%4==0) print $0}' | grep -v SCALE
grep NDOFNEW log/average_unitarity_1e-*.log | awk '{if (NR%4==3||NR%4==0) print $0}' | grep -e SCALE
grep "All: Sum of Multi-Quantity Average" log/*.log | awk '{if (NR%2==0) print $0}'
 
# grep "Tot Err:" log/average_unitarity_1e-4.log | tail -34 

