 tail -n +73 ../TauFit-alu/aluelab-show_new.log | awk '{if (NF==1) {s=$1}else{ print s," ", $0}}' | grep -e BaBar -e Belle | \
           sed -e 's/BaBar.Gamma3by5.submitted/132/g' \
               -e 's/BaBar.Gamma9by5.submitted/133/g' \
               -e 's/BaBar.Gamma10by5.submitted/134/g' \
               -e 's/BaBar.Gamma16.PUBLISHED/135/g' \
               -e 's/BaBar.Gamma35.preliminary/136/g' \
               -e 's/BaBar.Gamma40.preliminary/137/g' \
               -e 's/BaBar.Gamma60.published/138/g' \
               -e 's/BaBar.Gamma85.published/139/g' \
               -e 's/BaBar.Gamma93.published/140/g' \
               -e 's/BaBar.Gamma103.published/141/g' \
               -e 's/BaBar.Gamma96.published/142/g' \
               -e 's/BaBar.Gamma136.published/143/g' \
               -e 's/Belle.Gamma13.published/144/g' \
               -e 's/Belle.Gamma35.published/145/g' \
               -e 's/Belle.Gamma60.submitted/146/g' \
               -e 's/Belle.Gamma85.submitted/147/g' \
               -e 's/Belle.Gamma93.submitted/148/g' \
               -e 's/Belle.Gamma126.published/149/g' \
               -e 's/Belle.Gamma128.Published/150/g' \
               -e 's/Belle.Gamma130.Published/151/g' \
               -e 's/Belle.Gamma132.Published/152/g' \
               -e 's/Belle.Gamma96.submitted/153/g' | \
               awk '{if ($NF>0) print $0}' | awk '{if ($1<$2) print $0}' | awk '{printf "%% %4d %4d %14.10g\n",$1,$2,$3*100}'
