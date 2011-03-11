 tail -n +73 ../TauFit-alu/aluelab-show_new.log | awk '{if (NF==1) {s=$1}else{ print s," ", $0}}' | grep -e BaBar -e Belle | \
           sed -e 's/BaBar.Gamma3by5.submitted/135/g' \
               -e 's/BaBar.Gamma9by5.submitted/136/g' \
               -e 's/BaBar.Gamma10by5.submitted/137/g' \
               -e 's/BaBar.Gamma16.PUBLISHED/138/g' \
               -e 's/BaBar.Gamma35.preliminary/139/g' \
               -e 's/BaBar.Gamma40.preliminary/140/g' \
               -e 's/BaBar.Gamma60.published/141/g' \
               -e 's/BaBar.Gamma85.published/142/g' \
               -e 's/BaBar.Gamma93.published/143/g' \
               -e 's/BaBar.Gamma103.published/144/g' \
               -e 's/BaBar.Gamma96.published/145/g' \
               -e 's/BaBar.Gamma136.published/146/g' \
               -e 's/Belle.Gamma13.published/147/g' \
               -e 's/Belle.Gamma35.published/148/g' \
               -e 's/Belle.Gamma60.submitted/149/g' \
               -e 's/Belle.Gamma85.submitted/150/g' \
               -e 's/Belle.Gamma93.submitted/151/g' \
               -e 's/Belle.Gamma126.published/152/g' \
               -e 's/Belle.Gamma128.Published/153/g' \
               -e 's/Belle.Gamma130.Published/154/g' \
               -e 's/Belle.Gamma132.Published/155/g' \
               -e 's/Belle.Gamma96.submitted/156/g' | \
               awk '{if ($NF>0) print $0}' | awk '{if ($1<$2) print $0}' | awk '{printf "%% %4d %4d %14.10g\n",$1,$2,$3*100}'
