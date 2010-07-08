#!/bin/sh
# -*- ksh -*-

#
# in order to report in the plots both the single channel combination
# and the multiple channel combination, both plots should be produced
# in the correct order
# use this script in the directories:
# - TauTo3Prongs
# - TauTo1Prong
#

set -x
make aluplot
dirs=`/bin/ls -1 plot-*.input | /bin/sed -e 's/plot-\(.*\)[.]input/\1/'`
for dir in $dirs ; do
  cd "../TauTo$dir/"
  make aluplot plot link update_plot
done
