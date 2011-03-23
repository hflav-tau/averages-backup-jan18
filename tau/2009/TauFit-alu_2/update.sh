#!/bin/sh
# -*- ksh -*-

#
# update.sh
#
# update results from several different input files
#

usage ()
{
  echo 1>&2 "usage: `basename $0` [ -l ]"
  exit 1
}

log=0
while getopts hl opt; do
  case "$opt" in
    h) usage ;;
    l) log=1 ;;
    *) usage ;;
  esac
done
shift `expr $OPTIND - 1`

for input in \
 alucomb_average_constrained.input \
 alucomb_average_constrained_aleph_hcorr.input \
 alucomb_average_rescaled_constrained.input \
 alucomb_average_rescaled_constrained_aleph_hcorr.input \
 alucomb_average_unconstrained.input \
 alucomb_average_unconstrained_aleph_hcorr.input \
 alucomb_average_rescaled_unconstrained.input \
 alucomb_average_rescaled_unconstrained_aleph_hcorr.input \
  ; do

  echo "========"
  echo "processing $input"
  ../../../Common/bin/alucomb.sh $input
  name=`expr "$input" : '\(.*\)\..*'`
  /bin/rm "${name}.list"
  ../../../Common/bin/aluelab-list.r "${name}.rdata" | tee "${name}.list"
  if [ "$log" = "1" ] ; then
    /bin/cp -p "${name}.log" log/ 
    /bin/cp -p "${name}.list" log/ 
    : # /bin/cp -p "${name}.show" log/
    : # /bin/cp -p "${name}.res" log/
  fi
done

exit
