#!/bin/sh
# -*- ksh -*-

#
# aluelab-results2.sh
#
# execute alucomb2.r and print aluelab-results.r results
# obsolete, alucomb2.r now uses arguments
#

usage ()
{
  echo 1>&2 "usage: `basename $0`"
  exit 1
}

while getopts h opt; do
  case "$opt" in
    h) usage ;;
    *) usage ;;
  esac
done
shift `expr $OPTIND - 1`

file="$1"
if [ "$file" = "" ] ; then
  file="average.input"
fi

name=`expr "$file" : '\(.*\)\..*'`

../../../Common/bin/alucomb2.sh $file
if [ $? != 0 ] ; then
  exit 1
fi
fileres="${name}-results.log"
echo 1>&2 "producing file: ${fileres}"
/bin/rm -f ${fileres} ; ./aluelab-results.r "${name}.rdata" | tee ${fileres}

exit
