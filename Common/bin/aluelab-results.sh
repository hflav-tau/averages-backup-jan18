#!/bin/sh
# -*- ksh -*-

#
# aluelab-results.sh
#
# execute alucomb.r and print aluelab-results.r results
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

../../../Common/bin/alucomb.sh $file
if [ $? != 0 ] ; then
  exit 1
fi
fileres="${name}-results.log"
echo 1>&2 "producing file: ${fileres}"
/bin/rm -f ${fileres} ; ./aluelab-results.r | tee ${fileres}

exit
