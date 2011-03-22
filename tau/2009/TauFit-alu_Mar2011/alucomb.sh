#!/bin/sh
# -*- ksh -*-

#
# alucomb.sh
#

usage ()
{
  echo 1>&2 "usage: `basename $0` <input file>"
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
  usage
fi

name=`expr "$file" : '\(.*\)\..*'`

/bin/rm -f average.input
/bin/ln -s $file average.input

echo 1>&2 "make alucomb"
make alucomb
if [ $? != 0 ] ; then
  exit $?
fi

echo 1>&2 "produced file: ${name}.rdata"
echo 1>&2 "produced file: ${name}.log"
echo 1>&2 "produced file: ${name}.txt"

/bin/cp -p average_alucomb.rdata ${name}.rdata
/bin/cp -p average_alucomb.log ${name}.log
/bin/rm -f ${name}.txt ; ./aluelab-results.r | tee ${name}.txt

exit
