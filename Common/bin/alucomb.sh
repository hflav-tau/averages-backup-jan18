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

while getopts hl opt; do
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

if [ ! -r $file -o -z $file ] ; then
  echo "file '$file' must exists and be non-empty"
  usage
fi 

name=`expr "$file" : '\(.*\)\..*'`

/bin/rm -f average.input average_alucomb.rdata average_alucomb.log
/bin/ln -s $file average.input

echo 1>&2 "make alucomb"
make alucomb
if [ $? != 0 ] ; then
  exit $?
fi

echo 1>&2 "produced file: ${name}.rdata"
echo 1>&2 "produced file: ${name}.log"

/bin/cp -p average_alucomb.rdata ${name}.rdata
/bin/cp -p average_alucomb.log ${name}.log

exit
