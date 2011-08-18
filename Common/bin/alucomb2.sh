#!/bin/sh
# -*- ksh -*-

#
# alucomb2.sh
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

/bin/rm -f average.input average.rdata average.log

if [ ! -r $file -o -z $file ] ; then
  echo "file '$file' must exists and be non-empty"
  usage
fi 

name=`expr "$file" : '\(.*\)\..*'`

/bin/ln -s $file average.input

echo 1>&2 "make alucomb2"
make alucomb2
if [ $? != 0 ] ; then
  exit $?
fi

/bin/cp -p average.rdata ${name}.rdata
/bin/cp -p average.log ${name}.log

echo 1>&2 "produced file: ${name}.rdata"
echo 1>&2 "produced file: ${name}.log"

exit
