#!/bin/sh
# -*- ksh -*-

dir=`dirname "$0"`
if [ "$dir" = "" ] ; then
  dir="."
fi
logfile="$1"
if [ "$logfile" = "" ] ; then
  logfile="average_alucomb.log"
fi
if grep -q "aluelab-taufit.r:" "$logfile" ; then
  echo 1>&2 "aluelab-taufit.r info already present in $logfile"
else
  "${dir}/aluelab-taufit.r" "$logfile" >> "$logfile"
  echo 1>&2 "added aluelab-taufit.r to $logfile"
fi
