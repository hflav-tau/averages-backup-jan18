#!/bin/sh
# -*- ksh -*-

cd ..
plots=`find . -maxdepth 2 -name 'plot*.eps'`
pdfs=`echo $plots | sed 's/[.]eps/.pdf/g'`
for plot in $plots; do
  echo epstopdf "$plot"
  epstopdf "$plot"
done
echo pdftk $pdfs cat output /tmp/plots.pdf
pdftk $pdfs cat output /tmp/plots.pdf
echo "producing file /tmp/plots-3x4.pdf"
pdfnup --paper a4 --nup 3x4 --scale 0.95 /tmp/plots.pdf >/dev/null 2>&1
if [ $? != 0 ] ; then
  echo "bad return code from pdftk: $?"
else
  echo "done"
fi
