#!/bin/csh
# compile standalone : no MyLib.so needed in tmp directory
set NAME="$1"
g++ -o ${NAME}.exe ${NAME}.cc `root-config --cflags --libs` -g -Wno-deprecated -Llib -lHist -lGraf -lMatrix -lTree -lMinuit -lRIO \
                         -lMathCore -lFoam -lProof
exit

