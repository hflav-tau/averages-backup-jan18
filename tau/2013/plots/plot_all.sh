#!/bin/bash

for i in *.input
do
    root -b -q "aluAverPlot.cc+(\"$i\")"
done
