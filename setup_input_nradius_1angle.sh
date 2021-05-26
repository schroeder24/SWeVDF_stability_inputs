#!/bin/bash

angle=$1
shift
for i in $*;
do mkdir -p $i/$angle/distribution; 
cp input.dat $i/$angle/;
cd $i/$angle;
python ../../print_strahl_angle.py $i $angle;
cd ../..;
done;

