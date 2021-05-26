#!/bin/bash

radius=$1
shift
for i in $*;
do mkdir -p $radius/$i/distribution; 
cp input.dat $radius/$i/;
cd $radius/$i;
python ../../print_strahl_angle.py $radius $i;
cd ../..;
done;

