#!/bin/bash
angle=$1
shift
for i in $*;
do cd $i/$angle;
dsolve; # works if LEOPARD executable located in $PATH -- otherwise replace this line with path to LEOPARD executable  
cp omega.dat omega.dat_latestbackup;
cd ../..;
done;
