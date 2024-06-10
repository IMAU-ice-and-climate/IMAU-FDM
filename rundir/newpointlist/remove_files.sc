#!/bin/bash

outdir="$SCRATCH/data/output/era055/newrun/"
requestdir="/scratch/ms/nl/rumb/FDMtests/FGRN055/FGRN055_imp_run/requests/"
readydir="/scratch/ms/nl/rumb/FDMtests/FGRN055/FGRN055_imp_run/readydir/"
workexe="/scratch/ms/nl/rumb/FDMtests/FGRN055/FGRN055_imp_run/LocalCode/"

name_part1="ECMWF_FGRN055_imp_run_*_" 

echo "We remove all partly finished results. You have 5 seconds to cancel this script."

sleep 5

echo "Removing all unnecessary files."

while IFS= read -r line; do
    namefile="$outdir$name_part1$line.nc"
    #printf '%s\n' "$namefile"
    rm -f $namefile
done < removelist.txt

rm -f $requestdir/*
rm -f $readydir/*
rm -f $workexe/*

echo "All files removed."

