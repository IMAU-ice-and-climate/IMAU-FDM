#!/bin/bash


input_path="$SCRATCH/daily/output"
output_path="$SCRATCH/daily/output/processed_output"


for file in "$input_path"/FGRN055_era055_2Ddetail_*.nc; do
    filename=$(basename "$file")
    new_filename="$filename"
    #ncrename -d ind_t,time "$file"
    cdo timmean "$file" "$output_path/$new_filename"
    echo "Processed: $file -> $output_path/$new_filename"
done

