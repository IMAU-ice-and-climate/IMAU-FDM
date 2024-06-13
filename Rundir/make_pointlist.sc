#!/bin/bash

npoints=58265

ipoint=0

echo "We make a pointlist with $npoints points."

rm -f pointlist.txt
let "epoint=$npoints+$ipoint"
while [ $ipoint -lt $epoint ]; do
  let "ipoint+=1"
  echo $ipoint >> pointlist.txt
done

exit 0
  
