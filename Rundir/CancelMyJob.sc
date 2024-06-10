#!/bin/bash
echo "Cancel job FGRN055_tskrun_min1. If that is not what you want, you have five seconds left."
sleep 5

echo "Start canceling"
# post to all thread request files "fatal". To be sure that the distributor picks it up
for it in {0..32}; do
  requestfile=/scratch/ms/nl/rumb/FDMtests/FGRN055/FGRN055_tskrun_min1/requests/`printf "%.5d\n" $it`
  echo "fatal" > $requestfile
done  
echo "Done"
exit 0
