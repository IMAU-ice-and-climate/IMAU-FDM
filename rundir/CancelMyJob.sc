#!/bin/bash
echo "Cancel job FGRN055_cores. If that is not what you want, you have five seconds left."
sleep 5

echo "Start canceling"
# post to all thread request files "fatal". To be sure that the distributor picks it up
for it in {0..127}; do
  requestfile=/ec/res4/scratch/nld4814/cores/requests/`printf "%.5d\n" $it`
  echo "fatal" > $requestfile
done  
echo "Done"
exit 0
