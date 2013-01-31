#!/bin/bash
# @ output = $JOBDIR/output-$I.log
# @ error =  $JOBDIR/error-$I.log
# @ wall_clock_limit = 48:00:00
# @ class = large
# @ resources = ConsumableCpus(8) ConsumableMemory(2000mb) 
# @ queue

cd $ROOT

for i in `seq $iFROM $iTO`; do 
    ./fullerene < $JOBDIR/$i/input.inp > $JOBDIR/$i/output.log &
done

echo "Running $iFROM to $iTO. Cancel job when all programs are done."
while true; do 
    top -b -n1 | head -n 20;
    sleep 10000;
done
echo "This should never be reached."


