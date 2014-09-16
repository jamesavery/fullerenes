#!/bin/bash
# @ output = $JOBDIR/output-$I.log
# @ error =  $JOBDIR/error-$I.log
# @ wall_clock_limit = 72:00:00
# @ class = large
# @ resources = ConsumableCpus(8) ConsumableMemory(2000mb) 
# @ queue

cd $ROOT

for i in `seq $iFROM $iTO`; do 
    (./$EXECNAME `cat $JOBDIR/$i/input.inp` > $JOBDIR/$i/output.log) 2> $JOBDIR/$i/output.err &
done

hostname
echo "Running $iFROM to $iTO."
wait
echo "Range $iFROM to $iTO is done."
#  Cancel job when all programs are done."
# while true; do 
#     top -b -n1 | head -n 20;
#     sleep 10000;
# done
# echo "This should never be reached."


