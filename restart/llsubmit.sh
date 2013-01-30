#!/bin/sh
# @ output = $WORKDIR/output.log
# @ error =  $WORDKIR/error.log
# @ wall_clock_limit = 120:00:00
# @ class = large
# @ resources = ConsumableCpus(1) ConsumableMemory(500mb) 
# @ queue

cd $WORKDIR
../../../fullerene < input*.inp

