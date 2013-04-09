#!/bin/bash
#SBATCH -p octuplets 
#SBATCH -o $JOBDIR/output-$I.log
#SBATCH -D $ROOT
#SBATCH -J HamiltonPaths
#SBATCH -n $NCPUS
#SBATCH --mem-per-cpu=500

for i in `seq $iFROM $iTO`; do 
    ./fullerene < $JOBDIR/$i/input.inp > $JOBDIR/$i/output.log &
done

hostname
echo "Running $iFROM to $iTO. Cancel job when all programs are done."
while true; do 
    top -b -n1 | head -n 20;
    sleep 10000;
done
echo "This should never be reached."


