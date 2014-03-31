#!/bin/bash
#SBATCH -p low 
#SBATCH -o $JOBDIR/output-$I.log
#SBATCH -D $ROOT
#SBATCH -J $EXECNAME-$I
#SBATCH --nodes=1
#SBATCH --ntasks=$NCPUS
#SBATCH --cpus-per-task=1
#SBATCH --mem=5000

for i in `seq $iFROM $iTO`; do 
    (./$EXECNAME `cat $JOBDIR/$i/input.inp` > $JOBDIR/$i/output.log) 2> $JOBDIR/$i/output.err &
done

hostname
echo "Running $iFROM to $iTO. Cancel job when all programs are done."
while true; do 
    top -b -n1 | head -n 20;
    sleep 10000;
done
echo "This should never be reached."


