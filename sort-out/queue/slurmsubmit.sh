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
echo "Running $iFROM to $iTO."
wait
echo "Range $iFROM to $iTO is done."


