#!/bin/bash

echo $SLURM_LOCALID
env  > env.$SLURM_LOCALID.txt

cat /proc/cpuinfo > cpuinfo.$SLURM_LOCALID.txt

(sycl-ls > syclls.$SLURM_LOCALID.txt) 2> syclls.$SLURM_LOCALID.err

./sycltest
