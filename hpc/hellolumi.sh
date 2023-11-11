#!/bin/bash -l
#SBATCH --job-name=hellolumi   # Job name
#SBATCH --output=hellolumi.o%j # Name of stdout output file
#SBATCH --error=hellolumi.e%j  # Name of stderr error file
#SBATCH --partition=small-g  # Partition (queue) name
#SBATCH --nodes=2               # Total number of nodes 
#SBATCH --ntasks-per-node=8     # 8 MPI ranks per node, 16 total (2x8)
#SBATCH --gpus-per-node=8       # Allocate one gpu per MPI rank
#SBATCH --time=0-00:01:00       # Run time (d-hh:mm:ss)
#SBATCH --mail-type=all         # Send email at begin and end of job
#SBATCH --account=project_465000815  # Project for billing
#SBATCH --mail-user=avery@ece.au.dk

WD=${HOME}/fullerenes/hpc/

CPU_BIND="map_cpu:49,57,17,25,1,9,33,41"

export MPICH_GPU_SUPPORT_ENABLED=1

srun --cpu-bind=${CPU_BIND} ${WD}/select_gpu.sh ${WD}/hello.sh

