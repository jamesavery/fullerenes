#!/bin/bash -l
#SBATCH --job-name=benchmark    # Job name
#SBATCH --output=benchmark.o%j  # Name of stdout output file
#SBATCH --error=benchmark.e%j   # Name of stderr error file
#SBATCH --partition=standard-g  # Partition (queue) name
#SBATCH --nodes=2               # Total number of nodes 
#SBATCH --ntasks-per-node=8     # 8 MPI ranks per node, 16 total (2x8)
#SBATCH --gpus-per-node=8       # Allocate one gpu per MPI rank
#SBATCH --cpus-per-task=4 
#SBATCH --time=0-00:30:00       # Run time (d-hh:mm:ss)
#SBATCH --account=project_465000815  # Project for billing
#SBATCH --mail-user=avery@ece.au.dk
#SBATCH --export=ALL

SCRIPTDIR=${HOME}/fullerenes/hpc
RUNDIR=${HOME}/fullerenes/build

CPU_BIND="mask_cpu:7e000000000000,7e00000000000000"
CPU_BIND="${CPU_BIND},7e0000,7e000000"
CPU_BIND="${CPU_BIND},7e,7e00"
CPU_BIND="${CPU_BIND},7e00000000,7e0000000000"

if [ "$#" -eq 0 ]; then
    echo "Syntax: $0 <executable>"
    exit 1
fi


srun --cpu-bind=${CPU_BIND} ${SCRIPTDIR}/select_gpu.sh ${SCRIPTDIR}/run-bench.sh ${RUNDIR}/test/sycl-tests/sycl-pipeline-multinode $SLURM_LOCALID

