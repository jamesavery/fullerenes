#!/bin/bash -l
#SBATCH --job-name=@@JOBNAME@@    # Job name
#SBATCH --output=@@OUTPUTDIR@@/@@JOBNAME@@.o%j  # Name of stdout output file
#SBATCH --error=@@OUTPUTDIR@@/@@JOBNAME@@.e%j   # Name of stderr error file
#SBATCH --partition=standard-g  # Partition (queue) name
#SBATCH --nodes=@@NNODES@@               # Total number of nodes
#SBATCH --ntasks-per-node=@@NTASKSPERNODE@@    # 8 MPI ranks per node, 16 total (2x8)
#SBATCH --gpus-per-node=@@NTASKSPERNODE@@       # Allocate one gpu per MPI rank
#SBATCH --cpus-per-task=@@CPUSPERTASK@@
#SBATCH --time=0-00:30:00       # Run time (d-hh:mm:ss)
#SBATCH --account=project_465000815  # Project for billing
#SBATCH --export=ALL

export SCRIPTDIR=@@SCRIPTDIR@@ # Directory of this script
export JOBNAME=%j # Name of this job
export RUNDIR=@@RUNDIR@@ # CWD of the job

# TODO LUMI does not like this, even though it's their own example
# https://docs.lumi-supercomputer.eu/runjobs/scheduled-jobs/lumig-job/
CPU_BIND="mask_cpu:7e000000000000,7e00000000000000"
CPU_BIND="${CPU_BIND},7e0000,7e000000"
CPU_BIND="${CPU_BIND},7e,7e00"
CPU_BIND="${CPU_BIND},7e00000000,7e0000000000"

if [ "$#" -eq 0 ]; then
    echo "Syntax: $0 <executable>"
    exit 1
fi

export N_TASKS=@@NTASKS@@ # Total number of tasks

#srun --cpu-bind=${CPU_BIND} ${SCRIPTDIR}/select_gpu.sh ${SCRIPTDIR}/run-bench.sh $@
srun ${SCRIPTDIR}/select_gpu.sh ${SCRIPTDIR}/run-bench.sh $@
