#!/bin/bash

# Buildenv
export PROJECT_SCRATCH=/scratch/project_465000815
export PROJECT_APPL=/projappl/project_465000815
module load LUMI/23.03 partition/G
module load rocm
source $PROJECT_SCRATCH/intel/oneapi/setvars.sh

# Jobenv - didn't work from multinode-bench.sh
export MY_TASK_ID=$SLURM_PROCID # 0-indexed of which task out of N_TASKS this is

# cd into the CWD of the job
cd $RUNDIR

# Run the job
$@
