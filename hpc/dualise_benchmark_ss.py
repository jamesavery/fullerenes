import os

##### Benchmarking parameters / constants #####
# TODO change these to correspond to the benchmark!

# Paths
base_dir = os.path.dirname(os.path.realpath(__file__))
script_dir = f'{base_dir}'
run_dir = f'{base_dir}/../build'
gen_dir = f'{run_dir}/gen'
output_dir = f'{run_dir}/output'

# Program parameters
batch_size = 1000
n_isomer = 200

# HPC parameters
cpus_per_task = 4
n_gpus_to_test = [2**i for i in range(4)]
n_nodes_to_test = [] # [2**i for i in range(1,8)]

workers_per_task=1              # Increase for more buckyherd parallelism, but more IPC overhead
max_chunks_per_worker=1         # Increase for better load balancing, but more buckygen overhead

if(len(n_nodes_to_test)>1):
    N_chunks = n_gpus_to_test[-1]*n_nodes_to_test[-1]*workers_per_task*max_chunks_per_worker
else:
    N_chunks = n_gpus_to_test[-1]*workers_per_task*max_chunks_per_worker

#n_gpus_to_test = [2**i for i in range(1)]
#n_nodes_to_test = [2]#[2**i for i in range(1)]

###################################

def run_test(n_nodes, n_gpus, n_atoms=100):
    # Run the test
    jobname = f"test_{n_nodes}_{n_gpus}"
    sed_lut = {
        '@@JOBNAME@@': jobname,
        '@@NNODES@@': str(n_nodes),
        '@@NTASKSPERNODE@@': str(n_gpus),
        '@@CPUSPERTASK@@': str(cpus_per_task),
        '@@NTASKS@@': str(n_nodes * n_gpus),
        '@@SCRIPTDIR@@': script_dir,
        '@@RUNDIR@@': run_dir,
        '@@OUTPUTDIR@@': output_dir
    }
    job_path = f'{gen_dir}/{jobname}.sh'
    with open(f'{script_dir}/multinode-bench.sh', 'r') as fin:
        with open(job_path, 'w') as fout:
            for line in fin:
                for key, value in sed_lut.items():
                    line = line.replace(key, value)
                fout.write(line)
    os.system(f'chmod +x {job_path}')
    os.system(f'sbatch {job_path} ./benchmarks/sycl/full_pipeline_strong_scaling {n_atoms} {N_chunks} {workers_per_task}')
    #os.system(f"sbatch job.sh")

if __name__ == '__main__':
    # Create directories
    os.system(f'mkdir -p {gen_dir}')
    os.system(f'mkdir -p {output_dir}')
    #get N from command line
    import sys
    n_atoms = 100
    if len(sys.argv) > 1:
        n_atoms = int(sys.argv[1])

    # Run single node scaling
    for n_gpus in n_gpus_to_test:
        n_nodes = 1
        run_test(n_nodes, n_gpus, n_atoms=n_atoms)

    # Run multi node scaling
    for n_nodes in n_nodes_to_test:
        n_gpus = n_gpus_to_test[-1]
        run_test(n_nodes, n_gpus, n_atoms=n_atoms)
