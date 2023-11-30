import os

##### Benchmarking parameters #####

# Program parameters
batch_size = 1000
n_isomer = 200

# HPC parameters
cpus_per_task = 4
n_gpus_to_test = [2**i for i in range(4)]
n_nodes_to_test = [2**i for i in range(11)]

###################################

def run_test(n_nodes, n_gpus):
    # Run the test
    jobname = f"test_{n_nodes}_{n_gpus}"
    sed_lut = {
        '@@JOBNAME@@': jobname,
        '@@NNODES@@': str(n_nodes),
        '@@NTASKSPERNODE@@': str(n_gpus),
        '@@CPUSPERTASK@@': str(cpus_per_task),
        '@@NTASKS@@': str(n_nodes * n_gpus)
    }
    if not os.path.exists('gen'):
        os.makedirs('gen')
    with open('multinode-bench.sh', 'r') as fin:
        with open(f'gen/{jobname}.sh', 'w') as fout:
            for line in fin:
                for key, value in sed_lut.items():
                    line = line.replace(key, value)
                fout.write(line)
    #os.system(f"sbatch job.sh")

if __name__ == '__main__':
    for n_gpus in n_gpus_to_test:
        n_nodes = 1
        run_test(n_nodes, n_gpus)
    for n_nodes in n_nodes_to_test:
        n_gpus = n_gpus_to_test[-1]
        run_test(n_nodes, n_gpus)