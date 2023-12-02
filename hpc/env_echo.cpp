#include <iostream>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    std::cout << argv[0] << " says hello!" << std::endl;
    std::cout << "argc = " << argc << std::endl;
    for (int i = 0; i < argc; ++i) {
        std::cout << "argv[" << i << "] = " << argv[i] << std::endl;
    }

    auto rocm_ret = system("rocm-smi");
    std::cout << "rocm-smi returned " << rocm_ret << std::endl;
    auto nproc_ret = system("nproc --all");
    std::cout << "nproc --all returned " << nproc_ret << std::endl;

    std::vector<std::string> env_vars = {
        "SCRIPTDIR",
        "RUNDIR",
        "N_TASKS",
        "MY_TASK_ID",
        "ROCR_VISIBLE_DEVICES",
        "SLURM_PROCID",
        "SLURM_JOB_ID",
    };
    for (auto& env_var : env_vars) {
        auto env_val = getenv(env_var.c_str());
        if (env_val) {
            std::cout << env_var << " = " << env_val << std::endl;
        } else {
            std::cout << env_var << " is not set" << std::endl;
        }
    }
}