#include <iostream>
#include <string>

struct CmdArgs {
    size_t print_vector;
    size_t device_id;
    size_t natoms;
    size_t nfaces;
    size_t nruns;
    size_t nwarmup;
    size_t nisomers;
    size_t nlanczos;
    std::string device_type;
    std::string output_file;
    bool use_double_precision;;

    CmdArgs() : print_vector(0), device_id(0), natoms(20), nfaces(12), nruns(1), nwarmup(1), nisomers(1), nlanczos(50), device_type("gpu"), output_file("output.txt"), use_double_precision(false) {};
};

void parseArguments(int argc, char** argv, CmdArgs& args);
