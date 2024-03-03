#include <iostream>
#include <string>

struct CmdArgs {
    size_t natoms;
    size_t nfaces;
    size_t nruns;
    size_t nwarmup;
    size_t nisomers;
    size_t nlanczos;
    std::string device_type;
    std::string output_file;
    bool use_double_precision;;

    CmdArgs() = default;
};

void parseArguments(int argc, char** argv, CmdArgs& args);
