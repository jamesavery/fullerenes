#include <fullerenes/argparser.hh>
#include <iostream>
#include <string>
#include <algorithm>
#include <vector>

void parseArguments(int argc, char** argv, CmdArgs& args) {
    std::vector<std::string> argsVec(argv, argv + argc);
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "-natoms" || arg == "-n") args.natoms = (i + 1 < argc) ? std::stoul(argv[++i]) : args.natoms;
        else if (arg == "-nruns" || arg == "-nr") args.nruns = (i + 1 < argc) ? std::stoul(argv[++i]) : args.nruns;
        else if (arg == "-nwarmup" || arg == "-nw") args.nwarmup = (i + 1 < argc) ? std::stoul(argv[++i]) : args.nwarmup;
        else if (arg == "-nisomers" || arg == "-ni") args.nisomers = (i + 1 < argc) ? std::stoul(argv[++i]) : args.nisomers;
        else if (arg == "-nlanczos" || arg == "-nl") args.nlanczos = (i + 1 < argc) ? std::stoul(argv[++i]) : args.nlanczos;
        else if (arg == "-device_type" || arg == "-dt") args.device_type = (i + 1 < argc) ? argv[++i] : args.device_type;
        else if (arg == "-device_id" || arg == "-di") args.device_id = (i + 1 < argc) ? std::stoul(argv[++i]) : args.device_id;
        else if (arg == "-output_file" || arg == "-o") args.output_file = (i + 1 < argc) ? argv[++i] : args.output_file;
        else if (arg == "-print_vector" || arg == "-pv") args.print_vector = (i + 1 < argc) ? std::stoul(argv[++i]) : args.print_vector;
        else if (arg == "-use_double_precision" || arg == "-udp") args.use_double_precision = true;
        else if (arg == "-h") {
            std::cout << "Usage: " << argv[0] << " [-natoms <size_t>] [-nruns <size_t>] [-nwarmup <size_t>] [-nisomers <size_t>] [-nlanczos <size_t>] [-device_type <string>] [-output_file <string>] [-use_double_precision] [-print_vector <size_t>] [-device_id <size_t>]" << std::endl;
            std::cout << "Defaults: natoms = 20, nruns = 1, nwarmup = 1, nisomers = 1, nlanczos = 50, device_type = gpu, output_file = \"output.txt\", use_double_precision = false, print_vector = 0, device_id = 0" << std::endl;
            exit(0);
        }
    }
            if (argc == 1 || std::find(argsVec.begin(), argsVec.end(), "-h") != argsVec.end()) {
        std::cout << "Usage: " << argv[0] << " [-natoms <size_t>] [-nruns <size_t>] [-nwarmup <size_t>] [-nisomers <size_t>] [-nlanczos <size_t>] [-device_type <string>] [-output_file <string>] [-use_double_precision] [-device_id <size_t>] [-print_vector <size_t>]" << std::endl;
        std::cout << "Defaults: natoms = 20, nruns = 1, nwarmup = 1, nisomers = 1, nlanczos = 50, device_type = gpu, output_file = \"output.txt\", use_double_precision = false, print_vector = 0, " << std::endl;
    }
}