#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/gpu/isomer_queue.hh"
#include "fullerenes/gpu/isomer_batch.hh"
#include "fullerenes/gpu/kernels.hh"
#include "fullerenes/gpu/cuda_io.hh"
#include "fullerenes/gpu/cu_array.hh"
#include <filesystem>
#include <numeric>
#include <future>
#include <random>
const std::unordered_map<size_t,size_t> num_fullerenes = {{20,1},{22,0},{24,1},{26,1},{28,2},{30,3},{32,6},{34,6},{36,15},{38,17},{40,40},{42,45},{44,89},{46,116},{48,199},{50,271},{52,437},{54,580},{56,924},{58,1205},{60,1812},{62,2385},{64,3465},{66,4478},{68,6332},{70,8149},{72,11190},{74,14246},{76,19151},{78,24109},{80,31924},{82,39718},{84,51592},{86,63761},{88,81738},{90,99918},{92,126409},{94,153493},{96,191839},{98,231017},{100,285914},{102,341658},{104,419013},{106,497529},{108,604217},{110,713319},{112,860161},{114,1008444},{116,1207119},{118,1408553},{120,1674171},{122,1942929},{124,2295721},{126,2650866},{128,3114236},{130,3580637},{132,4182071},{134,4787715},{136,5566949},{138,6344698},{140,7341204},{142,8339033},{144,9604411},{146,10867631},{148,12469092},{150,14059174},{152,16066025},{154,18060979},{156,20558767},{158,23037594},{160,26142839},{162,29202543},{164,33022573},{166,36798433},{168,41478344},{170,46088157},{172,51809031},{174,57417264},{176,64353269},{178,71163452},{180,79538751},{182,87738311},{184,97841183},{186,107679717},{188,119761075},{190,131561744},{192,145976674},{194,159999462},{196,177175687},{198,193814658},{200,214127742},{202,233846463},{204,257815889},{206,281006325},{208,309273526},{210,336500830},{212,369580714},{214,401535955},{216,440216206},{218,477420176},{220,522599564},{222,565900181},{224,618309598},{226,668662698},{228,729414880},{230,787556069},{232,857934016},{234,925042498},{236,1006016526},{238,1083451816},{240,1176632247},{242,1265323971},{244,1372440782},{246,1474111053},{248,1596482232},{250,1712934069},{252,1852762875},{254,1985250572},{256,2144943655},{258,2295793276},{260,2477017558},{262,2648697036},{264,2854536850},{266,3048609900},{268,3282202941},{270,3501931260},{272,3765465341},{274,4014007928},{276,4311652376},{278,4591045471},{280,4926987377},{282,5241548270},{284,5618445787},{286,5972426835},{288,6395981131},{290,6791769082},{292,7267283603},{294,7710782991},{296,8241719706},{298,8738236515},{300,9332065811},{302,9884604767},{304,10548218751},{306,11164542762},{308,11902015724},{310,12588998862},{312,13410330482},{314,14171344797},{316,15085164571},{318,15930619304},{320,16942010457},{322,17880232383},{324,19002055537},{326,20037346408},{328,21280571390},{330,22426253115},{332,23796620378},{334,25063227406},{336,26577912084},{338,27970034826},{340,29642262229},{342,31177474996},{344,33014225318},{346,34705254287},{348,36728266430},{350,38580626759},{352,40806395661},{354,42842199753},{356,45278616586},{358,47513679057},{360,50189039868},{362,52628839448},{364,55562506886},{366,58236270451},{368,61437700788},{370,64363670678},{372,67868149215},{374,71052718441},{376,74884539987},{378,78364039771},{380,82532990559},{382,86329680991},{384,90881152117},{386,95001297565},{388,99963147805},{390,104453597992},{392,109837310021},{394,114722988623},{396,120585261143},{398,125873325588},{400,132247999328}};
using namespace gpu_kernels;
using namespace cuda_io;
#define SYNC LaunchPolicy::SYNC
#define ASYNC LaunchPolicy::ASYNC

int main(int ac, char** argv){
    int N_start                 = ac > 1 ? strtol(argv[1],0,0) : 20;     // Argument 1: Number of vertices N
    int N_end                   = ac > 2 ? strtol(argv[2],0,0) : N_start;     // Argument 1: Number of vertices N
    auto N_runs = 3;

    ofstream out_file("IsomerspaceOpt_V7_" + to_string(N_end) + ".txt");
    ofstream out_file_std("IsomerspaceOpt_V7_STD_" + to_string(N_end) + ".txt");
    for(int N = N_start; N < N_end + 1;  N+=2){
        if(N == 22) continue;
        auto n_fullerenes = num_fullerenes.find(N)->second;
        Graph G;
        auto Nf = N/2 +2;
        G.neighbours = neighbours_t(Nf, std::vector<node_t>(6));
        G.N = Nf;
        auto bucky = BuckyGen::start(N, false, false);

        std::string path = "isomerspace_samples/dual_layout_" + to_string(N) + "_seed_42";
        auto fsize = std::filesystem::file_size(path);
        auto n_samples = fsize / (Nf * 6 * sizeof(device_node_t));
        ifstream in_file(path,std::ios::binary);
        std::vector<device_node_t> dual_neighbours(n_samples * Nf * 6);
        in_file.read((char*)dual_neighbours.data(),n_samples * Nf * 6 * sizeof(device_node_t));
        auto optimal_batch_size = isomerspace_forcefield::optimal_batch_size(N);
        auto sample_size = min(optimal_batch_size,(int)n_samples);
        if (n_fullerenes < optimal_batch_size*2){
            sample_size = max((int)n_fullerenes/2,1);
        }else if(n_fullerenes >= optimal_batch_size*8){
            sample_size = optimal_batch_size*4;
        }else if (n_fullerenes >= optimal_batch_size*6){
            sample_size = optimal_batch_size*3;
        }else if (n_fullerenes >= optimal_batch_size*4){
            sample_size = optimal_batch_size*2;
        } else if (n_fullerenes >= optimal_batch_size*2){
            sample_size = optimal_batch_size;
        }

        std::vector<int> random_IDs(n_samples);
        std::iota(random_IDs.begin(), random_IDs.end(), 0);
        std::shuffle(random_IDs.begin(), random_IDs.end(), std::mt19937{42});
        auto id_range_end = min((int)n_samples, sample_size);
        std::vector<int> id_subset(random_IDs.begin(), random_IDs.begin()+n_samples);
        auto
            T_ends    = std::vector<std::chrono::nanoseconds>(N_runs, chrono::nanoseconds(1)),
            T_par    = std::vector<std::chrono::nanoseconds>(N_runs, chrono::nanoseconds(1)),
            T_io     = std::vector<std::chrono::nanoseconds>(N_runs, chrono::nanoseconds(1));
        
        auto finished_fullerenes = 0;
        //ACTUAL PIPELINE CODE
        for(int i = 0; i < N_runs; i++){
        finished_fullerenes = 0;
        LaunchCtx insert0_ctx(0);
        LaunchCtx insert1_ctx(1);
        LaunchCtx device0_ctx(0);
        LaunchCtx device1_ctx(1);

        IsomerQueue input0_queue(N, 0);
        IsomerQueue input1_queue(N, 1);
        IsomerQueue opt0_queue(N, 0);
        IsomerQueue opt1_queue(N, 1);
        IsomerQueue output0_queue(N,0);
        IsomerQueue output1_queue(N,1);
        input0_queue.resize(sample_size*2);
        input1_queue.resize(sample_size*2);
        output0_queue.resize(sample_size);
        output1_queue.resize(sample_size);
        opt0_queue.resize(n_samples*2);
        opt1_queue.resize(n_samples*2);

        IsomerBatch input0(N, sample_size*2, DEVICE_BUFFER, 0);
        IsomerBatch input1(N, sample_size*2, DEVICE_BUFFER, 1);
        IsomerBatch opt0(N, sample_size, DEVICE_BUFFER, 0);
        IsomerBatch opt1(N, sample_size, DEVICE_BUFFER, 1);
        
        int I_async = 0;
        auto generate_isomers = [&](int M){
            if(I_async == min(n_fullerenes,n_samples*10)) return false;
            for (int i = 0; i < M; i++){
                if (I_async < min(n_fullerenes,n_samples*10)){
                    for (size_t j = 0; j < Nf; j++){
                        G.neighbours[j].clear();
                        for (size_t k = 0; k < 6; k++) {
                            auto u = dual_neighbours[id_subset[I_async%n_samples]*Nf*6 + j*6 +k];
                            if(u != UINT16_MAX) G.neighbours[j].push_back(u);
                        }
                    }
                    input0_queue.insert(G,I_async, insert0_ctx, ASYNC);
                    input1_queue.insert(G,I_async, insert1_ctx, ASYNC);
                    I_async++;
                } 
            }
            
            input0_queue.refill_batch(input0, insert0_ctx, ASYNC);
            input1_queue.refill_batch(input1, insert1_ctx, ASYNC);
            isomerspace_dual::dualise(input0, insert0_ctx, ASYNC);
            isomerspace_dual::dualise(input1, insert1_ctx, ASYNC);
            isomerspace_tutte::tutte_layout(input0, 1000000, insert0_ctx, ASYNC);
            isomerspace_tutte::tutte_layout(input1, 1000000, insert1_ctx, ASYNC);
            isomerspace_X0::zero_order_geometry(input0, 4.0, insert0_ctx, ASYNC);
            isomerspace_X0::zero_order_geometry(input1, 4.0, insert1_ctx, ASYNC);
            insert0_ctx.wait(); insert1_ctx.wait();
            
            return I_async < min(size_t(ceil(n_fullerenes/2)),n_samples*10);
        };
        auto T1 = chrono::high_resolution_clock::now();
        generate_isomers(sample_size*2);
        
        opt0_queue.insert(input0, device0_ctx, ASYNC);
        opt1_queue.insert(input1, device1_ctx, ASYNC);
        opt0_queue.refill_batch(opt0, device0_ctx, ASYNC);
        opt1_queue.refill_batch(opt1, device1_ctx, ASYNC);
        device0_ctx.wait(); device1_ctx.wait();
        T_ends[i] += T1 - chrono::high_resolution_clock::now();
        bool more_to_do = true;
        bool more_to_generate = false;
        auto step = max(1, (int)N/2);
        
        while (more_to_do){
            bool optimise_more = true;
            auto generate_handle = std::async(std::launch::async,generate_isomers, opt0.isomer_capacity*2);
            while(optimise_more){
                auto T2 = chrono::high_resolution_clock::now();
                isomerspace_forcefield::optimise<PEDERSEN>(opt0,step, N*5, device0_ctx, ASYNC);
                isomerspace_forcefield::optimise<PEDERSEN>(opt1,step, N*5, device1_ctx, ASYNC);
                output0_queue.push(opt0, device0_ctx, ASYNC);
                output1_queue.push(opt1, device1_ctx, ASYNC);
                opt0_queue.refill_batch(opt0, device0_ctx, ASYNC);
                opt1_queue.refill_batch(opt1, device1_ctx, ASYNC);
                device0_ctx.wait(); device1_ctx.wait();
                T_par[i] += chrono::high_resolution_clock::now() - T2;
                finished_fullerenes += output0_queue.get_size() + output1_queue.get_size();
                output0_queue.clear(device0_ctx);
                output1_queue.clear(device1_ctx);
                optimise_more = opt0_queue.get_size() >= opt0.isomer_capacity;
            }
            auto T3 = chrono::high_resolution_clock::now();
            generate_handle.wait();
            more_to_generate = generate_handle.get();
            opt0_queue.insert(input0, device0_ctx, ASYNC);
            opt1_queue.insert(input1, device1_ctx, ASYNC);
            device0_ctx.wait(); device1_ctx.wait();
            finished_fullerenes += output0_queue.get_size() + output1_queue.get_size();
            auto T4 = chrono::high_resolution_clock::now();
            if(more_to_generate) T_par[i] += T4 - T3;
            if(!more_to_generate){
                while(opt0_queue.get_size() > 0){
                    isomerspace_forcefield::optimise<PEDERSEN>(opt0, step, N*5, device0_ctx, ASYNC);
                    isomerspace_forcefield::optimise<PEDERSEN>(opt1, step, N*5, device1_ctx, ASYNC);
                    output0_queue.push(opt0, device0_ctx, ASYNC);
                    output1_queue.push(opt1, device1_ctx, ASYNC);
                    opt0_queue.refill_batch(opt0, device0_ctx, ASYNC);
                    opt1_queue.refill_batch(opt1, device1_ctx, ASYNC);
                    device0_ctx.wait(); device1_ctx.wait();
                }
                for(int i = 0;  i <  N*5; i += step){
                    isomerspace_forcefield::optimise<PEDERSEN>(opt0, step, N*5, device0_ctx, ASYNC);
                    isomerspace_forcefield::optimise<PEDERSEN>(opt1, step, N*5, device1_ctx, ASYNC);
                }
                output0_queue.push(opt0, device0_ctx, ASYNC);
                output1_queue.push(opt1, device1_ctx, ASYNC);
                device0_ctx.wait(); device1_ctx.wait();
                more_to_do = false;
            }
            T_ends[i] += chrono::high_resolution_clock::now() - T4;

        }
        }
        if(n_fullerenes > n_samples*10){ 
            auto total = (float)(mean(T_io)/1ns + mean(T_par)/1ns);
            std::cout << std::fixed << std::setprecision(2) << N << ", "<< finished_fullerenes << ", " << (mean(T_par)/1ns)/total*100. << "%, " << (mean(T_io)/1ns)/total*100. << "%, " << (float)(mean(T_io)/1us+mean(T_par)/1us)/finished_fullerenes << "us/isomer\n";
            out_file << N << ", "<< n_fullerenes << ", "<< finished_fullerenes << ", " << mean(T_ends) /1ns << ", " << mean(T_par)/1ns << ", " << mean(T_io)/1ns << "\n";
            out_file_std << N << ", "<< n_fullerenes << ", "<< finished_fullerenes << ", " << sdev(T_ends) /1ns << ", " << sdev(T_par)/1ns << ", " << sdev(T_io)/1ns << "\n";}
        else{
            auto total = (float)(mean(T_io)/1ns + mean(T_par)/1ns + mean(T_ends)/1ns);
            std::cout << std::fixed << std::setprecision(2) << N << ", "<< n_fullerenes << ", " << (mean(T_par)/1ns)/total*100. << "%, " << (mean(T_io)/1ns)/total*100. << "%, " << (float)(mean(T_io)/1us+mean(T_par)/1us + mean(T_ends)/1us)/n_fullerenes << "us/isomer\n";
            out_file << N << ", "<< n_fullerenes << ", "<< n_fullerenes << ", " << mean(T_ends) /1ns << ", " << mean(T_par)/1ns + mean(T_ends)/1ns << ", " << mean(T_io)/1ns << "\n";
            out_file_std << N << ", "<< n_fullerenes << ", "<< n_fullerenes << ", " << sdev(T_ends) /1ns << ", " << sdev(T_par)/1ns + sdev(T_ends)/1ns << ", " << sdev(T_io)/1ns << "\n";}

        std::cout << (float)finished_fullerenes / (float)(mean(T_par)/1ms) << std::endl;
    }
    LaunchCtx::clear_allocations();
}