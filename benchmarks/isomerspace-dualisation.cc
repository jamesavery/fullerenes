#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/polyhedron.hh"
#include <chrono>
#include <fstream>
#include <iomanip> 		// Needed for std::setprecision

const std::unordered_map<size_t,size_t> num_fullerenes = {{20,1},{22,0},{24,1},{26,1},{28,2},{30,3},{32,6},{34,6},{36,15},{38,17},{40,40},{42,45},{44,89},{46,116},{48,199},{50,271},{52,437},{54,580},{56,924},{58,1205},{60,1812},{62,2385},{64,3465},{66,4478},{68,6332},{70,8149},{72,11190},{74,14246},{76,19151},{78,24109},{80,31924},{82,39718},{84,51592},{86,63761},{88,81738},{90,99918},{92,126409},{94,153493},{96,191839},{98,231017},{100,285914},{102,341658},{104,419013},{106,497529},{108,604217},{110,713319},{112,860161},{114,1008444},{116,1207119},{118,1408553},{120,1674171},{122,1942929},{124,2295721},{126,2650866},{128,3114236},{130,3580637},{132,4182071},{134,4787715},{136,5566949},{138,6344698},{140,7341204},{142,8339033},{144,9604411},{146,10867631},{148,12469092},{150,14059174},{152,16066025},{154,18060979},{156,20558767},{158,23037594},{160,26142839},{162,29202543},{164,33022573},{166,36798433},{168,41478344},{170,46088157},{172,51809031},{174,57417264},{176,64353269},{178,71163452},{180,79538751},{182,87738311},{184,97841183},{186,107679717},{188,119761075},{190,131561744},{192,145976674},{194,159999462},{196,177175687},{198,193814658},{200,214127742},{202,233846463},{204,257815889},{206,281006325},{208,309273526},{210,336500830},{212,369580714},{214,401535955},{216,440216206},{218,477420176},{220,522599564},{222,565900181},{224,618309598},{226,668662698},{228,729414880},{230,787556069},{232,857934016},{234,925042498},{236,1006016526},{238,1083451816},{240,1176632247},{242,1265323971},{244,1372440782},{246,1474111053},{248,1596482232},{250,1712934069},{252,1852762875},{254,1985250572},{256,2144943655},{258,2295793276},{260,2477017558},{262,2648697036},{264,2854536850},{266,3048609900},{268,3282202941},{270,3501931260},{272,3765465341},{274,4014007928},{276,4311652376},{278,4591045471},{280,4926987377},{282,5241548270},{284,5618445787},{286,5972426835},{288,6395981131},{290,6791769082},{292,7267283603},{294,7710782991},{296,8241719706},{298,8738236515},{300,9332065811},{302,9884604767},{304,10548218751},{306,11164542762},{308,11902015724},{310,12588998862},{312,13410330482},{314,14171344797},{316,15085164571},{318,15930619304},{320,16942010457},{322,17880232383},{324,19002055537},{326,20037346408},{328,21280571390},{330,22426253115},{332,23796620378},{334,25063227406},{336,26577912084},{338,27970034826},{340,29642262229},{342,31177474996},{344,33014225318},{346,34705254287},{348,36728266430},{350,38580626759},{352,40806395661},{354,42842199753},{356,45278616586},{358,47513679057},{360,50189039868},{362,52628839448},{364,55562506886},{366,58236270451},{368,61437700788},{370,64363670678},{372,67868149215},{374,71052718441},{376,74884539987},{378,78364039771},{380,82532990559},{382,86329680991},{384,90881152117},{386,95001297565},{388,99963147805},{390,104453597992},{392,109837310021},{394,114722988623},{396,120585261143},{398,125873325588},{400,132247999328}};
using namespace chrono;
using namespace chrono_literals;

#include "fullerenes/isomer_queue.hh"
#include "fullerenes/device_io.hh"
#include "fullerenes/gpu/kernels.hh"
#include "fullerenes/gpu/benchmark_functions.hh"
#include "numeric"
#include "random"

using namespace gpu_kernels;

#define VALIDATION 0

int main(int argc, char** argv){

    size_t N_start                = argc > 1 ? strtol(argv[1],0,0) : (size_t)20;        // Argument 1: Start of range of N
    size_t N_limit                = argc > 2 ? strtol(argv[2],0,0) : N_start;           // Argument 2: End of range of N
    size_t N_runs                 = argc > 3 ? strtol(argv[3],0,0) : 3;                 // Argument 3: Number of times to run experiment
    size_t warmup                 = argc > 4 ? strtol(argv[4],0,0) : 1;                 // Argument 4: Seconds to warmup GPU
    size_t batch_coeff            = argc > 5 ? strtol(argv[5],0,0) : 1;                 // Argument 5: Batch size coefficient
    size_t dual_version           = argc > 6 ? strtol(argv[6],0,0) : 0;                 // Argument 6: Dual version to use
    size_t Ngpu                 = argc > 7 ? strtol(argv[7],0,0) : LaunchCtx::get_device_count();                 // Argument 7: Number of GPUs to use        

    ofstream out_file   ("DualBenchmark" + to_string(N_limit) + ".txt");
    ofstream out_std    ("DualBenchmark_STD_" + to_string(N_limit) + ".txt");
    using namespace std::chrono;
    //if(dual_version >= 2) cuda_benchmark::warmup_kernel(warmup*200s);
    LaunchCtx ctx(0);
    for (size_t N = N_start; N < N_limit+1; N+=2)
    {   
        if(N == 22) continue;
        
        auto Nf = N/2 + 2;
        auto Nd = LaunchCtx::get_device_count();
        int mini_batch = isomerspace_dual::optimal_batch_size_2(Nf,0);
        auto batch_size = mini_batch*batch_coeff;
        if(dual_version >= 2) {mini_batch = batch_coeff; batch_size = mini_batch;}
        auto n_fullerenes = (int)num_fullerenes.find(N)->second;
        
        if (N == 22) continue;
        
        bool more_to_generate = true;
        
        std::vector<double> T_duals(N_runs), T_timing_diff(N_runs);


        FullereneDual F;
        F.neighbours = neighbours_t(Nf, std::vector<node_t>(6));
        F.N = Nf;


        auto path = "isomerspace_samples/dual_layout_" + to_string(N) + "_seed_42";
        ifstream isomer_sample(path,std::ios::binary);
        auto fsize = file_size(path);
        std::vector<device_node_t> input_buffer(fsize/sizeof(device_node_t));
        auto available_samples = fsize / (Nf*6*sizeof(device_node_t));
        isomer_sample.read(reinterpret_cast<char*>(input_buffer.data()), Nf*6*sizeof(device_node_t)*available_samples);

        std::vector<int> IDs(available_samples);
        std::iota(IDs.begin(), IDs.end(), 0);
        
        auto finished_isomers = 0;
        
            finished_isomers = 0;
            IsomerBatch<CPU> HBs[Nd] = {IsomerBatch<CPU>(N,mini_batch), IsomerBatch<CPU>(N, mini_batch)};
            IsomerBatch<GPU> DBs[Nd] = {IsomerBatch<GPU>(N,batch_size,0), IsomerBatch<GPU>(N, batch_size,1)};
            IsomerBatch<GPU> TempBatch[Nd] = {IsomerBatch<GPU>(N,mini_batch,0), IsomerBatch<GPU>(N, mini_batch,1)};
            #if VALIDATION
                IsomerBatch<CPU> ValBatch[Nd] = {IsomerBatch<CPU>(N,batch_size), IsomerBatch<CPU>(N, batch_size)};
                IsomerBatch<CPU> TestBatch[Nd] = {IsomerBatch<CPU>(N,batch_size), IsomerBatch<CPU>(N, batch_size)};
                IsomerBatch<CPU> TBs[Nd] = {IsomerBatch<CPU>(N,mini_batch), IsomerBatch<CPU>(N, mini_batch)};
            #endif
            device_io::IsomerQueue Q_test[Nd] = {device_io::IsomerQueue(N,0), device_io::IsomerQueue(N,1)};
            Q_test[0].resize(batch_size); Q_test[1].resize(batch_size);

            #if VALIDATION
                device_io::IsomerQueue Q_val[Nd] = {device_io::IsomerQueue(N,0), device_io::IsomerQueue(N,1)};
                Q_val[0].resize(batch_size); Q_val[1].resize(batch_size);
            #endif
            LaunchCtx ctxs[Nd] = {LaunchCtx(0),LaunchCtx(1)};

            for (int i = 0; i < mini_batch; i++){
                for (size_t j = 0; j < Nf; j++){
                    F.neighbours[j].clear();
                    for (size_t k = 0; k < 6; k++) {
                        auto u = input_buffer[IDs[i%available_samples]*Nf*6 + j*6 +k];
                        if(u != UINT16_MAX) F.neighbours[j].push_back(u);
                    }
                }
                for(auto j = 0; j < Nd; j++){
                    HBs[j].append(Graph(F), i);
                }
                #if VALIDATION

                    F.update();
                    auto FD = F.dual_graph();
                    for (auto j = 0; j < Nd; j++)
                    {
                        TBs[j].append(FD,i, false);
                    }
                #endif
            }
            for (size_t i = 0; i < Nd; i++) device_io::copy(TempBatch[i],HBs[i]);

            for (size_t i = 0; i < batch_coeff; i++) {
                for(auto j = 0; j < Nd; j++) {device_io::copy(TempBatch[j], HBs[j], ctxs[j], LaunchPolicy::ASYNC);Q_test[j].insert(TempBatch[j],ctxs[j],LaunchPolicy::ASYNC);} 
            }
            for(auto j = 0; j < Nd; j++) ctxs[j].wait();
            #if VALIDATION    
                for (size_t i = 0; i < Nd; i++) device_io::copy(TempBatch[i],TBs[i]);

                for (size_t i = 0; i < batch_coeff; i++) {
                    for(auto j = 0; j < Nd; j++) {device_io::copy(TempBatch[j], TBs[j], ctxs[j], LaunchPolicy::ASYNC);Q_val[j].insert(TempBatch[j],ctxs[j],LaunchPolicy::ASYNC);}
                }
            #endif
            
            //std::cout << "Queue Size: " << Q_test[0].get_size() << std::endl;
            for (size_t i = 0; i < Nd; i++) device_io::copy(DBs[i],Q_test[i].device_batch);
            for (size_t i = 0; i < Nd; i++) device_io::copy(Q_test[i].host_batch, Q_test[i].device_batch);
            
            
            for (size_t l = 0; l < N_runs + 1; l++){
                
                if (dual_version != 2) {for (size_t i = 0; i < Ngpu; i++) ctxs[(i+1)%Nd].start_timer();}
                auto T0 = steady_clock::now();
                switch (dual_version)
                {
                case 0:
                    for (size_t i = 0; i < Ngpu; i++) isomerspace_dual::dualise(DBs[(i+1)%Nd], ctxs[(i+1)%Nd], LaunchPolicy::ASYNC);
                    break;
                case 1:
                    for (size_t i = 0; i < Ngpu; i++) isomerspace_dual::dualise_2(DBs[(i+1)%Nd], ctxs[(i+1)%Nd], LaunchPolicy::ASYNC);
                    break;
                case 2:
                    for (size_t i = 0; i < Ngpu; i++) isomerspace_dual::dualise_3(Q_test[(i+1)%Nd].host_batch);
                    break;
                case 3:
                    for (size_t i = 0; i < Ngpu; i++) isomerspace_dual::dualise_4(Q_test[(i+1)%Nd].host_batch);
                    break;
                default:
                    break;
                }
                std::vector<nanoseconds> times(Nd);
                if (dual_version < 2) {for (size_t i = 0; i < Ngpu; i++) times[(i+1)%Nd] = ctxs[(i+1)%Nd].stop_timer();}
                auto T1 = steady_clock::now();
                //std::cout << "Dualisation time: " << duration<double, std::nano>(T1 - T0).count() << std::endl;
                if(l != 0) T_duals[l-1] = dual_version>=2 ? duration<double, std::nano>(T1 - T0).count() : duration<double, std::nano>(*max_element(times.begin(), times.end())).count();
                if(l != 0 && dual_version < 2) T_timing_diff[l-1]  = duration<double, std::nano>(T1 - T0).count() - duration<double, std::nano>(*max_element(times.begin(), times.end())).count();
            }
            //Q_test[0].host_batch.set_print_verbose();
            //std::cout << Q_test[0].host_batch << std::endl;
            #if VALIDATION
            //for(auto i = 0; i < Nd; i++) {ValBatch[i].set_print_verbose(); TestBatch[i].set_print_verbose();}
            //for(auto i = 0; i < Nd; i++) std::cout << ValBatch[i] << std::endl;
            //for(auto i = 0; i < Nd; i++) std::cout << TestBatch[i] << std::endl;
                for(auto i = 0; i < Nd; i++) device_io::copy(ValBatch[i], Q_val[i].device_batch);
                if (dual_version != 2) {for(auto i = 0; i < Nd; i++) device_io::copy(TestBatch[i], DBs[i]);}
                for(auto i = 0; i < Nd; i++) assert(ValBatch[i] == TestBatch[i]);
            #endif



        using namespace device_io;
        //Print out runtimes in us per isomer:
        std::cout << N << "  Dual: " << std::fixed << std::setprecision(2) << mean(T_duals)/float(batch_size*Ngpu) << " +- " << sdev(T_duals)/float(batch_size*Ngpu) <<  " Timing Diff: " << mean(T_timing_diff)/float(batch_size*Ngpu) << " +- " << sdev(T_timing_diff)/float(batch_size*Ngpu) << std::endl;

        //Print out what fraction of the runtime that each component took:
        //Remove the worst outlier
        std::sort(T_duals.begin(), T_duals.end());
        T_duals.pop_back();
        out_file << N << ", "<< batch_size*Ngpu << ", " << long(mean(T_duals)) << std::endl;
        out_std << N << ", "<< batch_size*Ngpu << ", "  << long(sdev(T_duals)) << std::endl;
     }
    
}
