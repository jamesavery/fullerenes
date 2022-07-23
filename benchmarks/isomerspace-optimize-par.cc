#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/polyhedron.hh"
#include <chrono>
#include <fstream>

const std::unordered_map<size_t,size_t> num_fullerenes = {{20,1},{22,0},{24,1},{26,1},{28,2},{30,3},{32,6},{34,6},{36,15},{38,17},{40,40},{42,45},{44,89},{46,116},{48,199},{50,271},{52,437},{54,580},{56,924},{58,1205},{60,1812},{62,2385},{64,3465},{66,4478},{68,6332},{70,8149},{72,11190},{74,14246},{76,19151},{78,24109},{80,31924},{82,39718},{84,51592},{86,63761},{88,81738},{90,99918},{92,126409},{94,153493},{96,191839},{98,231017},{100,285914},{102,341658},{104,419013},{106,497529},{108,604217},{110,713319},{112,860161},{114,1008444},{116,1207119},{118,1408553},{120,1674171},{122,1942929},{124,2295721},{126,2650866},{128,3114236},{130,3580637},{132,4182071},{134,4787715},{136,5566949},{138,6344698},{140,7341204},{142,8339033},{144,9604411},{146,10867631},{148,12469092},{150,14059174},{152,16066025},{154,18060979},{156,20558767},{158,23037594},{160,26142839},{162,29202543},{164,33022573},{166,36798433},{168,41478344},{170,46088157},{172,51809031},{174,57417264},{176,64353269},{178,71163452},{180,79538751},{182,87738311},{184,97841183},{186,107679717},{188,119761075},{190,131561744},{192,145976674},{194,159999462},{196,177175687},{198,193814658},{200,214127742},{202,233846463},{204,257815889},{206,281006325},{208,309273526},{210,336500830},{212,369580714},{214,401535955},{216,440216206},{218,477420176},{220,522599564},{222,565900181},{224,618309598},{226,668662698},{228,729414880},{230,787556069},{232,857934016},{234,925042498},{236,1006016526},{238,1083451816},{240,1176632247},{242,1265323971},{244,1372440782},{246,1474111053},{248,1596482232},{250,1712934069},{252,1852762875},{254,1985250572},{256,2144943655},{258,2295793276},{260,2477017558},{262,2648697036},{264,2854536850},{266,3048609900},{268,3282202941},{270,3501931260},{272,3765465341},{274,4014007928},{276,4311652376},{278,4591045471},{280,4926987377},{282,5241548270},{284,5618445787},{286,5972426835},{288,6395981131},{290,6791769082},{292,7267283603},{294,7710782991},{296,8241719706},{298,8738236515},{300,9332065811},{302,9884604767},{304,10548218751},{306,11164542762},{308,11902015724},{310,12588998862},{312,13410330482},{314,14171344797},{316,15085164571},{318,15930619304},{320,16942010457},{322,17880232383},{324,19002055537},{326,20037346408},{328,21280571390},{330,22426253115},{332,23796620378},{334,25063227406},{336,26577912084},{338,27970034826},{340,29642262229},{342,31177474996},{344,33014225318},{346,34705254287},{348,36728266430},{350,38580626759},{352,40806395661},{354,42842199753},{356,45278616586},{358,47513679057},{360,50189039868},{362,52628839448},{364,55562506886},{366,58236270451},{368,61437700788},{370,64363670678},{372,67868149215},{374,71052718441},{376,74884539987},{378,78364039771},{380,82532990559},{382,86329680991},{384,90881152117},{386,95001297565},{388,99963147805},{390,104453597992},{392,109837310021},{394,114722988623},{396,120585261143},{398,125873325588},{400,132247999328}};
using namespace chrono;
using namespace chrono_literals;

#include "fullerenes/gpu/batch_queue.hh"
#include "fullerenes/gpu/cuda_io.hh"
#include "fullerenes/gpu/kernels.hh"
#include "numeric"
#include "random"
#include "filesystem"
using namespace gpu_kernels;

int main(int argc, char** argv){

    size_t N_start                = argc > 1 ? strtol(argv[1],0,0) : (size_t)20;     // Argument 1: Number of vertices N
    size_t N_limit                = argc > 2 ? strtol(argv[2],0,0) : N_start;     // Argument 1: Number of vertices N
    
    ofstream out_file("ParBenchmark_" + to_string(N_limit) + ".txt");
    ofstream out_std("ParBenchmark_STD_" + to_string(N_limit) + ".txt");
    out_file << "Generate, Samples, Update, Dual, Tutte, X0, Optimize \n";
    for (size_t N = N_start; N < N_limit+1; N+=2)
    {   
        BuckyGen::buckygen_queue Q = BuckyGen::start(N,false,false);  

        FullereneDual F;

        auto batch_size = isomerspace_forcefield::optimal_batch_size(N);
        auto sample_size = min(batch_size*1,(int)num_fullerenes.find(N)->second);
        if (N == 22) continue;
        std::cout << sample_size << endl;
        IsomerBatch batch0(N,sample_size,DEVICE_BUFFER);
        IsomerBatch batch1(N,sample_size,DEVICE_BUFFER);
        IsomerBatch h_batch(N,sample_size,HOST_BUFFER);
        cuda_io::IsomerQueue isomer_q(N);
        cuda_io::IsomerQueue isomer_q_cubic(N);
        bool more_to_generate = true;

        auto T0 = high_resolution_clock::now();
        auto

            Tgen_mean    = high_resolution_clock::now()-T0,
            Tgen_std    = high_resolution_clock::now()-T0,
            Tupdate_mean = high_resolution_clock::now()-T0,
            Tupdate_std = high_resolution_clock::now()-T0,
            Tdual_mean   = high_resolution_clock::now()-T0,
            Tdual_std   = high_resolution_clock::now()-T0,    
            Ttutte_mean  = high_resolution_clock::now()-T0,
            Ttutte_std  = high_resolution_clock::now()-T0,
            TX0_mean     = high_resolution_clock::now()-T0,
            TX0_std     = high_resolution_clock::now()-T0,
            Topt_mean    = high_resolution_clock::now()-T0,
            Topt_std    = high_resolution_clock::now()-T0;
        
        std::vector<std::chrono::nanoseconds> 
            T_gens(10),
            T_duals(10),
            T_tuttes(10),
            T_X0s(10),
            T_opts(10);

        auto path = "isomerspace_samples/" + to_string(N) + "_seed_42";
        ifstream isomer_sample(path,std::ios::binary);
        auto fsize = std::filesystem::file_size(path);
        std::vector<device_node_t> input_buffer(fsize/sizeof(device_node_t));
        isomer_sample.read(reinterpret_cast<char*>(input_buffer.data()), fsize);
        auto available_samples = fsize / (N*3*sizeof(device_node_t));

        std::vector<int> random_IDs(available_samples);
        std::iota(random_IDs.begin(), random_IDs.end(), 0);
        std::shuffle(random_IDs.begin(), random_IDs.end(), std::mt19937{42});
        std::vector<int> IDs(random_IDs.begin(), random_IDs.begin()+sample_size*T_gens.size());

        auto I = 0;
        while (more_to_generate && I < sample_size*T_gens.size()){
                auto T1 = high_resolution_clock::now();
            more_to_generate &= BuckyGen::next_fullerene(Q, F);
            if (!more_to_generate){break;}
                auto T2 = high_resolution_clock::now(); T_gens[I / sample_size] += T2 - T1; 
            isomer_q.insert(Graph(F), I);
            PlanarGraph PG(F);
            PG.neighbours = std::vector(N, std::vector<node_t>(3));
                auto T3 = high_resolution_clock::now(); 
                for(int j = 0; j < N; j++){
                    for (int k = 0; k < 3; k++){
                        PG.neighbours[j][k] = input_buffer[N*IDs[I]*3 + j*3+k];
                    }
                } 
            isomer_q_cubic.insert(PG,I,LaunchCtx(),LaunchPolicy::SYNC,false);
            ++I;
        }

        auto mean = [&](std::vector<std::chrono::nanoseconds> &input){
            auto size = min((int)input.size() , (int)num_fullerenes.find(N)->second/sample_size + 1);
            auto result = std::chrono::nanoseconds(0);
            for (int i = 0; i < size; i++){
                result += input[i];
            }
            return result / size;
        };

        auto standard_deviation = [&](std::vector<std::chrono::nanoseconds> &input){
            auto size = min((int)input.size() , (int)num_fullerenes.find(N)->second/sample_size + 1);
            auto mean_ = mean(input);
            auto result = std::chrono::nanoseconds(0);
            for (int i = 0; i < size; i++){
                result += (input[i] - mean_) *  (input[i] - mean_).count();
            }
            return std::chrono::nanoseconds((int)std::sqrt( (result / size).count()));
        };
        for (size_t i = 0; i < T_gens.size(); i++)
        {
            isomer_q.refill_batch(batch0);
            isomer_q_cubic.refill_batch(batch1);
            auto TDual = isomerspace_dual::time_spent();
            isomerspace_dual::cubic_layout(batch0);
            auto TTutte = isomerspace_tutte::time_spent(); T_duals[i] = isomerspace_dual::time_spent() - TDual;
            isomerspace_tutte::tutte_layout(batch1);
            auto TX0 = isomerspace_X0::time_spent(); T_tuttes[i] = isomerspace_tutte::time_spent() - TTutte;
            isomerspace_X0::zero_order_geometry(batch1, 4.0);
            T_X0s[i] = isomerspace_X0::time_spent() - TX0;
            cuda_io::reset_convergence_statuses(batch1);
            auto TFF = isomerspace_forcefield::time_spent();
            isomerspace_forcefield::optimize_batch(batch1,N*4,N*4);
            T_opts[i] = isomerspace_forcefield::time_spent() - TFF;
        }

        Tgen_mean = mean(T_gens);
        Tgen_std = standard_deviation(T_gens);
        Tdual_mean = mean(T_duals);
        Tdual_std = standard_deviation(T_duals);
        Ttutte_mean = mean(T_tuttes);
        Ttutte_std = standard_deviation(T_tuttes);
        TX0_mean = mean(T_X0s);
        TX0_std = standard_deviation(T_X0s);
        Topt_mean = mean(T_opts);
        Topt_std = standard_deviation(T_opts);

        out_file << N << ", "<< sample_size << ", " << Tgen_mean/1ns << ", " << Tdual_mean/1ns <<", " <<  TX0_mean/1ns <<", " << Ttutte_mean/1ns<< ", " << Topt_mean/1ns << "\n";
        out_std << N << ", "<< sample_size << ", " << Tgen_std/1ns << ", "  << Tdual_std/1ns <<", " <<  TX0_std/1ns <<", " << Ttutte_std/1ns<< ", " << Topt_std/1ns << "\n";
     }

}