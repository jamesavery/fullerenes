#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/gpu/cuda_definitions.h"
#include <chrono>
#include <fstream>
#include "filesystem"
#include "random"
#include "numeric"
#include "fullerenes/gpu/isomer_queue.hh"
#include "fullerenes/gpu/cuda_io.hh"
#include "fullerenes/gpu/kernels.hh"
#include <stdio.h>

const std::unordered_map<size_t,size_t> num_fullerenes = {{20,1},{22,0},{24,1},{26,1},{28,2},{30,3},{32,6},{34,6},{36,15},{38,17},{40,40},{42,45},{44,89},{46,116},{48,199},{50,271},{52,437},{54,580},{56,924},{58,1205},{60,1812},{62,2385},{64,3465},{66,4478},{68,6332},{70,8149},{72,11190},{74,14246},{76,19151},{78,24109},{80,31924},{82,39718},{84,51592},{86,63761},{88,81738},{90,99918},{92,126409},{94,153493},{96,191839},{98,231017},{100,285914},{102,341658},{104,419013},{106,497529},{108,604217},{110,713319},{112,860161},{114,1008444},{116,1207119},{118,1408553},{120,1674171},{122,1942929},{124,2295721},{126,2650866},{128,3114236},{130,3580637},{132,4182071},{134,4787715},{136,5566949},{138,6344698},{140,7341204},{142,8339033},{144,9604411},{146,10867631},{148,12469092},{150,14059174},{152,16066025},{154,18060979},{156,20558767},{158,23037594},{160,26142839},{162,29202543},{164,33022573},{166,36798433},{168,41478344},{170,46088157},{172,51809031},{174,57417264},{176,64353269},{178,71163452},{180,79538751},{182,87738311},{184,97841183},{186,107679717},{188,119761075},{190,131561744},{192,145976674},{194,159999462},{196,177175687},{198,193814658},{200,214127742},{202,233846463},{204,257815889},{206,281006325},{208,309273526},{210,336500830},{212,369580714},{214,401535955},{216,440216206},{218,477420176},{220,522599564},{222,565900181},{224,618309598},{226,668662698},{228,729414880},{230,787556069},{232,857934016},{234,925042498},{236,1006016526},{238,1083451816},{240,1176632247},{242,1265323971},{244,1372440782},{246,1474111053},{248,1596482232},{250,1712934069},{252,1852762875},{254,1985250572},{256,2144943655},{258,2295793276},{260,2477017558},{262,2648697036},{264,2854536850},{266,3048609900},{268,3282202941},{270,3501931260},{272,3765465341},{274,4014007928},{276,4311652376},{278,4591045471},{280,4926987377},{282,5241548270},{284,5618445787},{286,5972426835},{288,6395981131},{290,6791769082},{292,7267283603},{294,7710782991},{296,8241719706},{298,8738236515},{300,9332065811},{302,9884604767},{304,10548218751},{306,11164542762},{308,11902015724},{310,12588998862},{312,13410330482},{314,14171344797},{316,15085164571},{318,15930619304},{320,16942010457},{322,17880232383},{324,19002055537},{326,20037346408},{328,21280571390},{330,22426253115},{332,23796620378},{334,25063227406},{336,26577912084},{338,27970034826},{340,29642262229},{342,31177474996},{344,33014225318},{346,34705254287},{348,36728266430},{350,38580626759},{352,40806395661},{354,42842199753},{356,45278616586},{358,47513679057},{360,50189039868},{362,52628839448},{364,55562506886},{366,58236270451},{368,61437700788},{370,64363670678},{372,67868149215},{374,71052718441},{376,74884539987},{378,78364039771},{380,82532990559},{382,86329680991},{384,90881152117},{386,95001297565},{388,99963147805},{390,104453597992},{392,109837310021},{394,114722988623},{396,120585261143},{398,125873325588},{400,132247999328}};
using namespace chrono;
using namespace chrono_literals;

int main(int argc, char** argv){
    const size_t N_limit                = argc>1 ? strtol(argv[1],0,0) : 200;     // Argument 1: Number of vertices N
    const int generate_cpu_stats        = argc>2 ? strtol(argv[2],0,0) : 0;     // Argument 2: Boolean to generate CPU stats
    ofstream BONDS_CUDA_Wirz("Stats/IsomerspaceBondRMS_CUDA_Wirz" + to_string(N_limit) + ".txt");
    ofstream BONDS_CUDA_Buster("Stats/IsomerspaceBondRMS_CUDA_Buster" + to_string(N_limit) + ".txt");
    ofstream BONDS_CUDA_Flat("Stats/IsomerspaceBondRMS_CUDA_Flat" + to_string(N_limit) + ".txt");
    ofstream BONDS_Wirz("Stats/IsomerspaceBondRMS_Wirz" + to_string(N_limit) + ".txt");

    ofstream ANGLES_CUDA_Wirz("Stats/IsomerspaceAngleRMS_CUDA_Wirz" + to_string(N_limit) + ".txt");
    ofstream ANGLES_CUDA_Buster("Stats/IsomerspaceAngleRMS_CUDA_Buster" + to_string(N_limit) + ".txt");
    ofstream ANGLES_CUDA_Flat("Stats/IsomerspaceAngleRMS_CUDA_Flat" + to_string(N_limit) + ".txt");
    ofstream ANGLES_Wirz("Stats/IsomerspaceAngleRMS_Wirz" + to_string(N_limit) + ".txt");

    ofstream DIHEDRALS_CUDA_Wirz("Stats/IsomerspaceDihedralRMS_CUDA_Wirz" + to_string(N_limit) + ".txt");
    ofstream DIHEDRALS_CUDA_Buster("Stats/IsomerspaceDihedralRMS_CUDA_Buster" + to_string(N_limit) + ".txt");
    ofstream DIHEDRALS_CUDA_Flat("Stats/IsomerspaceDihedralRMS_CUDA_Flat" + to_string(N_limit) + ".txt");
    ofstream DIHEDRALS_Wirz("Stats/IsomerspaceDihedralRMS_Wirz" + to_string(N_limit) + ".txt");

    ofstream FLATNESS_CUDA_Wirz("Stats/IsomerspaceFlatness_CUDA_Wirz" + to_string(N_limit) + ".txt");
    ofstream FLATNESS_CUDA_Buster("Stats/IsomerspaceFlatness_CUDA_Buster" + to_string(N_limit) + ".txt");
    ofstream FLATNESS_CUDA_Flat("Stats/IsomerspaceFlatness_CUDA_Flat" + to_string(N_limit) + ".txt");
    ofstream FLATNESS_Wirz("Stats/IsomerspaceFlatness_Wirz" + to_string(N_limit) + ".txt");
    
    ofstream ENERGY_CUDA_Wirz("Stats/IsomerspaceEnergy_CUDA_Wirz" + to_string(N_limit) + ".txt");
    ofstream ENERGY_CUDA_Buster("Stats/IsomerspaceEnergy_CUDA_Buster" + to_string(N_limit) + ".txt");
    ofstream ENERGY_CUDA_Flat("Stats/IsomerspaceEnergy_CUDA_Flat" + to_string(N_limit) + ".txt");
    ofstream ENERGY_Wirz("Stats/IsomerspaceEnergy_Wirz" + to_string(N_limit) + ".txt");
    auto max_sample_size = 1000;


    for (size_t N = 20; N < N_limit+1; N+=2)
    {   
        if (N == 22) continue;
        BuckyGen::buckygen_queue Q = BuckyGen::start(N,false,false);  
        auto sample_size = min(max_sample_size,(int)num_fullerenes.find(N)->second);
        
        IsomerBatch batch0(N,sample_size,DEVICE_BUFFER);
        auto Nf = N/2 + 2;
        FullereneDual G;
        G.neighbours = neighbours_t(Nf, std::vector<node_t>(6));
        G.N = Nf;
        bool more_to_generate = true;

        auto path = "isomerspace_samples/dual_layout_" + to_string(N) + "_seed_42";
        ifstream isomer_sample(path,std::ios::binary);
        auto fsize = std::filesystem::file_size(path);
        std::vector<device_node_t> input_buffer(fsize/sizeof(device_node_t));
        auto available_samples = fsize / (Nf*6*sizeof(device_node_t));
        isomer_sample.read(reinterpret_cast<char*>(input_buffer.data()), Nf*6*sizeof(device_node_t)*available_samples);

        std::vector<int> random_IDs(available_samples);
        std::iota(random_IDs.begin(), random_IDs.end(), 0);
        std::shuffle(random_IDs.begin(), random_IDs.end(), std::mt19937{42});
        std::vector<int> id_subset(random_IDs.begin(), random_IDs.begin()+sample_size);

        using namespace cuda_io;
        IsomerQueue OptimisedQueue(N,0);
        IsomerQueue InputQueue(N,0);

        IsomerBatch GPUBatch(N,sample_size,DEVICE_BUFFER);
        for (int i = 0; i < sample_size; ++i){
                for (size_t j = 0; j < Nf; j++){
                    G.neighbours[j].clear();
                    for (size_t k = 0; k < 6; k++) {
                        auto u = input_buffer[id_subset[i]*Nf*6 + j*6 +k];
                        if(u != UINT16_MAX) G.neighbours[j].push_back(u);
                    }
                }
                InputQueue.insert(Graph(G),id_subset[i]);
                if(generate_cpu_stats > 0){
                    G.update();
                    PlanarGraph pG = G.dual_graph();
                    pG.layout2d  = pG.tutte_layout();
                    Polyhedron      P(pG);
                    P.points    = P.zero_order_geometry();
                    P.optimise();
                    OptimisedQueue.insert(P,id_subset[i]);
                }

        }


        InputQueue.refill_batch(GPUBatch);

        gpu_kernels::isomerspace_dual::dualise(GPUBatch);
        gpu_kernels::isomerspace_tutte::tutte_layout(GPUBatch, N*10);
        gpu_kernels::isomerspace_X0::zero_order_geometry(GPUBatch,4.0f);
        reset_convergence_statuses(GPUBatch);
        IsomerBatch WirzBatch(N,sample_size,DEVICE_BUFFER);
        IsomerBatch FlatBatch(N, sample_size, DEVICE_BUFFER);
        copy(WirzBatch,GPUBatch);
        copy(FlatBatch,GPUBatch);

        gpu_kernels::isomerspace_forcefield::optimise<PEDERSEN>(GPUBatch,N*5,N*5);
        gpu_kernels::isomerspace_forcefield::optimise<WIRZ>(WirzBatch,N*5,N*5);
        gpu_kernels::isomerspace_forcefield::optimise<FLATNESS_ENABLED>(FlatBatch,N*5,N*5);


        CuArray<float> RMSBonds_CUDA(max_sample_size);
        CuArray<float> RMSBonds_Wirz(max_sample_size);
        CuArray<float> RMSBonds_Flat(max_sample_size);

        CuArray<float> RMSAngles_CUDA(max_sample_size);
        CuArray<float> RMSAngles_Wirz(max_sample_size);
        CuArray<float> RMSAngles_Flat(max_sample_size);

        CuArray<float> RMSDihedrals_CUDA(max_sample_size);
        CuArray<float> RMSDihedrals_Wirz(max_sample_size);
        CuArray<float> RMSDihedrals_Flat(max_sample_size);

        CuArray<float> RMSFlatness_CUDA(max_sample_size);
        CuArray<float> RMSFlatness_Wirz(max_sample_size);
        CuArray<float> RMSFlatness_Flat(max_sample_size);

        CuArray<float> Energy_CUDA_WIRZ(max_sample_size);
        CuArray<float> Energy_CUDA_FLAT(max_sample_size);
        CuArray<float> Energy_CUDA_PEDERSEN(max_sample_size);
        CuArray<float> Energy_WIRZ(max_sample_size);


        gpu_kernels::isomerspace_forcefield::get_bond_rrmse<PEDERSEN>(GPUBatch,RMSBonds_CUDA);
        gpu_kernels::isomerspace_forcefield::get_angle_rrmse<PEDERSEN>(GPUBatch,RMSAngles_CUDA);
        gpu_kernels::isomerspace_forcefield::get_dihedral_rrmse<PEDERSEN>(GPUBatch,RMSDihedrals_CUDA);
        gpu_kernels::isomerspace_forcefield::get_flat_rmse<FLATNESS_ENABLED>(GPUBatch,RMSFlatness_CUDA);

        gpu_kernels::isomerspace_forcefield::get_bond_rrmse<WIRZ>(WirzBatch,RMSBonds_Wirz);
        gpu_kernels::isomerspace_forcefield::get_angle_rrmse<WIRZ>(WirzBatch,RMSAngles_Wirz);
        gpu_kernels::isomerspace_forcefield::get_dihedral_rrmse<WIRZ>(WirzBatch,RMSDihedrals_Wirz);
        gpu_kernels::isomerspace_forcefield::get_flat_rmse<FLATNESS_ENABLED>(WirzBatch,RMSFlatness_Wirz);

        gpu_kernels::isomerspace_forcefield::get_bond_rrmse<FLATNESS_ENABLED>(FlatBatch,RMSBonds_Flat);
        gpu_kernels::isomerspace_forcefield::get_angle_rrmse<FLATNESS_ENABLED>(FlatBatch,RMSAngles_Flat);
        gpu_kernels::isomerspace_forcefield::get_dihedral_rrmse<FLATNESS_ENABLED>(FlatBatch,RMSDihedrals_Flat);
        gpu_kernels::isomerspace_forcefield::get_flat_rmse<FLATNESS_ENABLED>(FlatBatch,RMSFlatness_Flat);

        gpu_kernels::isomerspace_forcefield::get_energies<WIRZ>(GPUBatch,Energy_CUDA_PEDERSEN);
        gpu_kernels::isomerspace_forcefield::get_energies<WIRZ>(WirzBatch,Energy_CUDA_WIRZ);
        gpu_kernels::isomerspace_forcefield::get_energies<WIRZ>(FlatBatch,Energy_CUDA_FLAT);


        BONDS_CUDA_Buster << N << ", " << sample_size << ", ";
        ANGLES_CUDA_Buster << N << ", " << sample_size << ", ";
        DIHEDRALS_CUDA_Buster << N << ", " << sample_size << ", ";
        FLATNESS_CUDA_Buster << N << ", " << sample_size << ", ";

        BONDS_CUDA_Wirz << N << ", " << sample_size << ", ";
        ANGLES_CUDA_Wirz << N << ", " << sample_size << ", ";
        DIHEDRALS_CUDA_Wirz << N << ", " << sample_size << ", ";
        FLATNESS_CUDA_Wirz << N << ", " << sample_size << ", ";

        BONDS_CUDA_Flat << N << ", " << sample_size << ", ";
        ANGLES_CUDA_Flat << N << ", " << sample_size << ", ";
        DIHEDRALS_CUDA_Flat << N << ", " << sample_size << ", ";
        FLATNESS_CUDA_Flat << N << ", " << sample_size << ", ";

        BONDS_Wirz << N << ", " << sample_size << ", ";
        ANGLES_Wirz << N << ", " << sample_size << ", ";
        DIHEDRALS_Wirz << N << ", " << sample_size << ", ";
        FLATNESS_Wirz << N << ", " << sample_size << ", ";

        ENERGY_CUDA_Buster << N << ", " << sample_size << ", ";
        ENERGY_CUDA_Wirz << N << ", " << sample_size << ", ";
        ENERGY_CUDA_Flat << N << ", " << sample_size << ", ";
        ENERGY_Wirz << N << ", " << sample_size << ", ";


        for (size_t j = 0; j < max_sample_size; j++){
            if (j != max_sample_size-1) {
                BONDS_CUDA_Buster << RMSBonds_CUDA[j] << ", ";
                ANGLES_CUDA_Buster << RMSAngles_CUDA[j] << ", ";
                DIHEDRALS_CUDA_Buster << RMSDihedrals_CUDA[j] << ", ";
                FLATNESS_CUDA_Buster << RMSFlatness_CUDA[j] << ", ";

                BONDS_CUDA_Wirz << RMSBonds_Wirz[j] << ", ";
                ANGLES_CUDA_Wirz << RMSAngles_Wirz[j] << ", ";
                DIHEDRALS_CUDA_Wirz << RMSDihedrals_Wirz[j] << ", ";
                FLATNESS_CUDA_Wirz << RMSFlatness_Wirz[j] << ", ";

                BONDS_CUDA_Flat << RMSBonds_Flat[j] << ", ";
                ANGLES_CUDA_Flat << RMSAngles_Flat[j] << ", ";
                DIHEDRALS_CUDA_Flat << RMSDihedrals_Flat[j] << ", ";
                FLATNESS_CUDA_Flat << RMSFlatness_Flat[j] << ", ";

                ENERGY_CUDA_Buster << Energy_CUDA_PEDERSEN[j] << ", ";
                ENERGY_CUDA_Wirz << Energy_CUDA_WIRZ[j] << ", ";
                ENERGY_CUDA_Flat << Energy_CUDA_FLAT[j] << ", ";
            }
            else{
                BONDS_CUDA_Buster << RMSBonds_CUDA[j] << "\n";
                ANGLES_CUDA_Buster << RMSAngles_CUDA[j] << "\n";
                DIHEDRALS_CUDA_Buster << RMSDihedrals_CUDA[j] << "\n";
                FLATNESS_CUDA_Buster << RMSFlatness_CUDA[j] << "\n";

                BONDS_CUDA_Wirz << RMSBonds_Wirz[j] << "\n";
                ANGLES_CUDA_Wirz << RMSAngles_Wirz[j] << "\n";
                DIHEDRALS_CUDA_Wirz << RMSDihedrals_Wirz[j] << "\n";
                FLATNESS_CUDA_Wirz << RMSFlatness_Wirz[j] << "\n";

                BONDS_CUDA_Flat << RMSBonds_Flat[j] << "\n";
                ANGLES_CUDA_Flat << RMSAngles_Flat[j] << "\n";
                DIHEDRALS_CUDA_Flat << RMSDihedrals_Flat[j] << "\n";
                FLATNESS_CUDA_Flat << RMSFlatness_Flat[j] << "\n";

                ENERGY_CUDA_Buster << Energy_CUDA_PEDERSEN[j] << "\n";
                ENERGY_CUDA_Wirz << Energy_CUDA_WIRZ[j] << "\n";
                ENERGY_CUDA_Flat << Energy_CUDA_FLAT[j] << "\n";
            }
        }

        if (generate_cpu_stats > 0){
            CuArray<float> RMSBonds_Fortran(max_sample_size);
            CuArray<float> RMSAngles_Fortran(max_sample_size);
            CuArray<float> RMSDihedrals_Fortran(max_sample_size);
            CuArray<float> RMSFlatness_Fortran(max_sample_size);
            CuArray<float> Energy_Fortran(max_sample_size);
            cuda_io::copy(OptimisedQueue.device_batch, OptimisedQueue.host_batch);
            gpu_kernels::isomerspace_forcefield::get_bond_rrmse<PEDERSEN>(OptimisedQueue.device_batch,RMSBonds_Fortran);
            gpu_kernels::isomerspace_forcefield::get_angle_rrmse<PEDERSEN>(OptimisedQueue.device_batch,RMSAngles_Fortran);
            gpu_kernels::isomerspace_forcefield::get_dihedral_rrmse<PEDERSEN>(OptimisedQueue.device_batch,RMSDihedrals_Fortran);
            gpu_kernels::isomerspace_forcefield::get_flat_rmse<FLATNESS_ENABLED>(OptimisedQueue.device_batch,RMSFlatness_Fortran);
            gpu_kernels::isomerspace_forcefield::get_energies<WIRZ>(OptimisedQueue.device_batch,Energy_Fortran);

            for (size_t j = 0; j < max_sample_size; j++){
            if (j != max_sample_size-1) {
                BONDS_Wirz << RMSBonds_Fortran[j] << ", ";
                ANGLES_Wirz << RMSAngles_Fortran[j] << ", ";
                DIHEDRALS_Wirz << RMSDihedrals_Fortran[j] << ", ";
                FLATNESS_Wirz << RMSFlatness_Fortran[j] << ", ";
                ENERGY_Wirz << Energy_Fortran[j] << ", ";
            }
            else{
                BONDS_Wirz << RMSBonds_Fortran[j] << "\n";
                ANGLES_Wirz << RMSAngles_Fortran[j] << "\n";
                DIHEDRALS_Wirz << RMSDihedrals_Fortran[j] << "\n";
                FLATNESS_Wirz << RMSFlatness_Fortran[j] << "\n";
                ENERGY_Wirz << Energy_Fortran[j] << "\n";
            } 
            }

        }
     }

}