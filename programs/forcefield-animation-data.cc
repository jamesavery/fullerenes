#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/polyhedron.hh"
#include <chrono>
#include <fstream>

const std::unordered_map<size_t,size_t> num_fullerenes = {{20,1},{22,0},{24,1},{26,1},{28,2},{30,3},{32,6},{34,6},{36,15},{38,17},{40,40},{42,45},{44,89},{46,116},{48,199},{50,271},{52,437},{54,580},{56,924},{58,1205},{60,1812},{62,2385},{64,3465},{66,4478},{68,6332},{70,8149},{72,11190},{74,14246},{76,19151},{78,24109},{80,31924},{82,39718},{84,51592},{86,63761},{88,81738},{90,99918},{92,126409},{94,153493},{96,191839},{98,231017},{100,285914},{102,341658},{104,419013},{106,497529},{108,604217},{110,713319},{112,860161},{114,1008444},{116,1207119},{118,1408553},{120,1674171},{122,1942929},{124,2295721},{126,2650866},{128,3114236},{130,3580637},{132,4182071},{134,4787715},{136,5566949},{138,6344698},{140,7341204},{142,8339033},{144,9604411},{146,10867631},{148,12469092},{150,14059174},{152,16066025},{154,18060979},{156,20558767},{158,23037594},{160,26142839},{162,29202543},{164,33022573},{166,36798433},{168,41478344},{170,46088157},{172,51809031},{174,57417264},{176,64353269},{178,71163452},{180,79538751},{182,87738311},{184,97841183},{186,107679717},{188,119761075},{190,131561744},{192,145976674},{194,159999462},{196,177175687},{198,193814658},{200,214127742},{202,233846463},{204,257815889},{206,281006325},{208,309273526},{210,336500830},{212,369580714},{214,401535955},{216,440216206},{218,477420176},{220,522599564},{222,565900181},{224,618309598},{226,668662698},{228,729414880},{230,787556069},{232,857934016},{234,925042498},{236,1006016526},{238,1083451816},{240,1176632247},{242,1265323971},{244,1372440782},{246,1474111053},{248,1596482232},{250,1712934069},{252,1852762875},{254,1985250572},{256,2144943655},{258,2295793276},{260,2477017558},{262,2648697036},{264,2854536850},{266,3048609900},{268,3282202941},{270,3501931260},{272,3765465341},{274,4014007928},{276,4311652376},{278,4591045471},{280,4926987377},{282,5241548270},{284,5618445787},{286,5972426835},{288,6395981131},{290,6791769082},{292,7267283603},{294,7710782991},{296,8241719706},{298,8738236515},{300,9332065811},{302,9884604767},{304,10548218751},{306,11164542762},{308,11902015724},{310,12588998862},{312,13410330482},{314,14171344797},{316,15085164571},{318,15930619304},{320,16942010457},{322,17880232383},{324,19002055537},{326,20037346408},{328,21280571390},{330,22426253115},{332,23796620378},{334,25063227406},{336,26577912084},{338,27970034826},{340,29642262229},{342,31177474996},{344,33014225318},{346,34705254287},{348,36728266430},{350,38580626759},{352,40806395661},{354,42842199753},{356,45278616586},{358,47513679057},{360,50189039868},{362,52628839448},{364,55562506886},{366,58236270451},{368,61437700788},{370,64363670678},{372,67868149215},{374,71052718441},{376,74884539987},{378,78364039771},{380,82532990559},{382,86329680991},{384,90881152117},{386,95001297565},{388,99963147805},{390,104453597992},{392,109837310021},{394,114722988623},{396,120585261143},{398,125873325588},{400,132247999328}};
using namespace chrono;
using namespace chrono_literals;

#include "fullerenes/isomer_queue.hh"
#include "fullerenes/device_io.hh"
#include "fullerenes/gpu/kernels.hh"
#include "fullerenes/gpu/benchmark_functions.hh"
#include "fullerenes/progress_bar.hh"
#include "numeric"
#include "random"
using namespace gpu_kernels;

int main(int argc, char** argv){

    size_t N                = argc > 1 ? strtol(argv[1],0,0) : (size_t)20;        // Argument 1: Start of range of N
    size_t isomers          = argc > 2 ? strtol(argv[2],0,0) : (size_t)num_fullerenes.find(N)->second;    // Argument 2: Number of isomers to generate
    Graph G;
    auto Nf = N/2 + 2;
    G.neighbours = neighbours_t(Nf, std::vector<node_t>(6));
    G.N = Nf;
    ProgressBar progress_bar('=', 70);

    auto batch_size = isomerspace_forcefield::optimal_batch_size(N);
    auto n_fullerenes = (int)num_fullerenes.find(N)->second;
    auto sample_size = min((int)isomers,n_fullerenes);
    
    bool more_to_generate = true;
    ofstream out0; out0.open("AnimationCoordinates.bin", std::ios::binary);
    ofstream out1; out1.open("AnimationGraphs.bin", std::ios::binary);
    ofstream out2; out2.open("Animation2DCoordinates.bin", std::ios::binary); //Tutte layouts
    
    auto path = "isomerspace_samples/dual_layout_" + to_string(N) + "_seed_42";
    ifstream isomer_sample(path,std::ios::binary);
    auto fsize = file_size(path);
    std::vector<device_node_t> input_buffer(fsize/sizeof(device_node_t));
    auto available_samples = fsize / (Nf*6*sizeof(device_node_t));
    isomer_sample.read(reinterpret_cast<char*>(input_buffer.data()), Nf*6*sizeof(device_node_t)*available_samples);

    std::vector<int> random_IDs(available_samples);
    std::iota(random_IDs.begin(), random_IDs.end(), 0);
    std::shuffle(random_IDs.begin(), random_IDs.end(), std::mt19937{42});
    std::vector<int> id_subset(random_IDs.begin(), random_IDs.begin()+sample_size);

    
    auto finished_isomers = 0;
    IsomerBatch<GPU> batch0(N,sample_size);
    IsomerBatch<GPU> batch1(N,sample_size);
    IsomerBatch<CPU> h_batch(N,sample_size);
    device_io::IsomerQueue Q0(N);
    Q0.resize(min(n_fullerenes,sample_size));
    for (int i = 0; i < sample_size; i++){
        for (size_t j = 0; j < Nf; j++){
            G.neighbours[j].clear();
            for (size_t k = 0; k < 6; k++) {
                auto u = input_buffer[id_subset[i]*Nf*6 + j*6 +k];
                if(u != UINT16_MAX) G.neighbours[j].push_back(u);
            }
        }
        Q0.insert(G, i);
    }

    Q0.refill_batch(batch0);
    isomerspace_dual::dualize(batch0);


    std::cout << std::endl;
    for (size_t i = 0; i < 500; i++)
    {
        isomerspace_tutte::tutte_layout(batch0, i*N/20);
        device_io::copy(h_batch, batch0);
        out2.write(reinterpret_cast<char*>(h_batch.xys), N*2*sizeof(device_real_t)*sample_size);
        out2.flush();
        progress_bar.update_progress((i+1)/500.f);
    }

    isomerspace_X0::zero_order_geometry(batch0, 4.0);

    for (int i = 0; i < 100; i++){
        isomerspace_forcefield::optimize<PEDERSEN>(batch0,int(N*5/100),N*5);
        device_io::copy(h_batch, batch0);
        out0.write(reinterpret_cast<char*>(h_batch.X), N*3*sizeof(device_real_t)*sample_size); 
        out0.flush();
        progress_bar.update_progress((i+1)/100.f);
    }

    CuArray<device_real_t> Ecc(sample_size);
    isomerspace_properties::eccentricities(batch0, Ecc);
    out1.write(reinterpret_cast<char*>(h_batch.cubic_neighbours), N*3*sizeof(device_node_t)*sample_size);
    vector<int> fail_ids;
    for (int i = 0; i < sample_size; i++){
        if (h_batch.statuses[i] == IsomerStatus::FAILED){
            fail_ids.push_back(i);
        }
    }
    std::vector<int> ecc_ids(sample_size);
    std::iota(ecc_ids.begin(), ecc_ids.end(), 0);
    std::sort(ecc_ids.begin(), ecc_ids.end(), [&Ecc](int i1, int i2) {return Ecc[i1] < Ecc[i2];});

    std::cout << "\nEccentricities ids: " << ecc_ids << std::endl;


    std::cout <<  "\nFailed isomers: " << fail_ids << std::endl;

    out0.close();
    out1.close();
    out2.close();
}
