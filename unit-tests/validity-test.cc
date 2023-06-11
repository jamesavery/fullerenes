#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/gpu/isomer_queue.hh"
#include "fullerenes/gpu/isomer_batch.hh"
#include "fullerenes/gpu/kernels.hh"
#include "fullerenes/gpu/cuda_io.hh"
const std::unordered_map<size_t,size_t> num_fullerenes = {{20,1},{22,0},{24,1},{26,1},{28,2},{30,3},{32,6},{34,6},{36,15},{38,17},{40,40},{42,45},{44,89},{46,116},{48,199},{50,271},{52,437},{54,580},{56,924},{58,1205},{60,1812},{62,2385},{64,3465},{66,4478},{68,6332},{70,8149},{72,11190},{74,14246},{76,19151},{78,24109},{80,31924},{82,39718},{84,51592},{86,63761},{88,81738},{90,99918},{92,126409},{94,153493},{96,191839},{98,231017},{100,285914},{102,341658},{104,419013},{106,497529},{108,604217},{110,713319},{112,860161},{114,1008444},{116,1207119},{118,1408553},{120,1674171},{122,1942929},{124,2295721},{126,2650866},{128,3114236},{130,3580637},{132,4182071},{134,4787715},{136,5566949},{138,6344698},{140,7341204},{142,8339033},{144,9604411},{146,10867631},{148,12469092},{150,14059174},{152,16066025},{154,18060979},{156,20558767},{158,23037594},{160,26142839},{162,29202543},{164,33022573},{166,36798433},{168,41478344},{170,46088157},{172,51809031},{174,57417264},{176,64353269},{178,71163452},{180,79538751},{182,87738311},{184,97841183},{186,107679717},{188,119761075},{190,131561744},{192,145976674},{194,159999462},{196,177175687},{198,193814658},{200,214127742},{202,233846463},{204,257815889},{206,281006325},{208,309273526},{210,336500830},{212,369580714},{214,401535955},{216,440216206},{218,477420176},{220,522599564},{222,565900181},{224,618309598},{226,668662698},{228,729414880},{230,787556069},{232,857934016},{234,925042498},{236,1006016526},{238,1083451816},{240,1176632247},{242,1265323971},{244,1372440782},{246,1474111053},{248,1596482232},{250,1712934069},{252,1852762875},{254,1985250572},{256,2144943655},{258,2295793276},{260,2477017558},{262,2648697036},{264,2854536850},{266,3048609900},{268,3282202941},{270,3501931260},{272,3765465341},{274,4014007928},{276,4311652376},{278,4591045471},{280,4926987377},{282,5241548270},{284,5618445787},{286,5972426835},{288,6395981131},{290,6791769082},{292,7267283603},{294,7710782991},{296,8241719706},{298,8738236515},{300,9332065811},{302,9884604767},{304,10548218751},{306,11164542762},{308,11902015724},{310,12588998862},{312,13410330482},{314,14171344797},{316,15085164571},{318,15930619304},{320,16942010457},{322,17880232383},{324,19002055537},{326,20037346408},{328,21280571390},{330,22426253115},{332,23796620378},{334,25063227406},{336,26577912084},{338,27970034826},{340,29642262229},{342,31177474996},{344,33014225318},{346,34705254287},{348,36728266430},{350,38580626759},{352,40806395661},{354,42842199753},{356,45278616586},{358,47513679057},{360,50189039868},{362,52628839448},{364,55562506886},{366,58236270451},{368,61437700788},{370,64363670678},{372,67868149215},{374,71052718441},{376,74884539987},{378,78364039771},{380,82532990559},{382,86329680991},{384,90881152117},{386,95001297565},{388,99963147805},{390,104453597992},{392,109837310021},{394,114722988623},{396,120585261143},{398,125873325588},{400,132247999328}};

template <typename T>
std::pair<std::string,bool> validate_data(const std::string& validation_path, const std::vector<T>& array, const int num_isomers, const int elements_per_isomer){
    std::vector<T> validation(num_isomers * elements_per_isomer * sizeof(T));
    ifstream file(validation_path, ios::binary);
    file.read((char*)validation.data(), num_isomers * elements_per_isomer * sizeof(T));
    std::string fail_string = "Failed in Isomers: [";
    bool passed = true;

    for (size_t i = 0; i < num_isomers; i++){   
        auto isomer_passed = true;
        for (size_t j = 0; j < elements_per_isomer; j++){
            isomer_passed &= array[i*elements_per_isomer + j ] == validation[i*elements_per_isomer + j] || isnan(array[i*elements_per_isomer + j ]) && isnan(validation[i*elements_per_isomer + j]);
        }
        if(!isomer_passed) fail_string += to_string(i) + ", ";
        passed &= isomer_passed;
    }

    if(passed) {return {"Passed Succesfully!", true};}
    else {return {fail_string + "]", false};}
}


int main(int ac, char **argv){
    
    int N                = strtol(argv[1],0,0);     // Argument 1: Number of vertices N
    auto n_isomers = num_fullerenes.find(N)->second;

    IsomerBatch d_batch(N, n_isomers, DEVICE_BUFFER);
    IsomerBatch h_batch(N, n_isomers, HOST_BUFFER);
    IsomerBatch h_x0_batch(N, n_isomers, HOST_BUFFER);

    BuckyGen::buckygen_queue Q = BuckyGen::start(N,false, false);
    FullereneDual FD;
    bool more_to_generate = true;
    queue<std::tuple<Polyhedron, size_t, IsomerStatus>> out_queue;
    queue<std::tuple<Polyhedron, size_t, IsomerStatus>> x0_out_queue;

    map<int, Polyhedron> out_map;
    map<int, Polyhedron> x0_map;

    cuda_io::IsomerQueue BQ(N);
    int I = 0;


    while(more_to_generate)
    {
        more_to_generate &= BuckyGen::next_fullerene(Q,FD);
        if(!more_to_generate)break;
        FD.update();
        PlanarGraph G = FD.dual_graph();
        BQ.insert(G,I, LaunchCtx(), LaunchPolicy::SYNC, false);
        ++I;

    }
    BQ.refill_batch(d_batch);
    using namespace gpu_kernels;
    isomerspace_tutte::tutte_layout(d_batch);
    isomerspace_X0::zero_order_geometry(d_batch, 4.0f);
    cuda_io::copy(h_x0_batch, d_batch);
    cuda_io::reset_convergence_statuses(d_batch);
    
    isomerspace_forcefield::optimise<PEDERSEN>(d_batch, N*5, N*5);
    cuda_io::copy(h_batch,d_batch);
    cuda_io::output_to_queue(out_queue,h_batch,true);
    cuda_io::output_to_queue(x0_out_queue, h_x0_batch, true);
    while (!out_queue.empty()){   
        auto& [P, ID, status] = out_queue.front();
        auto& [Px0, IDx0, statusx0] = x0_out_queue.front();
        out_map.insert({ID, P});
        x0_map.insert({IDx0, Px0});
        out_queue.pop();
        x0_out_queue.pop();

    } 

    std::vector<int> neighbours(n_isomers*N*3);
    std::vector<float> xys(n_isomers*N*2);
    std::vector<float> X0(n_isomers*N*3);
    std::vector<float> X(n_isomers*N*3);

    for (size_t i = 0; i < n_isomers; i++)
    for (size_t j = 0; j < N; j++){
        xys[i*N*2 + j*2 + 0] = out_map.at(i).layout2d[j].first;
        xys[i*N*2 + j*2 + 1] = out_map.at(i).layout2d[j].second;
        for (size_t k = 0; k < 3; k++){
            neighbours[i*N*3 + j*3 + k] = out_map.at(i).neighbours[j][k];
            X0[i*N*3 + j*3 + k] = x0_map.at(i).points[j][k];
            X[i*N*3 + j*3 + k]  = out_map.at(i).points[j][k];
    }}
    
    auto [graph_string, graph_passed]   = validate_data("isomerspace_"+ to_string(N) +"_neighbours.bin", neighbours, n_isomers, N*3);
    auto [tutte_string, tutte_passed]   = validate_data("isomerspace_"+ to_string(N) +"_xys.bin", xys, n_isomers, N*2);
    auto [X0_string, X0_passed]         = validate_data("isomerspace_"+ to_string(N) +"_X0.bin", X0, n_isomers, N*3);
    auto [Xopt_string, Xopt_passed]     = validate_data("isomerspace_"+ to_string(N) +"_X_opt.bin", X, n_isomers, N*3);

    std::cout << "Graph Match " << graph_string << std::endl;
    std::cout << "Tutte Match " << tutte_string << std::endl;
    std::cout << "X0 Match " << X0_string << std::endl;
    std::cout << "Xopt Match " << Xopt_string << std::endl;
}