#include "filesystem"
#include "numeric"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/gpu/benchmark_functions.hh"
#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/gpu/kernels.hh"
#include "fullerenes/gpu/cuda_io.hh"
#include "fullerenes/gpu/isomer_queue.hh"
const std::unordered_map<size_t,size_t> num_fullerenes = {{20,1},{22,0},{24,1},{26,1},{28,2},{30,3},{32,6},{34,6},{36,15},{38,17},{40,40},{42,45},{44,89},{46,116},{48,199},{50,271},{52,437},{54,580},{56,924},{58,1205},{60,1812},{62,2385},{64,3465},{66,4478},{68,6332},{70,8149},{72,11190},{74,14246},{76,19151},{78,24109},{80,31924},{82,39718},{84,51592},{86,63761},{88,81738},{90,99918},{92,126409},{94,153493},{96,191839},{98,231017},{100,285914},{102,341658},{104,419013},{106,497529},{108,604217},{110,713319},{112,860161},{114,1008444},{116,1207119},{118,1408553},{120,1674171},{122,1942929},{124,2295721},{126,2650866},{128,3114236},{130,3580637},{132,4182071},{134,4787715},{136,5566949},{138,6344698},{140,7341204},{142,8339033},{144,9604411},{146,10867631},{148,12469092},{150,14059174},{152,16066025},{154,18060979},{156,20558767},{158,23037594},{160,26142839},{162,29202543},{164,33022573},{166,36798433},{168,41478344},{170,46088157},{172,51809031},{174,57417264},{176,64353269},{178,71163452},{180,79538751},{182,87738311},{184,97841183},{186,107679717},{188,119761075},{190,131561744},{192,145976674},{194,159999462},{196,177175687},{198,193814658},{200,214127742},{202,233846463},{204,257815889},{206,281006325},{208,309273526},{210,336500830},{212,369580714},{214,401535955},{216,440216206},{218,477420176},{220,522599564},{222,565900181},{224,618309598},{226,668662698},{228,729414880},{230,787556069},{232,857934016},{234,925042498},{236,1006016526},{238,1083451816},{240,1176632247},{242,1265323971},{244,1372440782},{246,1474111053},{248,1596482232},{250,1712934069},{252,1852762875},{254,1985250572},{256,2144943655},{258,2295793276},{260,2477017558},{262,2648697036},{264,2854536850},{266,3048609900},{268,3282202941},{270,3501931260},{272,3765465341},{274,4014007928},{276,4311652376},{278,4591045471},{280,4926987377},{282,5241548270},{284,5618445787},{286,5972426835},{288,6395981131},{290,6791769082},{292,7267283603},{294,7710782991},{296,8241719706},{298,8738236515},{300,9332065811},{302,9884604767},{304,10548218751},{306,11164542762},{308,11902015724},{310,12588998862},{312,13410330482},{314,14171344797},{316,15085164571},{318,15930619304},{320,16942010457},{322,17880232383},{324,19002055537},{326,20037346408},{328,21280571390},{330,22426253115},{332,23796620378},{334,25063227406},{336,26577912084},{338,27970034826},{340,29642262229},{342,31177474996},{344,33014225318},{346,34705254287},{348,36728266430},{350,38580626759},{352,40806395661},{354,42842199753},{356,45278616586},{358,47513679057},{360,50189039868},{362,52628839448},{364,55562506886},{366,58236270451},{368,61437700788},{370,64363670678},{372,67868149215},{374,71052718441},{376,74884539987},{378,78364039771},{380,82532990559},{382,86329680991},{384,90881152117},{386,95001297565},{388,99963147805},{390,104453597992},{392,109837310021},{394,114722988623},{396,120585261143},{398,125873325588},{400,132247999328}};
using namespace gpu_kernels;
using namespace isomerspace_dual;
using namespace isomerspace_eigen;
using namespace isomerspace_forcefield;
using namespace isomerspace_hessian;
using namespace isomerspace_X0;
using namespace isomerspace_tutte;
int main(int argc, char** argv){
    const size_t N                = argc>1 ? strtol(argv[1],0,0) : 60;     // Argument 1: Number of vertices 
    //const size_t Mlanczos_steps   = argc>2 ? strtol(argv[2],0,0) : 40;     // Argument 2: Number of Lanczos steps
    float reldelta                = argc>2 ? strtof(argv[2],0) : 1e-5;    // Argument 3: Relative delta
    int isomer_num                = argc>3 ? strtol(argv[3],0,0) : 0;     // Argument 4: Isomer number
    std::string spiral_           = argc>4 ? argv[4] : "C60-[1,7,9,11,13,15,18,20,22,24,26,32]-fullerene";          // Argument 2: Spiral
    std::string name_             = argc>5 ? argv[5] : "C60ih";          // Argument 2: Spiral
    std::string filename          = argc>6 ? argv[6] : "hessian_validation";        // Argument 2: Filename

    std::ofstream hess_analytical(filename + "_analytical.csv");
    std::ofstream hess_numerical(filename + "_numerical.csv");
    std::ofstream hess_cols(filename + "_cols.csv");
    std::ofstream cubic_graph;//("C" + to_string(N) + "_CubicGraph_" + to_string(isomer_num) + ".bin");
    std::ofstream geometry;//("C" + to_string(N) + "_Geometry_" + to_string(isomer_num) + ".bin");
    if (isomer_num == -1){
        cubic_graph = std::ofstream(name_ + "_CubicGraph.bin");
        geometry = std::ofstream(name_ + "_Geometry.bin");
    }

    bool more_to_do = true;
    auto n_isomers = num_fullerenes.find(N)->second;
    Graph G;
    auto Nf = N/2 +2;
    G.neighbours = neighbours_t(Nf, std::vector<node_t>(6));
    G.neighbours.resize(Nf);
    G.N = Nf;
    PlanarGraph Pg;
    auto batch_size = min(isomerspace_forcefield::optimal_batch_size(N)*16, (int)n_isomers);
    //BuckyGen::buckygen_queue BuckyQ = BuckyGen::start(N,0,0);  
    IsomerBatch Bhost(N,batch_size,HOST_BUFFER);
    int Nd = LaunchCtx::get_device_count();
    IsomerBatch Bdev(N,batch_size,DEVICE_BUFFER,0);
    if (isomer_num == -1) {
        spiral_nomenclature C60name(spiral_);    
        Triangulation C60dual(C60name);
        auto C60cubic = C60dual.dual_graph();
        for(int i = 0;  i < batch_size; i++) Bhost.append(C60cubic,0, false);
    } else {
        BuckyGen::buckygen_queue BuckyQ = BuckyGen::start(N, false, false);  
        for (size_t i = 0; i < batch_size; i++)
        {   
            int C90_problem_isomers[3] = {66530,67259,83863}; 
            bool success = BuckyGen::next_fullerene(BuckyQ, G);
            if(success) Bhost.append(G,i);
            else break;
        }
    }
    LaunchCtx ctx(0);
    LaunchPolicy policy = LaunchPolicy::SYNC;
    cuda_io::copy(Bdev, Bhost);

    dualise(Bdev);
    tutte_layout(Bdev, (int)10*N);
    zero_order_geometry(Bdev, 4.0);
    optimise<PEDERSEN>(Bdev, 5*N, 5*N);
    CuArray<device_real_t> hessians(N*3*3*10 * batch_size);
    CuArray<device_node_t> cols(N*3*3*10 * batch_size);
    CuArray<device_real_t> lambda_mins(batch_size);
    CuArray<device_real_t> lambda_maxs(batch_size);

    auto T0 = std::chrono::steady_clock::now();
    compute_hessians<PEDERSEN>(Bdev, hessians, cols);
    //std::cout << hessians << std::endl;
    auto T1 = std::chrono::steady_clock::now();
    std::cout << "Hessian time: " << std::chrono::duration_cast<std::chrono::microseconds>(T1-T0).count()/ (float)batch_size << " us/ isomer" << std::endl;
    
    hess_analytical.write((char*)hessians.data, N*3*3*10* batch_size*sizeof(device_real_t));
    hess_cols.write((char*)cols.data, N*3*3*10* batch_size*sizeof(device_node_t));


    //compute_hessians_fd<PEDERSEN>(Bdev, hessians, cols, reldelta);
    //hess_numerical.write((char*)hessians.data, N*3*3*10* batch_size*sizeof(device_real_t));
    
    //std::cout << "Eigenvalues: " << eigenvalues << std::endl;
    //cuda_io::copy(Bhost, Bdev); // Copy back to host
    //auto FG = Bhost.get_isomer(0).value();
    //cubic_graph.write((char*)Bhost.cubic_neighbours, N*3*sizeof(device_node_t));
    //geometry.write((char*)Bhost.X, N*3*sizeof(device_real_t));
    
    //eigensolve_cusolver(Bdev, hessians, cols, eigenvalues);
    //eigensolve_cusolver(Bdev, hessians, cols, eigenvalues);
    //std::cout << "cuSOLVE Eigensolver took: " << (time/1us)/(float)batch_size << " us / graph" << std::endl;
    CuArray<device_real_t> Q(N*3*N*3*batch_size);
    CuArray<device_real_t> eigs(N*3*batch_size);
    //lambda_max(Bdev, hessians, cols, lambda_maxs);
    auto start = std::chrono::steady_clock::now();
    spectrum_ends(Bdev, hessians, cols, lambda_mins, lambda_maxs, 50);
    auto time = std::chrono::steady_clock::now() - start;
    std::cout << "Spectrum Ends took: " << (time/1us)/(float)batch_size << " us / graph" << std::endl;

    start = std::chrono::steady_clock::now();
    eigensolve(Bdev, Q, hessians, cols, eigs);
    ofstream eigs_out("eigs.float32", ios::binary); eigs_out.write((char*)eigs.data, N*3*batch_size*sizeof(device_real_t)); eigs_out.close();
    time = std::chrono::steady_clock::now() - start;
    std::cout << "Full eigensolve took: " << (time/1us)/(float)batch_size << " us / graph" << std::endl;

    std::vector<device_real_t> min_eigs_ref(batch_size), max_eigs_ref(batch_size);
    
    for(int i = 0; i < batch_size; i++) {
        std::sort(eigs.data + i*N*3, eigs.data + (i+1)*N*3);
        min_eigs_ref[i] = eigs[i*N*3 + 6];
        max_eigs_ref[i] = eigs[(i+1)*N*3 - 1];
    }
    std::vector<device_real_t> rel_err_min(batch_size), abs_err_min(batch_size), rel_err_max(batch_size), abs_err_max(batch_size), rel_err_width(batch_size), abs_err_width(batch_size);
    std::vector<int> nan_or_inf;
    for(int i = 0; i < batch_size; i++) {
        if (std::isnan(lambda_mins[i]) || std::isinf(lambda_mins[i]) || std::isnan(lambda_maxs[i]) || std::isinf(lambda_maxs[i]) || std::isnan(min_eigs_ref[i]) || std::isinf(min_eigs_ref[i]) || std::isnan(max_eigs_ref[i]) || std::isinf(max_eigs_ref[i])) {
            nan_or_inf.push_back(i);
            continue;
        }
        device_real_t epsilon = numeric_limits<device_real_t>::epsilon()*1e1;
        rel_err_min[i] = fabs(lambda_mins[i] - min_eigs_ref[i] + epsilon)/fabs(min_eigs_ref[i] + epsilon);
        abs_err_min[i] = fabs(lambda_mins[i] - min_eigs_ref[i] + epsilon);
        rel_err_max[i] = fabs(lambda_maxs[i] - max_eigs_ref[i] + epsilon)/fabs(max_eigs_ref[i] + epsilon);
        abs_err_max[i] = fabs(lambda_maxs[i] - max_eigs_ref[i] + epsilon);
        rel_err_width[i] = fabs( (lambda_maxs[i] - lambda_mins[i]) - (max_eigs_ref[i] - min_eigs_ref[i]) + epsilon)/fabs(max_eigs_ref[i] - min_eigs_ref[i] + epsilon);
        abs_err_width[i] = fabs( (lambda_maxs[i] - lambda_mins[i]) - (max_eigs_ref[i] - min_eigs_ref[i]) + epsilon);
    }

    ofstream file("MinLambdaError_Relative.float32", ios::out | ios::binary);
    file.write((char*)rel_err_min.data(), batch_size*sizeof(device_real_t));
    file.close();
    ofstream file2("MinLambdaError_Absolute.float32", ios::out | ios::binary);
    file2.write((char*)abs_err_min.data(), batch_size*sizeof(device_real_t));
    file2.close();
    ofstream file3("MaxLambdaError_Relative.float32", ios::out | ios::binary);
    file3.write((char*)rel_err_max.data(), batch_size*sizeof(device_real_t));
    file3.close();
    ofstream file4("MaxLambdaError_Absolute.float32", ios::out | ios::binary);
    file4.write((char*)abs_err_max.data(), batch_size*sizeof(device_real_t));
    file4.close();
    ofstream file5("WidthError_Relative.float32", ios::out | ios::binary);
    file5.write((char*)rel_err_width.data(), batch_size*sizeof(device_real_t));
    file5.close();
    ofstream file6("WidthError_Absolute.float32", ios::out | ios::binary);
    file6.write((char*)abs_err_width.data(), batch_size*sizeof(device_real_t));
    file6.close();

    ofstream file7("Q.float32", ios::out | ios::binary);
    file7.write((char*)Q.data, Q.size()*sizeof(device_real_t));
    file7.close();

    //std::cout << "Reference min eigenvalues: " << min_eigs_ref << std::endl;
    //std::cout << lambda_maxs << std::endl;
    
        //eigensolve(Bdev, Q, hessians, cols, eigenvalues);
        //start = std::chrono::steady_clock::now();
        //spectrum_ends(Bdev, hessians, cols, lambda_min, lambda_max);
        //time = std::chrono::steady_clock::now() - start;
        

        // Find the 20 smallest elements in lambda_maxs
        /* std::vector<device_real_t> lambda_maxs_vec(lambda_maxs.data, lambda_maxs.data + batch_size);
        std::vector<device_real_t> smallest_40(40);
        std::vector<size_t> smallest_40_indices(40);
        std::vector<size_t> indices(lambda_maxs_vec.size());
        std::iota(indices.begin(), indices.end(), 0);
        std::sort(indices.begin(), indices.end(), [&lambda_maxs_vec](size_t i1, size_t i2) {return lambda_maxs_vec[i1] < lambda_maxs_vec[i2];});
        for (int i = 0; i < 40; i++) {
            smallest_40[i] = lambda_maxs_vec[indices[i]];
            smallest_40_indices[i] = indices[i];
        }
        std::cout << "40 smallest lambda_maxs: " << smallest_40 << std::endl;
        std::cout << "40 smallest indices: " << smallest_40_indices << std::endl;
        std::sort(indices.begin(), indices.end(), [&lambda_maxs_vec](size_t i1, size_t i2) {return lambda_maxs_vec[i1] > lambda_maxs_vec[i2];});
        std::vector<device_real_t> largest_40(40);
        std::vector<size_t> largest_40_indices(40);
        for (int i = 0; i < 40; i++) {
            largest_40[i] = lambda_maxs_vec[indices[i]];
            largest_40_indices[i] = indices[i];
        }
        std::cout << "40 largest lambda_maxs: " << largest_40 << std::endl;
        std::cout << "40 largest indices: " << largest_40_indices << std::endl;
        std::ofstream output_tridiags("tridiags.bin", std::ios::binary);
        std::ofstream output_Q("Q.bin", std::ios::binary); */
        //output_tridiags.write((char*)eigenvalues.data, N*3*2* batch_size*sizeof(device_real_t));
        //output_Q.write((char*)Q.data, N*3*N*3* batch_size*sizeof(device_real_t));

        for(int i = N*2; i < N*3; i++) {
            //std::cout << "Q[" << i  << "] : " << std::vector<device_real_t>(Q.data + i*N*3, Q.data + (i+1)*N*3) << std::endl;
        }
        //std::cout << "Eigenvalues: " << eigenvalues << std::endl;
        //FILE* f = fopen((name_ + "_Geometry.mol2").c_str(), "w");
        //FG.to_mol2(FG,f);

    //    std::cout << Bdev << std::endl;
    }
