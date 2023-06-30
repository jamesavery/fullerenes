#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/gpu/isomer_queue.hh"
#include "fullerenes/gpu/isomer_batch.hh"
#include "fullerenes/gpu/kernels.hh"
#include "fullerenes/gpu/cuda_io.hh"
#include "fullerenes/gpu/cu_array.hh"
#include <numeric>
#include <future>
#include <random>
const std::unordered_map<size_t,size_t> num_fullerenes = {{20,1},{22,0},{24,1},{26,1},{28,2},{30,3},{32,6},{34,6},{36,15},{38,17},{40,40},{42,45},{44,89},{46,116},{48,199},{50,271},{52,437},{54,580},{56,924},{58,1205},{60,1812},{62,2385},{64,3465},{66,4478},{68,6332},{70,8149},{72,11190},{74,14246},{76,19151},{78,24109},{80,31924},{82,39718},{84,51592},{86,63761},{88,81738},{90,99918},{92,126409},{94,153493},{96,191839},{98,231017},{100,285914},{102,341658},{104,419013},{106,497529},{108,604217},{110,713319},{112,860161},{114,1008444},{116,1207119},{118,1408553},{120,1674171},{122,1942929},{124,2295721},{126,2650866},{128,3114236},{130,3580637},{132,4182071},{134,4787715},{136,5566949},{138,6344698},{140,7341204},{142,8339033},{144,9604411},{146,10867631},{148,12469092},{150,14059174},{152,16066025},{154,18060979},{156,20558767},{158,23037594},{160,26142839},{162,29202543},{164,33022573},{166,36798433},{168,41478344},{170,46088157},{172,51809031},{174,57417264},{176,64353269},{178,71163452},{180,79538751},{182,87738311},{184,97841183},{186,107679717},{188,119761075},{190,131561744},{192,145976674},{194,159999462},{196,177175687},{198,193814658},{200,214127742},{202,233846463},{204,257815889},{206,281006325},{208,309273526},{210,336500830},{212,369580714},{214,401535955},{216,440216206},{218,477420176},{220,522599564},{222,565900181},{224,618309598},{226,668662698},{228,729414880},{230,787556069},{232,857934016},{234,925042498},{236,1006016526},{238,1083451816},{240,1176632247},{242,1265323971},{244,1372440782},{246,1474111053},{248,1596482232},{250,1712934069},{252,1852762875},{254,1985250572},{256,2144943655},{258,2295793276},{260,2477017558},{262,2648697036},{264,2854536850},{266,3048609900},{268,3282202941},{270,3501931260},{272,3765465341},{274,4014007928},{276,4311652376},{278,4591045471},{280,4926987377},{282,5241548270},{284,5618445787},{286,5972426835},{288,6395981131},{290,6791769082},{292,7267283603},{294,7710782991},{296,8241719706},{298,8738236515},{300,9332065811},{302,9884604767},{304,10548218751},{306,11164542762},{308,11902015724},{310,12588998862},{312,13410330482},{314,14171344797},{316,15085164571},{318,15930619304},{320,16942010457},{322,17880232383},{324,19002055537},{326,20037346408},{328,21280571390},{330,22426253115},{332,23796620378},{334,25063227406},{336,26577912084},{338,27970034826},{340,29642262229},{342,31177474996},{344,33014225318},{346,34705254287},{348,36728266430},{350,38580626759},{352,40806395661},{354,42842199753},{356,45278616586},{358,47513679057},{360,50189039868},{362,52628839448},{364,55562506886},{366,58236270451},{368,61437700788},{370,64363670678},{372,67868149215},{374,71052718441},{376,74884539987},{378,78364039771},{380,82532990559},{382,86329680991},{384,90881152117},{386,95001297565},{388,99963147805},{390,104453597992},{392,109837310021},{394,114722988623},{396,120585261143},{398,125873325588},{400,132247999328}};
using namespace gpu_kernels;
using namespace cuda_io;
#define SYNC LaunchPolicy::SYNC
#define ASYNC LaunchPolicy::ASYNC

int main(int ac, char **argv){
    int N_start                 = ac > 1 ? strtol(argv[1],0,0) : 20;     // Argument 1: Number of vertices N
    int N_end                   = ac > 2 ? strtol(argv[2],0,0) : 200;     // Argument 1: Number of vertices N
    ofstream file_bond_adiff("async_validation/adiff_bond.txt");
    ofstream file_angle_adiff("async_validation/adiff_angle.txt");
    ofstream file_dihedral_adiff("async_validation/adiff_dihedral.txt");
    ofstream file_bond_rdiff("async_validation/rdiff_bond.txt");
    ofstream file_angle_rdiff("async_validation/rdiff_angle.txt");
    ofstream file_dihedral_rdiff("async_validation/rdiff_dihedral.txt");
    ofstream file_iterations("async_validation/iterations.txt");

    for (size_t N = N_start; N < N_end + 1; N+=2)
    {
    
    
    if(N==22) continue;
    auto n_isomers = num_fullerenes.find(N)->second;
    Graph G;
    auto Nf = N/2 +2;
    G.neighbours = neighbours_t(Nf, std::vector<node_t>(6));
    G.N = Nf;
    auto bucky = BuckyGen::start(N, false, false);

    std::string path = "isomerspace_samples/dual_layout_" + to_string(N) + "_seed_42";
    auto fsize = file_size(path);
    auto n_samples = fsize / (Nf * 6 * sizeof(device_node_t));
    ifstream in_file(path,std::ios::binary);
    std::vector<device_node_t> dual_neighbours(n_samples * Nf * 6);
    in_file.read((char*)dual_neighbours.data(),n_samples * Nf * 6 * sizeof(device_node_t));
    
    int sample_size = min((int)n_samples,isomerspace_forcefield::optimal_batch_size(N,0));

    std::vector<int> random_IDs(n_samples);
    std::iota(random_IDs.begin(), random_IDs.end(), 0);
    std::shuffle(random_IDs.begin(), random_IDs.end(), std::mt19937{42});
    auto id_range_end = min((int)n_samples, sample_size);
    std::vector<int> id_subset(random_IDs.begin(), random_IDs.begin()+n_samples);

    IsomerQueue input_queue(N, 0);
    IsomerQueue output_queue(N,0);
    IsomerQueue opt_queue(N, 0);
    input_queue.resize(sample_size*2);
    output_queue.resize(n_samples*2);
    opt_queue.resize(n_samples*2);
    


    IsomerBatch input_test(N, sample_size*2, DEVICE_BUFFER, 0);
    IsomerBatch h_input_test(N, sample_size*2, HOST_BUFFER, 0);
    IsomerBatch opt_test(N, sample_size, DEVICE_BUFFER, 0);
    IsomerBatch h_opt_test(N, sample_size, HOST_BUFFER, 0);

    IsomerBatch d_control(N,n_samples, DEVICE_BUFFER, 0);
    IsomerBatch h_control(N,n_samples, HOST_BUFFER, 0);
    
    LaunchCtx insert_ctx(0);

    for (int i = 0; i < n_samples; ++i){
        for (size_t j = 0; j < Nf; j++){
            G.neighbours[j].clear();
            for (size_t k = 0; k < 6; k++) {
                auto u = dual_neighbours[id_subset[i]*Nf*6 + j*6 +k];
                if(u != UINT16_MAX) G.neighbours[j].push_back(u);
            }
            for (size_t k = 0; k < G.neighbours[j].size(); ++k){
                h_control.dual_neighbours[i*Nf*6 + j*6 + k] = G.neighbours[j][k];
                h_control.face_degrees[i*Nf + j] = G.neighbours[j].size();
                h_control.statuses[i] = IsomerStatus::NOT_CONVERGED;
                h_control.IDs[i] = i;
            }
        }
    }
    //===================== SERIAL PIPELINE =====================
    cuda_io::copy(d_control,h_control);
    isomerspace_dual::dualise(d_control);
    isomerspace_tutte::tutte_layout(d_control, 1000000);
    isomerspace_X0::zero_order_geometry(d_control, 4.0);
    reset_convergence_statuses(d_control);

    isomerspace_forcefield::optimise<PEDERSEN>(d_control, N*50, N*50);
    //===================== END of SERIAL =====================

    int I_async = 0;
    auto generate_isomers = [&](int M){
        if(I_async == n_samples) return false;
        for (int i = 0; i < M; i++){
            if (I_async < n_samples){
                for (size_t j = 0; j < Nf; j++){
                    G.neighbours[j].clear();
                    for (size_t k = 0; k < 6; k++) {
                        auto u = dual_neighbours[id_subset[I_async]*Nf*6 + j*6 +k];
                        if(u != UINT16_MAX) G.neighbours[j].push_back(u);
                    }
                }
                input_queue.insert(G,I_async, insert_ctx, ASYNC);
                I_async++;
            } 
        }
        
        input_queue.refill_batch(input_test, insert_ctx, ASYNC);
        isomerspace_dual::dualise(input_test, insert_ctx, ASYNC);
        isomerspace_tutte::tutte_layout(input_test, 1000000, insert_ctx, ASYNC);
        isomerspace_X0::zero_order_geometry(input_test, 4.0, insert_ctx, ASYNC);
        reset_convergence_statuses(input_test, insert_ctx, ASYNC);
        insert_ctx.wait();
        
        return I_async < n_samples;
    };

    generate_isomers(sample_size*2);
    LaunchCtx device0(0);
    
    opt_queue.insert(input_test, device0, SYNC);
    opt_queue.refill_batch(opt_test, device0, SYNC);
    bool more_to_do = true;
    bool more_to_generate = false;
    auto step = max(1, (int)N/10);
    while (more_to_do){
        bool optimise_more = true;
        auto generate_handle = std::async(std::launch::async,generate_isomers, opt_test.isomer_capacity*2);
        while(optimise_more){
            isomerspace_forcefield::optimise<PEDERSEN>(opt_test,step, N*50, device0, ASYNC);
            output_queue.push_done(opt_test, device0, ASYNC);
            opt_queue.refill_batch(opt_test, device0, ASYNC);
            device0.wait();
            optimise_more = opt_queue.get_size() >= opt_test.isomer_capacity;
        }
        device0.wait();
        //cuda_io::copy(h_opt_test, opt_test);
        generate_handle.wait();
        more_to_generate = generate_handle.get();
        opt_queue.insert(input_test, device0, ASYNC);
        device0.wait();

        if(!more_to_generate){
            while(opt_queue.get_size() > 0){
                isomerspace_forcefield::optimise<PEDERSEN>(opt_test, step, N*50, device0, ASYNC);
                output_queue.push_done(opt_test, device0, ASYNC);
                device0.wait();
            
                opt_queue.refill_batch(opt_test, device0, ASYNC);
            }
            for(int i = 0;  i <  N*50; i += step){
                isomerspace_forcefield::optimise<PEDERSEN>(opt_test,step, N*50, device0, SYNC);
            }
            output_queue.push_done(opt_test, device0, SYNC);
            more_to_do = false;
        }

    }

    copy(output_queue.host_batch , output_queue.device_batch);
    IsomerBatch& output = output_queue.host_batch;



    CuArray<float> control_bonds(d_control.isomer_capacity*N*3);
    CuArray<float> control_angles(d_control.isomer_capacity*N*3);
    CuArray<float> control_dihedrals(d_control.isomer_capacity*N*3);
    CuArray<float> test_bonds(d_control.isomer_capacity*N*3);
    CuArray<float> test_angles(d_control.isomer_capacity*N*3);
    CuArray<float> test_dihedrals(d_control.isomer_capacity*N*3);

    output_queue.device_batch.shrink_to_fit();

    cuda_io::sort(d_control, IDS);
    cuda_io::sort(output_queue.device_batch, IDS);
    
    isomerspace_forcefield::get_bonds(d_control, control_bonds);
    isomerspace_forcefield::get_bonds(output_queue.device_batch, test_bonds);
    isomerspace_forcefield::get_angles(d_control, control_angles);
    isomerspace_forcefield::get_angles(output_queue.device_batch, test_angles);
    isomerspace_forcefield::get_dihedrals(d_control, control_dihedrals);
    isomerspace_forcefield::get_dihedrals(output_queue.device_batch, test_dihedrals);

    std::vector<float> rdiffs_bond(max((int)d_control.isomer_capacity,10000));
    std::vector<float> rdiffs_angle(max((int)d_control.isomer_capacity,10000));
    std::vector<float> rdiffs_dihedral(max((int)d_control.isomer_capacity,10000));
    std::vector<float> adiffs_bond(max((int)d_control.isomer_capacity,10000));
    std::vector<float> adiffs_angle(max((int)d_control.isomer_capacity,10000));
    std::vector<float> adiffs_dihedral(max((int)d_control.isomer_capacity,10000));

    cuda_io::copy(h_control, d_control);
    output_queue.host_batch.shrink_to_fit();
    cuda_io::copy(output_queue.host_batch, output_queue.device_batch);
    std::vector<int> iterations(max( (int)d_control.isomer_capacity, 10000));
    {
    for (size_t i = 0; i < output_queue.device_batch.isomer_capacity; i++)
        iterations[i] = output_queue.host_batch.iterations[i];
    }
    

    for (size_t i = 0; i < d_control.isomer_capacity; i++)
    {
        auto [n_correct_b,adiff_b, rdiff_b ] = cuda_io::compare_isomer_arrays(test_bonds.data + i*N*3, control_bonds.data +i*N*3,1, N*3, 1e-4);
        auto [n_correct_a,adiff_a, rdiff_a ] = cuda_io::compare_isomer_arrays(test_angles.data + i*N*3, control_angles.data +i*N*3,1, N*3, 1e-4);
        auto [n_correct_d,adiff_d, rdiff_d ] = cuda_io::compare_isomer_arrays(test_dihedrals.data + i*N*3, control_dihedrals.data +i*N*3,1, N*3, 1e-4);
        adiffs_bond[i] = adiff_b; adiffs_angle[i] = adiff_a; adiffs_dihedral[i] = adiff_d;
        rdiffs_bond[i] = rdiff_b; rdiffs_angle[i] = rdiff_a; rdiffs_dihedral[i] = rdiff_d;
    }
    file_bond_adiff << N << "," << count_batch_status(d_control, IsomerStatus::CONVERGED) + count_batch_status(d_control,IsomerStatus::FAILED) << ",";
    file_angle_adiff <<  N << "," << count_batch_status(d_control, IsomerStatus::CONVERGED) + count_batch_status(d_control,IsomerStatus::FAILED)  << ",";
    file_dihedral_adiff <<  N << "," << count_batch_status(d_control, IsomerStatus::CONVERGED) + count_batch_status(d_control,IsomerStatus::FAILED)  << ",";
    file_bond_rdiff <<  N << "," << count_batch_status(d_control, IsomerStatus::CONVERGED) + count_batch_status(d_control,IsomerStatus::FAILED)  << ",";
    file_angle_rdiff <<  N << "," << count_batch_status(d_control, IsomerStatus::CONVERGED) + count_batch_status(d_control,IsomerStatus::FAILED)  << ",";
    file_dihedral_rdiff <<  N << "," << count_batch_status(d_control, IsomerStatus::CONVERGED) + count_batch_status(d_control,IsomerStatus::FAILED)  << ",";
    file_iterations <<  N << "," << count_batch_status(d_control, IsomerStatus::CONVERGED) + count_batch_status(d_control,IsomerStatus::FAILED)  << ",";
    
    
    for (int i = 0; i < h_control.isomer_capacity; i++){
        if(h_control.IDs[i] != output.IDs[i]){
            output.print(IDS, {i-50, i+50});
            std::cout << i << "," << output.IDs[i] << endl;
            assert(h_control.IDs[i] == output.IDs[i]);
        }
    }


        for (int i = 0; i < adiffs_bond.size(); i++){
            if (i != adiffs_bond.size()-1) {
                file_bond_adiff << adiffs_bond[i] << ",";
                file_angle_adiff << adiffs_angle[i] << ",";
                file_dihedral_adiff << adiffs_dihedral[i] << ",";
                file_bond_rdiff << rdiffs_bond[i] << ",";
                file_angle_rdiff << rdiffs_angle[i] << ",";
                file_dihedral_rdiff << rdiffs_dihedral[i] << ",";
                file_iterations << iterations[i] <<  ",";
            }
            else{
                file_bond_adiff << adiffs_bond[i] << "\n";
                file_angle_adiff << adiffs_angle[i] << "\n";
                file_dihedral_adiff << adiffs_dihedral[i] << "\n";
                file_bond_rdiff << rdiffs_bond[i] << "\n";
                file_angle_rdiff << rdiffs_angle[i] << "\n";
                file_dihedral_rdiff << rdiffs_dihedral[i] << "\n";
                file_iterations << iterations[i] <<  "\n";
            } 
        }
        
    }
    
}
