#include "fullerenes/polyhedron.hh"
#include "fullerenes/buckygen-wrapper.hh"

#include "fullerenes/gpu/kernels.hh"
#include "fullerenes/device_io.hh"
#include "fullerenes/isomer_queue.hh"
#include "fullerenes/isomerdb.hh"

int main(int argc, char** argv){
    const size_t N_start                = argc>1 ? strtol(argv[1],0,0) : 20;     // Argument 1: Number of vertices 
    const size_t N_end                = argc>2 ? strtol(argv[2],0,0) : N_start;     // Argument 1: Number of vertices 
    const size_t N_samples                = argc>3 ? strtol(argv[3],0,0) : 10000;     // Argument 1: Number of vertices 
    
    for (size_t N = N_start; N < N_end+1; N+=2)
    {    
        bool more_to_do = true;
        auto n_isomers = IsomerDB::number_isomers(N);
        FullereneDual G;
        auto Nf = N/2 +2;

        std::string path = "isomerspace_samples/dual_layout_" + to_string(N) + "_seed_42";
        auto fsize = file_size(path);
        auto n_samples = fsize / (Nf * 6 * sizeof(device_node_t));
        ifstream in_file(path,std::ios::binary);
        std::vector<device_node_t> dual_neighbours(n_samples * Nf * 6);
        in_file.read((char*)dual_neighbours.data(),n_samples * Nf * 6 * sizeof(device_node_t));

        G.neighbours = neighbours_t(Nf, std::vector<node_t>(6));
        G.neighbours.resize(Nf);
        G.N = Nf;

        auto batch_size = min((int)n_samples, (int)n_isomers);
        device_io::IsomerQueue test_queue(N);
        device_io::IsomerQueue validation_queue(N);
        IsomerBatch<CPU> h_test(N,batch_size);
        
        IsomerBatch<CPU> h_validation(N,batch_size);
        IsomerBatch<GPU> d_test(N,batch_size);
        IsomerBatch<GPU> d_validation(N,batch_size);
        

        for (size_t I = 0; I < n_samples; I++)
        {
            for (size_t j = 0; j < Nf; j++)
            {
                G.neighbours[j].clear();
                for (size_t k = 0; k < 6; k++) {
                    auto u = dual_neighbours[I*Nf*6 + j*6 +k];
                    if(u != UINT16_MAX) G.neighbours[j].push_back(u);
                }
            }
            test_queue.insert(Graph(G),I);
            G.update();
            auto FD = G.dual_graph();
            validation_queue.insert(FD,I, LaunchCtx(), LaunchPolicy::SYNC, false);
        }
        //std::cout << F.neighbours << "\n";
        test_queue.refill_batch(d_test);                        validation_queue.refill_batch(d_validation);
        device_io::copy(h_test, d_test);
        gpu_kernels::isomerspace_dual::dualise_4(h_test);    

        //device_io::copy(h_test, d_test); 
        device_io::copy(h_validation, d_validation);
        //h_test.set_print_verbose(); h_validation.set_print_verbose();

        if(h_test == h_validation) std::cout << "Test passed!\n";
        else std::cout << "Test Failed" << std::endl;

    
    }

}
