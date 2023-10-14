#include "numeric"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/buckygen-wrapper.hh"

#include "fullerenes/gpu/kernels.hh"
#include "fullerenes/device_io.hh"
#include "fullerenes/isomer_queue.hh"
#include "fullerenes/isomerdb.hh"

int main(int argc, char** argv){
    const size_t N_start                =   argc > 1 ? strtol(argv[1],0,0) : 20;     // Argument 1: Number of vertices 
    const size_t N_end                  =   argc > 2 ? strtol(argv[2],0,0) : N_start;     // Argument 1: Number of vertices 

    ofstream rel_file("tutte_validation_rel_" + to_string(N_start) + "_" + to_string(N_end) + ".txt");
    ofstream abs_file("tutte_validation_abs_" + to_string(N_start) + "_" + to_string(N_end) + ".txt");
    ofstream mag_file("tutte_validation_mag_" + to_string(N_start) + "_" + to_string(N_end) + ".txt");

    for (size_t n = N_start; n < N_end + 1; n+=2)
    {
        bool more_to_do = true;
        auto n_isomers = IsomerDB::number_isomers(n);
        auto bucky_queue  = BuckyGen::start(n, false, false);
        auto max_samples = 10000;
        auto batch_size = min(max_samples, (int)n_isomers);
        device_io::IsomerQueue test_queue(n);
        device_io::IsomerQueue validation_queue(n);
        IsomerBatch<CPU> h_test(n,batch_size);
        IsomerBatch<CPU> h_validation(n,batch_size);
        IsomerBatch<GPU> d_test(n,batch_size);
        IsomerBatch<GPU> d_validation(n,batch_size);
        FullereneDual F;
        auto I = 0;
        while(more_to_do && I < batch_size){

            more_to_do &= BuckyGen::next_fullerene(bucky_queue,F);
            if(!more_to_do)break;
            F.update();
            auto FD = F.dual_graph();
            test_queue.insert(FD,I, LaunchCtx(), LaunchPolicy::SYNC, false);
            FD.layout2d =  FD.tutte_layout();
            validation_queue.insert(FD,I, LaunchCtx(), LaunchPolicy::SYNC, true);
            ++I;
        }
        test_queue.refill_batch(d_test);                        validation_queue.refill_batch(d_validation);
 
        gpu_kernels::isomerspace_tutte::tutte_layout(d_test);

        device_io::copy(h_test, d_test); device_io::copy(h_validation, d_validation);
        device_io::sort(h_test); device_io::sort(h_validation);
        //std::cout << h_test;
        //std::cout << h_validation;
        std::vector<float> rdiffs(max_samples);
        std::vector<float> adiffs(max_samples);
        for (size_t j = 0; j < batch_size; j++)
        {
            auto [a, b, rdiff] = device_io::compare_isomer_arrays(h_test.xys + 2*n*j, h_validation.xys + 2*n*j, 1, n);
            rdiffs[j] = rdiff;
            adiffs[j] = b;
        }

        float magnitude = 0.0;
        for (size_t j = 0; j < batch_size*n*2; j++)
        {
            magnitude += std::abs(h_validation.xys[j]);
        }
        magnitude /= batch_size*n*2;
        
        std::cout << "Isomerspace: " << n << " Magnitutde:" << magnitude << "\n";

        std::cout << "Isomerspace: " << n << " RErr:" << sum(rdiffs)/(float)batch_size << "\n";
        std::cout << "Isomerspace: " << n << " AErr:" << sum(adiffs)/(float)batch_size << "\n";

        rel_file << batch_size << ", ";
        abs_file << batch_size << ", ";

        for (size_t j = 0; j < max_samples; j++){
            if (j != max_samples-1) {
                rel_file << rdiffs[j] << ", ";
                abs_file << adiffs[j] << ", ";
            }
            else{
                rel_file << rdiffs[j] << "\n";
                abs_file << adiffs[j] << "\n";
            } 
        }
        
    }
    
    
}
