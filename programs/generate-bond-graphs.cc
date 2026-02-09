#include <fullerenes/argparser.hh>
#include <fullerenes/sycl-headers/sycl-fullerene-structs.hh>
#include <fullerenes/sycl-headers/sycl-fullerene-structs.hh>
#include <fullerenes/sycl-headers/sycl-device-queue.hh>
#include <fullerenes/sycl-headers/sycl-vector.hh>
#include <fullerenes/kernel-headers/all-kernels.hh>
#include <fullerenes/buckygen-wrapper.hh>

int main(int argc, char *argv[]) {
    // Parse command line arguments
    CmdArgs input_args;
    parseArguments(argc, argv, input_args);

    int N = input_args.natoms; // Number of atoms in the fullerene
    auto capacity = input_args.nisomers; // Capacity of the queue
    std::string output_file = input_args.output_file;

    // Create a FullereneQueue with the specified capacity
    FullereneQueue<double, uint16_t> queue(N, capacity);
    auto buckyQ = BuckyGen::start(N, false, false);
    // Generate bond graphs and push them into the queue
    bool more = true;
    Graph g(N/2 + 2);
    std::cout << "Generating.." << std::endl;
    while (more && queue.size() < queue.capacity()) {
        more = BuckyGen::next_fullerene(buckyQ, g);
        queue.push_back(g);
    }
    FullereneBatch<double, uint16_t> batch(N, capacity);
    std::cout << "Dualizing.." << std::endl;
    
    DualizeFunctor<double, uint16_t> dualizer;
    std::cout << "Computing Tutte Embeddings.." << std::endl;
    TutteFunctor<double, uint16_t> tutter;
    std::cout << "Computing Projections.." << std::endl;
    SphericalProjectionFunctor<double, uint16_t> projector;
    std::cout << "Optimizing Forcefield.." << std::endl;
    ForcefieldOptimizeFunctor<ForcefieldType::PEDERSEN, double, uint16_t> forcefield_optimizer;

    SyclQueue ctx(Device::default_device(), true);
    QueueUtil::push(ctx, batch, queue);
    dualizer(ctx, batch, LaunchPolicy::SYNC);
    tutter(ctx, batch, LaunchPolicy::SYNC);
    projector(ctx, batch, LaunchPolicy::SYNC);
    forcefield_optimizer(ctx, batch, LaunchPolicy::SYNC, 5*N, 5*N);
    
    ctx.wait();
    std::cout << count_if(batch.d_.X_cubic_.begin(), batch.d_.X_cubic_.end(), [](const auto number) {
        return std::isnan(number[0]) || std::isnan(number[1]) || std::isnan(number[2]);
    }) << " NaNs found" << std::endl;

    std::cout << "Writing coordinates to binary file: " << output_file << std::endl;
    std::ofstream out(output_file, std::ios::binary);
    if (!out) {
        std::cerr << "Error opening output file: " << output_file << std::endl;
        return 1;
    }

    // Get the number of structures in the queue
    size_t num_structures = queue.size();
    out.write(reinterpret_cast<const char*>(batch.d_.X_cubic_.data()), batch.d_.X_cubic_.to_span().size_bytes());
    out.close();


    return 0;
}
