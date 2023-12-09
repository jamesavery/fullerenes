#include "string"
#include "algorithm"
#include "iostream"
#include "numeric"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/isomerdb.hh"
#include "fullerenes/gpu/benchmark_functions.hh"
#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/gpu/kernels.hh"
#include "fullerenes/device_io.hh"
#include "fullerenes/isomer_queue.hh"

using namespace gpu_kernels;
using namespace isomerspace_dual;
using namespace isomerspace_eigen;
using namespace isomerspace_forcefield;
using namespace isomerspace_hessian;
using namespace isomerspace_X0;
using namespace isomerspace_tutte;

int main(int argc, char** argv){
    size_t N = argc > 1 ? atoi(argv[1]) : 60;
    auto n_isomers = IsomerDB::number_isomers(N);
    Graph G;
    auto Nf = N/2 +2;
    G.neighbours = neighbours_t(Nf, std::vector<node_t>(6));
    G.neighbours.resize(Nf);
    G.N = Nf;
    PlanarGraph Pg;
    auto batch_size = min(isomerspace_forcefield::optimal_batch_size(N)*8, (int)n_isomers);
    //BuckyGen::buckygen_queue BuckyQ = BuckyGen::start(N,0,0);  
    IsomerBatch<CPU> Bhost(N,batch_size);
    IsomerBatch<CPU> Bstart(N,batch_size);
    IsomerBatch<GPU> Bdev(N,batch_size);
    for (size_t i = 0; i < batch_size; i++)
    {
        auto ID = cuda_benchmark::random_isomer("isomerspace_samples/dual_layout_"+to_string(N)+"_seed_42", G);
        Bhost.append(G,ID);
    }
    device_io::copy(Bdev, Bhost);
    dualize(Bdev);
    tutte_layout(Bdev, N*10);
    zero_order_geometry(Bdev, 4.0);
    device_io::copy(Bstart, Bdev);
    optimize<PEDERSEN>(Bdev, N*5, N*5);
    std::string Before_ = "Before_"+to_string(N);
    std::string After_ = "After_"+to_string(N);
    std::string Most_ = "Most_Eccentric_"+to_string(N);
    //FILE* f = fopen((Before_ + "_Geometry.mol2").c_str(), "w");
  //  FILE* g = fopen((After_ + "_Geometry.mol2").c_str(), "w");
    device_io::copy(Bhost, Bdev);
    CuArray<device_real_t> orthogonalities(batch_size);
    CuArray<device_real_t> eigs(batch_size*3);
    CuArray<device_real_t> eigvecs(batch_size*9);
    CuArray<device_real_t> inertia(batch_size*9);
    isomerspace_properties::debug_function(Bdev, eigs, eigvecs, inertia, orthogonalities);

    //isomerspace_properties::eccentricities(Bdev, ecce);
    //isomerspace_properties::transform_coordinates(Bdev);
    device_io::copy(Bhost, Bdev);
    std::vector<int> sorted_indices(batch_size);
    std::iota(sorted_indices.begin(), sorted_indices.end(), 0);
    std::sort(sorted_indices.begin(), sorted_indices.end(), [&orthogonalities](int i1, int i2) {return orthogonalities[i1] < orthogonalities[i2];});
    ofstream badeigs("Eigenvals_" + to_string(N) + ".float32", std::ios::binary);
    ofstream badeigvecs("Eigenvecs_" + to_string(N) + ".float32", std::ios::binary);
    ofstream badinertia("Inertia_" + to_string(N) + ".float32", std::ios::binary);
    ofstream optgeom("OptGeom_" + to_string(N) + ".float32", std::ios::binary);
    ofstream graphs("Graphs_" + to_string(N) + ".uint16", std::ios::binary);
    ofstream startgeom("StartGeom_" + to_string(N) + ".float32", std::ios::binary);
    badeigs.write(reinterpret_cast<char*>(eigs.data), batch_size*3*sizeof(device_real_t));
    badeigvecs.write(reinterpret_cast<char*>(eigvecs.data), batch_size*9*sizeof(device_real_t));
    badinertia.write(reinterpret_cast<char*>(inertia.data), batch_size*9*sizeof(device_real_t));
    optgeom.write(reinterpret_cast<char*>(Bhost.X), batch_size*N*3*sizeof(device_real_t));
    startgeom.write(reinterpret_cast<char*>(Bstart.X), batch_size*N*3*sizeof(device_real_t));
    graphs.write(reinterpret_cast<char*>(Bhost.cubic_neighbours), batch_size*N*3*sizeof(device_node_t));

    auto P = Bhost.get_isomer(0).value();
    P.volume_divergence();
    CuArray<device_real_t> VD(batch_size);
    isomerspace_properties::volume_divergences(Bdev, VD);
    std::cout << "Volume Divergence: " << VD << std::endl;
    std::cout << "Volume Divergence: " << P.volume_divergence() << std::endl;

    CuArray<device_real_t> SA(batch_size);
    isomerspace_properties::surface_areas(Bdev, SA);
    std::cout << "Surface Area: " << SA << std::endl;
    std::cout << "Surface Area: " << P.surface_area() << std::endl;
    
    
    //std::cout << "Worst ID: " << sorted_indices[batch_size-1] << std::endl;
    //std::cout << "Wors Inertia" << std::vector<double>(inertia.data + sorted_indices[batch_size-1]*9, inertia.data + sorted_indices[batch_size-1]*9 + 9) << std::endl;
    //std::cout << "Worst Orthogonality: " << orthogonalities[sorted_indices[batch_size-1]] << std::endl;
    //std::sort(orthogonalities.data, orthogonalities.data+batch_size);
    //std::cout << "Orthogonality:\n " << orthogonalities << std::endl;


    //std::cout << eigs << std::endl;
    //std::cout << "Indices:\n " << sorted_indices << std::endl;
}
