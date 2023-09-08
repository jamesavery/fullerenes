#include "numeric"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/gpu/benchmark_functions.hh"
#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/gpu/kernels.hh"
#include "fullerenes/device_io.hh"
#include "fullerenes/isomer_queue.hh"
#include "fullerenes/progress_bar.hh"
#include "fullerenes/isomerdb.hh"

using namespace gpu_kernels;
using namespace isomerspace_dual;
using namespace isomerspace_eigen;
using namespace isomerspace_forcefield;
using namespace isomerspace_hessian;
using namespace isomerspace_X0;
using namespace isomerspace_tutte;

constexpr device_real_t carbon_mass = 1.9944733e-26/*kg*/, aangstrom_length = 1e-10/*m*/;
int main(int argc, char** argv){
    const size_t N                = argc>1 ? strtol(argv[1],0,0) : 60;     // Argument 1: Number of vertices 
    const size_t Mlanczos_steps   = argc>2 ? strtol(argv[2],0,0) : 40;     // Argument 2: Number of Lanczos steps
    float reldelta                = argc>3 ? strtof(argv[3],0) : 1e-5;    // Argument 3: Relative delta
    int isomer_num                = argc>4 ? strtol(argv[4],0,0) : 0;     // Argument 4: Isomer number
    std::string spiral_           = argc>5 ? argv[5] : "C60-[1,7,9,11,13,15,18,20,22,24,26,32]-fullerene";          // Argument 2: Spiral
    std::string name_             = argc>6 ? argv[6] : "C60ih";          // Argument 2: Spiral
    std::string filename          = argc>7 ? argv[7] : "hessian_validation";        // Argument 2: Filename
    
    std::string type = "float32";
    std::ofstream hess_analytical(filename + "_analytical." + type);
    std::ofstream hess_numerical(filename + "_numerical." + type);
    std::ofstream hess_cols(filename + "_cols.uint16");
    std::ofstream cubic_graph;//("C" + to_string(N) + "_CubicGraph_" + to_string(isomer_num) + ".bin");
    std::ofstream geometry;//("C" + to_string(N) + "_Geometry_" + to_string(isomer_num) + ".bin");
    

    bool more_to_do = true;
    auto n_isomers = IsomerDB::number_isomers(N);
    Graph G;
    auto Nf = N/2 +2;
    G.neighbours = neighbours_t(Nf, std::vector<node_t>(6));
    G.neighbours.resize(Nf);
    G.N = Nf;
    auto batch_size = min(2000, (int)n_isomers);
    if (isomer_num == -1) batch_size = 1;

    CuArray<device_real_t> hessians(N*3*3*10 * batch_size);
    CuArray<device_real_t> hessians_fd(N*3*3*10 * batch_size);
    CuArray<device_node_t> cols(N*3*3*10 * batch_size);
    CuArray<device_real_t> lambda_mins(batch_size);
    CuArray<device_real_t> lambda_maxs(batch_size);
    CuArray<device_real_t> Q(N*3*N*3*batch_size);
    CuArray<device_real_t> eigs(N*3*batch_size);
    CuArray<device_real_t> min_eigvects(N*3*batch_size);
    CuArray<device_real_t> max_eigvects(N*3*batch_size); 

    PlanarGraph Pg;
    //BuckyGen::buckygen_queue BuckyQ = BuckyGen::start(N,0,0);  
    //IsomerQueue Q0(N);
    IsomerBatch<CPU> Bhost(N,batch_size);

    int Nd = LaunchCtx::get_device_count();
    IsomerBatch<GPU> Bdev(N,batch_size,0);
    BuckyGen::buckygen_queue BuckyQ = BuckyGen::start(N, false, false);  
    if (isomer_num == -1) {
        spiral_nomenclature C60name(spiral_);    
        Triangulation C60dual(C60name);
        PlanarGraph C60cubic = C60dual.dual_graph();
        for(int i = 0;  i < batch_size; i++) Bhost.append(C60cubic,0, false);
    } else {
        
        for (size_t i = 0; i < batch_size; i++)
        {   
            bool success = BuckyGen::next_fullerene(BuckyQ, G);
            if(success) Bhost.append(G,i);
            else break;
        }
    }

    
    LaunchCtx ctx(0);
    LaunchPolicy policy = LaunchPolicy::SYNC;
    ProgressBar progbar('=', 50);
    //std::cout << vector<device_real_t>(Bhost.cubic_neighbours, Bhost.cubic_neighbours + N*3) << std::endl;
    device_io::copy(Bdev, Bhost);
    //device_io::copy(Bhost, Bdev);
    //ifstream geometry_in("X.float64", std::ios::binary); geometry_in.read((char*)Bhost.X, N*3*batch_size*sizeof(float));
/*     for (size_t i = 0; i < batch_size; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            for (size_t k = 0; k < 3; k++)
            {
                Bhost.X[i*N*3 + j*3 + k] = device_real_t(X0[i*N*3 + j*3 + k]);
            }
        }
    } */
    ofstream fullspectrum_file("spectrum." + type, std::ios::binary);
    ofstream lambda_mins_file("spectrum_mins." + type, std::ios::binary);
    ofstream lambda_maxs_file("spectrum_maxs." + type, std::ios::binary);
    if(isomer_num != -1) dualise(Bdev);
    tutte_layout(Bdev, (int)20*N);
    zero_order_geometry(Bdev, 4.0);
    optimise<PEDERSEN>(Bdev, 5*N, 6*N);
    optimise<PEDERSEN>(Bdev, 1*N, 6*N);
    isomerspace_properties::transform_coordinates(Bdev);
    compute_hessians<PEDERSEN>(Bdev, hessians, cols);
    CuArray<device_real_t> vols(batch_size);
    isomerspace_properties::eccentricities(Bdev, vols);
    std::vector<device_real_t> vols_host(vols.data, vols.data + vols.size());
    std::vector<device_node_t> indices(vols.size());
    device_io::copy(Bhost, Bdev); 
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&vols_host, &Bhost](const auto& a, const auto& b) {
        return std::abs(vols_host[a]) < std::abs(vols_host[b]) && Bhost.statuses[a] == IsomerStatus::CONVERGED;
    });
    auto min_index = indices[0];
    
    std::cout << "Min volume: " << vols[min_index] << std::endl;
    std::cout << "Min volume index: " << min_index << std::endl;

    Polyhedron P = Bhost.get_isomer(min_index).value();
    Polyhedron::to_file(P, "MinVol.mol2");





    hess_analytical.write((char*)hessians.data, hessians.size()*sizeof(device_real_t));
    hess_cols.write((char*)cols.data, cols.size()*sizeof(device_node_t));
    int lanczos_max = std::min(300, (int)N*3 - 6);
    //spectrum_ends(Bdev, hessians, cols, lambda_mins, lambda_maxs, min_eigvects, max_eigvects, 40);
    for (int i = 10; i < lanczos_max; i++){
        
        spectrum_ends(Bdev, hessians, cols, lambda_mins, lambda_maxs, min_eigvects, max_eigvects, i);
        lambda_mins_file.write((char*)lambda_mins.data, lambda_mins.size()*sizeof(device_real_t));
        lambda_maxs_file.write((char*)lambda_maxs.data, lambda_maxs.size()*sizeof(device_real_t));
        progbar.update_progress((i-9)/float(lanczos_max-10));
    }

    eigensolve(Bdev, Q, hessians, cols, eigs);
    //auto eigs_min_ecc = vector<device_real_t>(eigs.data + N*3 * min_index, eigs.data + N*3 * min_index + N*3); 
    /* std::sort(eigs_min_ecc.begin(), eigs_min_ecc.end(), std::less<device_real_t>());
    std::transform(eigs_min_ecc.begin(), eigs_min_ecc.end(), eigs_min_ecc.begin(), [](auto& x){return  33.356 * sqrt(x/carbon_mass)/(2*M_PI)*1e-12;   });
    std::cout << eigs_min_ecc << std::endl; */
    //std::cout << "From spectrum ends: " << 33.356 * sqrt(lambda_mins[min_index]/carbon_mass)/(2*M_PI)*1e-12 << std::endl;
    fullspectrum_file.write((char*)eigs.data, eigs.size()*sizeof(device_real_t));
    //isomerspace_properties::transform_coordinates(Bdev);
    
    if (isomer_num == -1){
        std::cout << eigs << endl;
        ofstream file(spiral_ + "_hessian." + type); file.write((char*)hessians.data, 30*N*3*sizeof(device_real_t));
        ofstream file2(spiral_ + "_cols.uint16" ); file2.write((char*)cols.data, 30*N*3*sizeof(device_node_t));
        ofstream file3(spiral_ + "_eigs." + type); file3.write((char*)eigs.data, N*3*sizeof(device_real_t));
        ofstream file4(spiral_ + "_eigvects." + type); file4.write((char*)Q.data, N*3*N*3*sizeof(device_real_t));
        ofstream file5(spiral_ + "_X." + type); file5.write((char*)Bhost.X, N*3*sizeof(device_real_t));
        ofstream file6(spiral_ + "_adjacency." + type); file6.write((char*)Bhost.cubic_neighbours, N*3*sizeof(device_node_t));
    }

    //compute_hessians_fd<PEDERSEN>(Bdev, hessians_fd, cols, reldelta);

    
    hess_analytical.write((char*)hessians.data, hessians.size()*sizeof(device_real_t));
    hess_cols.write((char*)cols.data, cols.size()*sizeof(device_node_t));
    hess_numerical.write((char*)hessians_fd.data, hessians_fd.size()*sizeof(device_real_t));


    
    //std::cout << "Eigenvalues: " << eigenvalues << std::endl;
    //device_io::copy(Bhost, Bdev); // Copy back to host
    //auto FG = Bhost.get_isomer(0).value();
    //cubic_graph.write((char*)Bhost.cubic_neighbours, N*3*sizeof(device_node_t));
    //geometry.write((char*)Bhost.X, N*3*sizeof(device_real_t));
    
    //eigensolve_cusolver(Bdev, hessians, cols, eigenvalues);
    //eigensolve_cusolver(Bdev, hessians, cols, eigenvalues);
    //std::cout << "cuSOLVE Eigensolver took: " << (time/1us)/(float)batch_size << " us / graph" << std::endl;
       
    //lambda_max(Bdev, hessians, cols, lambda_maxs);
    auto start = std::chrono::steady_clock::now();
    
    auto time = std::chrono::steady_clock::now() - start;
    std::cout << "Spectrum Ends took: " << (time/1us)/(float)batch_size << " us / graph" << std::endl;

    start = std::chrono::steady_clock::now();
    //eigensolve(Bdev, Q, hessians, cols, eigs);
    //ofstream eigs_out("eigs.float32", ios::binary); eigs_out.write((char*)eigs.data, N*3*batch_size*sizeof(device_real_t)); eigs_out.close();
    time = std::chrono::steady_clock::now() - start;
    std::cout << "Full eigensolve took: " << (time/1us)/(float)batch_size << " us / graph" << std::endl;
    device_io::copy(Bhost, Bdev); // Copy back to host

    std::vector<device_real_t> min_eigs_ref(batch_size), max_eigs_ref(batch_size);
    
    for(int i = 0; i < batch_size; i++) {
        std::sort(eigs.data + i*N*3, eigs.data + (i+1)*N*3);
        min_eigs_ref[i] = eigs[i*N*3 + 6];
        max_eigs_ref[i] = eigs[(i+1)*N*3 - 1];
    }
    std::vector<device_real_t> rel_err_min(batch_size), abs_err_min(batch_size), rel_err_max(batch_size), abs_err_max(batch_size), rel_err_width(batch_size), abs_err_width(batch_size);
    std::vector<int> nan_or_inf;
    for(int i = 0; i < batch_size; i++) {
        if (std::isnan(lambda_mins[i]) || std::isinf(lambda_mins[i]) || std::isnan(lambda_maxs[i]) || std::isinf(lambda_maxs[i]) || std::isnan(min_eigs_ref[i]) || std::isinf(min_eigs_ref[i]) || std::isnan(max_eigs_ref[i]) || std::isinf(max_eigs_ref[i]) || Bhost.statuses[i] == IsomerStatus::FAILED) {
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
    std::cout << "MinLambdaError_Relative: " << rel_err_min << std::endl;
    std::cout << "MaxLambdaError_Relative: " << rel_err_max << std::endl;
    //std::cout << "MinLambdaError_Absolute: " << abs_err_min << std::endl;

    std::cout << "Failed: " << nan_or_inf << std::endl;
    if (isomer_num == -1){
        isomerspace_properties::transform_coordinates(Bdev);
        device_io::copy(Bhost, Bdev); // Copy back to host
        Polyhedron::to_file(Bhost.get_isomer(0).value(), spiral_ + ".mol2");
    }

    

    /* ofstream fileX("X."+ type, ios::out | ios::binary);
    fileX.write((char*)Bhost.X, batch_size*3*N*sizeof(float));

    ofstream fileA("A.uint16", ios::out | ios::binary);
    fileA.write((char*)Bhost.cubic_neighbours, batch_size*3*N*sizeof(device_node_t));

    ofstream file("MinLambdaError_Relative."+ type, ios::out | ios::binary);
    file.write((char*)rel_err_min.data(), batch_size*sizeof(device_real_t));

    ofstream file2("MinLambdaError_Absolute."+ type, ios::out | ios::binary);
    file2.write((char*)abs_err_min.data(), batch_size*sizeof(device_real_t));

    ofstream file3("MaxLambdaError_Relative."+ type, ios::out | ios::binary);
    file3.write((char*)rel_err_max.data(), batch_size*sizeof(device_real_t));

    ofstream file4("MaxLambdaError_Absolute."+ type, ios::out | ios::binary);
    file4.write((char*)abs_err_max.data(), batch_size*sizeof(device_real_t));

    ofstream file5("WidthError_Relative."+ type, ios::out | ios::binary);
    file5.write((char*)rel_err_width.data(), batch_size*sizeof(device_real_t));

    ofstream file6("WidthError_Absolute."+ type, ios::out | ios::binary);
    file6.write((char*)abs_err_width.data(), batch_size*sizeof(device_real_t));

    ofstream file7("Q."+ type, ios::out | ios::binary);
    file7.write((char*)Q.data, Q.size()*sizeof(device_real_t));

    ofstream file8("MinEigVects."+ type, ios::out | ios::binary);
    file8.write((char*)min_eigvects.data, min_eigvects.size()*sizeof(device_real_t));

    ofstream file9("MaxEigVects."+ type, ios::out | ios::binary);
    file9.write((char*)max_eigvects.data, max_eigvects.size()*sizeof(device_real_t));

    ofstream file10("MinEigVals."+ type, ios::out | ios::binary);
    file10.write((char*)lambda_mins.data, lambda_mins.size()*sizeof(device_real_t));

    ofstream file11("MaxEigVals."+ type, ios::out | ios::binary);
    file11.write((char*)lambda_maxs.data, lambda_maxs.size()*sizeof(device_real_t)); */

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
        //LaunchCtx::clear_allocations();
    }