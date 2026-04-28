#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <fullerenes/argparser.hh>
#include <fullerenes/buckygen-wrapper.hh>
#include <fullerenes/graph.hh>

#include <fullerenes/batch/batch.hh>
#include <fullerenes/batch/batch_state.hh>
#include <fullerenes/batch/utilities.hh>
#include <fullerenes/kernel-headers/all-kernels.hh>
#include <fullerenes/sycl-headers/sycl-device-queue.hh>
#include <fullerenes/sycl-headers/sycl-vector.hh>

int main(int argc, char *argv[]) {
    CmdArgs input_args;
    parseArguments(argc, argv, input_args);

    using real_t  = double;
    using index_t = std::uint16_t;
    using RSR     = Spanify::RSRAdjacencyView<index_t>;

    const int N        = input_args.natoms;
    const int Nf       = N/2 + 2;
    const int capacity = input_args.nisomers;
    const std::string output_file = input_args.output_file;

    // ---- Allocate batch buffers ----------------------------------------
    batch::Batch<RSR>   src_dual (Nf, capacity, /*dmax*/6);
    batch::Batch<RSR>   dst_cubic(N,  capacity, /*dmax*/3);
    batch::BatchState   st(capacity);

    SyclVector<std::array<index_t,6>> faces_cubic_buf(capacity * Nf);
    SyclVector<std::array<index_t,3>> faces_dual_buf (capacity * N);
    SyclVector<std::array<real_t,2>>  layout2d       (capacity * N);
    SyclVector<std::array<real_t,3>>  xyz            (capacity * N);

    // ---- Fill src_dual from BuckyGen ------------------------------------
    std::cout << "Generating.." << std::endl;

    auto BQ = BuckyGen::start(N, false, false);
    Graph G(Nf, std::vector<node_t>(6, -1));
    int pushed = 0;
    for (int i = 0; i < capacity; ++i) {
        if (!BuckyGen::next_fullerene(BQ, G)) break;
        batch::push_back(src_dual, G);
        st.push_back(uint64_t(i),
                     StatusFlag(int(StatusEnum::DUAL_INITIALIZED)));
        ++pushed;
    }
    BuckyGen::stop(BQ);
    std::cout << "Generated " << pushed << " dual graphs" << std::endl;

    // ---- Pipeline --------------------------------------------------------
    SyclQueue ctx(Device::default_device(), true);
    auto src_view = src_dual .view_capacity();
    auto dst_view = dst_cubic.view_capacity();

    std::cout << "Dualizing.." << std::endl;
    dualize<real_t, index_t>(ctx, src_view, dst_view, st.view(),
        std::span<std::array<index_t,6>>(faces_cubic_buf.data(), faces_cubic_buf.size()),
        std::span<std::array<index_t,3>>(faces_dual_buf .data(), faces_dual_buf .size())
    ).wait();
    std::cout << "Computing Tutte Embeddings.." << std::endl;
    tutte_layout<real_t, index_t>(ctx, dst_view,
        std::span<std::array<real_t,2>>(layout2d.data(), layout2d.size()),
        st.view()
    ).wait();
    std::cout << "Computing Projections.." << std::endl;
    spherical_projection<real_t, index_t>(ctx, dst_view,
        std::span<std::array<real_t,2>>(layout2d.data(), layout2d.size()),
        std::span<std::array<real_t,3>>(xyz     .data(), xyz     .size()),
        st.view()
    ).wait();
    std::cout << "Optimizing Forcefield.." << std::endl;
    forcefield_optimize<ForcefieldType::PEDERSEN, real_t, index_t>(ctx, dst_view,
        std::span<std::array<real_t,3>>(xyz.data(), xyz.size()),
        st.view(), /*batch_iters=*/5*N, /*max_iters=*/5*N
    ).wait();
    ctx.wait();

    int nan_count = 0;
    for (size_t i = 0; i < xyz.size(); ++i) {
        const auto& c = xyz[i];
        if (std::isnan(c[0]) || std::isnan(c[1]) || std::isnan(c[2])) ++nan_count;
    }
    std::cout << nan_count << " NaNs found" << std::endl;

    std::cout << "Writing coordinates to binary file: " << output_file << std::endl;
    std::ofstream out(output_file, std::ios::binary);
    if (!out) {
        std::cerr << "Error opening output file: " << output_file << std::endl;
        return 1;
    }
    out.write(reinterpret_cast<const char*>(xyz.data()),
              xyz.size() * sizeof(std::array<real_t,3>));
    out.close();

    return 0;
}
