#include <array>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <fullerenes/buckygen-wrapper.hh>
#include <fullerenes/dense_graph.hh>
#include <fullerenes/fullerenegraph.hh>
#include <fullerenes/graph.hh>
#include <fullerenes/planargraph.hh>
#include <fullerenes/polyhedron.hh>

#include <fullerenes/batch/batch.hh>
#include <fullerenes/batch/batch_state.hh>
#include <fullerenes/batch/utilities.hh>
#include <fullerenes/kernel-headers/dualize.hh>
#include <fullerenes/kernel-headers/forcefield-optimize.hh>
#include <fullerenes/kernel-headers/spherical-projection.hh>
#include <fullerenes/kernel-headers/tutte.hh>
#include <fullerenes/sycl-headers/sycl-device-queue.hh>
#include <fullerenes/sycl-headers/sycl-vector.hh>

namespace {

using real_t  = double;
using index_t = std::uint16_t;

constexpr int kBatchSize = 2;

std::string output_name(const std::string& prefix, int n, int index)
{
    return prefix + "-C" + std::to_string(n) + "-iso" + std::to_string(index) + ".mol2";
}

} // namespace

int main(int argc, char** argv)
{
    try {
        const int n        = (argc > 1) ? std::stoi(argv[1]) : 60;
        const std::string prefix = (argc > 2) ? argv[2] : "functor-test";
        const int Nf       = n / 2 + 2;
        const int capacity = kBatchSize;

        using RSR = Spanify::RSRAdjacencyView<index_t>;

        // ---- Allocate per-pipeline buffers -----------------------------------
        batch::Batch<RSR> src_dual (Nf, capacity, /*dmax*/6);
        batch::Batch<RSR> dst_cubic(n,  capacity, /*dmax*/3);
        batch::BatchState st(capacity);

        SyclVector<std::array<index_t,6>> faces_cubic_buf(capacity * Nf);
        SyclVector<std::array<index_t,3>> faces_dual_buf (capacity * n);
        SyclVector<std::array<real_t,2>>  layout2d       (capacity * n);
        SyclVector<std::array<real_t,3>>  xyz            (capacity * n);

        // ---- Fill src_dual from BuckyGen -------------------------------------
        auto BQ = BuckyGen::start(n, false, false);
        Graph G(Nf, std::vector<node_t>(6, -1));
        for (int i = 0; i < capacity; ++i) {
            if (!BuckyGen::next_fullerene(BQ, G))
                throw std::runtime_error("BuckyGen exhausted before filling batch");
            batch::push_back(src_dual, G);
            st.push_back(uint64_t(i),
                         StatusFlag(int(StatusEnum::DUAL_INITIALIZED)));
        }
        BuckyGen::stop(BQ);

        // ---- Pipeline --------------------------------------------------------
        SyclQueue Q(Device::default_device(), true);

        auto src_view = src_dual .view_capacity();
        auto dst_view = dst_cubic.view_capacity();

        dualize<real_t, index_t>(Q, src_view, dst_view, st.view(),
            std::span<std::array<index_t,6>>(faces_cubic_buf.data(), faces_cubic_buf.size()),
            std::span<std::array<index_t,3>>(faces_dual_buf .data(), faces_dual_buf .size())
        ).wait();
        tutte_layout<real_t, index_t>(Q, dst_view,
            std::span<std::array<real_t,2>>(layout2d.data(), layout2d.size()),
            st.view()
        ).wait();
        spherical_projection<real_t, index_t>(Q, dst_view,
            std::span<std::array<real_t,2>>(layout2d.data(), layout2d.size()),
            std::span<std::array<real_t,3>>(xyz     .data(), xyz     .size()),
            st.view()
        ).wait();
        forcefield_optimize<ForcefieldType::PEDERSEN, real_t, index_t>(Q, dst_view,
            std::span<std::array<real_t,3>>(xyz.data(), xyz.size()),
            st.view(), /*batch_iters=*/5*n, /*max_iters=*/5*n
        ).wait();
        Q.wait();

        // ---- Write mol2 per isomer -------------------------------------------
        auto dst_cap = dst_cubic.view_capacity();
        for (int i = 0; i < capacity; ++i) {
            auto view_i = dst_cap[std::size_t(i)];
            Graph cubicG = batch::to_graph(view_i);
            std::vector<coord3d> pts = batch::to_points(
                std::span<const std::array<real_t,3>>(xyz.data() + std::size_t(i) * n, n));
            
            std::cout << cubicG.owned_neighbours << '\n';
            Polyhedron polyhedron(PlanarGraph(cubicG), pts);
            std::cout << PlanarGraph(cubicG).owned_neighbours << '\n';
            std::cout << polyhedron.owned_neighbours << '\n';
            const std::string path = output_name(prefix, n, i);
            if (!Polyhedron::to_file(polyhedron, path))
                throw std::runtime_error("failed to write mol2 file: " + path);
            std::cout << "wrote " << path << '\n';
        }

        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "fatal: " << ex.what() << '\n';
        return 1;
    }
}
