#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include <fullerenes/buckygen-wrapper.hh>
#include <fullerenes/triangulation.hh>
#include <fullerenes/graph.hh>
#include <fullerenes/kernel-headers/all-kernels.hh>
#include <fullerenes/planargraph.hh>
#include <fullerenes/polyhedron.hh>

#include <fullerenes/batch/batch.hh>
#include <fullerenes/batch/batch_state.hh>
#include <fullerenes/sycl-headers/sycl-device-queue.hh>
#include <fullerenes/sycl-headers/sycl-vector.hh>

template<typename T>
T stddev(const std::vector<T>& data)
{
    if (data.size() <= 1) return 0;
    T mn = mean(data);
    T sum_of_squares = 0.0;
    for (const T& value : data) {
        T diff = value - mn;
        sum_of_squares += diff * diff;
    }
    T variance = sum_of_squares / (data.size() - 1);
    return std::sqrt(variance);
}

int main(int argc, char** argv) {
    using real_t = float;
    using node_t = uint16_t;
    using RSR    = Spanify::RSRAdjacencyView<node_t>;

    size_t N  = argc > 1 ? strtol(argv[1], 0, 0) : 80;
    size_t Nf = N/2 + 2;
    std::string device_type = argc > 2 ? argv[2] : "gpu";
    size_t Nruns = argc > 3 ? strtol(argv[3], 0, 0) : 10;

    size_t BatchSize = 1280;

    SyclQueue Q{Device(device_type)};

    // ---- Per-run buffers (layout/xyz scratch reused across runs) ------
    SyclVector<std::array<real_t,2>> layout2d(BatchSize * N);
    SyclVector<std::array<real_t,3>> xyz     (BatchSize * N);

    std::span<std::array<real_t,2>> layout_span(layout2d.data(), layout2d.size());
    std::span<std::array<real_t,3>> xyz_span   (xyz.data(), xyz.size());

    Graph G(Nf, true);

    auto fill_and_dualize_host = [&](batch::Batch<RSR>& dst_cubic,
                                     batch::BatchState& st,
                                     double& filltime, double& dualtime)
    {
        auto spans = dst_cubic.view_capacity().spans();
        auto& dst_adj = std::get<0>(spans);
        auto& dst_deg = std::get<1>(spans);

        BuckyGen::buckygen_queue BuckyQ = BuckyGen::start(N, 0, 0);
        double ftime = 0, dtime = 0;
        st.clear();
        for (size_t ii = 0; ii < BatchSize; ++ii) {
            auto start = std::chrono::steady_clock::now();
            auto more = BuckyGen::next_fullerene(BuckyQ, G);
            FullereneDual FD(G);
            auto T0 = std::chrono::steady_clock::now(); ftime += std::chrono::duration<double, std::nano>(T0 - start).count();
            PlanarGraph pG = FD.dual_graph();
            // Copy the cubic adjacency directly into dst_cubic (skip GPU dualize).
            for (size_t u = 0; u < N; ++u) {
                for (size_t k = 0; k < 3; ++k)
                    dst_adj[ii*N*3 + u*3 + k] = node_t(pG.neighbours[u*pG.dmax + k]);
                dst_deg[ii*N + u] = uint8_t(pG.deg[u]);
            }
            st.push_back(uint64_t(ii),
                         StatusFlag(int(StatusEnum::FULLERENEGRAPH_PREPARED)));
            auto T1 = std::chrono::steady_clock::now(); dtime += std::chrono::duration<double, std::nano>(T1 - T0).count();
            if (!more) break;
        }
        BuckyGen::stop(BuckyQ);
        filltime = ftime;
        dualtime = dtime;
    };

    std::vector<double> times_generate(Nruns), times_memcpy(Nruns);
    std::vector<double> times_dual(Nruns), times_tutte(Nruns);
    std::vector<double> times_project(Nruns), times_opt(Nruns);

    TutteFunctor<real_t, node_t>                        tutte_layout;
    SphericalProjectionFunctor<real_t, node_t>          spherical_projection;
    ForcefieldOptimizeFunctor<PEDERSEN, real_t, node_t> forcefield_optimize;

    for (int i = 0; i < (int)Nruns; ++i) {
        batch::Batch<RSR> dst_cubic(N, BatchSize, /*dmax*/3);
        batch::BatchState st(BatchSize);

        fill_and_dualize_host(dst_cubic, st, times_generate[i], times_dual[i]);
        auto dst_view = dst_cubic.view_capacity();

        auto T0 = std::chrono::steady_clock::now();
        auto T1 = std::chrono::steady_clock::now(); times_memcpy[i] = std::chrono::duration<double, std::nano>(T1 - T0).count();
        tutte_layout.compute(Q, dst_view, layout_span, st.view()).wait();
        auto T2 = std::chrono::steady_clock::now(); times_tutte[i] = std::chrono::duration<double, std::nano>(T2 - T1).count();
        spherical_projection.compute(Q, dst_view, layout_span, xyz_span, st.view()).wait();
        auto T3 = std::chrono::steady_clock::now(); times_project[i] = std::chrono::duration<double, std::nano>(T3 - T2).count();
        forcefield_optimize.compute(Q, dst_view, xyz_span, st.view(), 5*N, 5*N).wait();
        auto T4 = std::chrono::steady_clock::now(); times_opt[i] = std::chrono::duration<double, std::nano>(T4 - T3).count();
    }

    std::cout << "N, Nf, BatchSize, device_type, generate, memcpy, dual, tutte, project, opt" << std::endl;
    std::cout << "N: " << N << ", Nf: " << Nf << ", BatchSize: " << BatchSize << ", device_type: " << device_type << "\n";
    std::cout << "Generate: " << mean(times_generate)/BatchSize << ", " << stddev(times_generate)/BatchSize << " ns \n";
    std::cout << "Memcpy: "   << mean(times_memcpy)  /BatchSize << ", " << stddev(times_memcpy)  /BatchSize << " ns \n";
    std::cout << "Dual: "     << mean(times_dual)    /BatchSize << ", " << stddev(times_dual)    /BatchSize << " ns \n";
    std::cout << "Tutte: "    << mean(times_tutte)   /BatchSize << ", " << stddev(times_tutte)   /BatchSize << " ns \n";
    std::cout << "Project: "  << mean(times_project) /BatchSize << ", " << stddev(times_project) /BatchSize << " ns \n";
    std::cout << "Opt: "      << mean(times_opt)     /BatchSize << ", " << stddev(times_opt)     /BatchSize << " ns \n";
    return 0;
}
