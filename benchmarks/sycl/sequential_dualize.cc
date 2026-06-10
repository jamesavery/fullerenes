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

    size_t N        = argc > 1 ? strtol(argv[1], 0, 0) : 80;
    size_t Nf       = N/2 + 2;
    size_t NumNodes = argc > 2 ? strtol(argv[2], 0, 0) : 2000000;
    std::string device_type = argc > 3 ? argv[3] : "gpu";
    size_t Nruns    = argc > 4 ? strtol(argv[4], 0, 0) : 10;

    size_t BatchSize = std::ceil((real_t)NumNodes / (real_t)N);

    SyclQueue Q(device_type);
    (void)Q;

    batch::Batch<RSR> src_dual (Nf, BatchSize, /*dmax*/6);
    batch::Batch<RSR> dst_cubic(N,  BatchSize, /*dmax*/3);
    batch::BatchState st(BatchSize);

    Graph G(N);
    auto fill_and_dualize = [&](double& filltime, double& dualtime) {
        auto src_spans = src_dual.view_capacity().spans();
        auto& src_adj  = std::get<0>(src_spans);
        auto& src_deg  = std::get<1>(src_spans);
        auto dst_spans = dst_cubic.view_capacity().spans();
        auto& dst_adj  = std::get<0>(dst_spans);
        auto& dst_deg  = std::get<1>(dst_spans);

        st.clear();
        double ftime = 0, dtime = 0;
        BuckyGen::buckygen_queue BuckyQ = BuckyGen::start(N, 0, 0);
        for (size_t ii = 0; ii < BatchSize; ++ii) {
            auto start = std::chrono::steady_clock::now();
            auto more = BuckyGen::next_fullerene(BuckyQ, G);
            for (size_t j = 0; j < Nf; ++j) {
                const auto deg_j = G.degree(j);
                for (size_t k = 0; k < deg_j; ++k)
                    src_adj[ii*Nf*6 + j*6 + k] = node_t(G.nbrs(j)[k]);
                if (deg_j == 5) {
                    src_adj[ii*Nf*6 + j*6 + 5] = std::numeric_limits<node_t>::max();
                    src_deg[ii*Nf + j] = 5;
                } else {
                    src_deg[ii*Nf + j] = 6;
                }
            }
            FullereneDual FD(G);
            auto T0 = std::chrono::steady_clock::now(); ftime += std::chrono::duration<double, std::nano>(T0 - start).count();
            PlanarGraph pG = FD.dual_graph();
            for (size_t j = 0; j < N; ++j) {
                for (size_t k = 0; k < 3; ++k)
                    dst_adj[ii*N*3 + j*3 + k] = node_t(pG.nbrs(j)[k]);
                dst_deg[ii*N + j] = 3;
            }
            auto T1 = std::chrono::steady_clock::now(); dtime += std::chrono::duration<double, std::nano>(T1 - T0).count();
            if (!more) break;
            st.push_back(uint64_t(ii),
                         StatusFlag(int(StatusEnum::DUAL_INITIALIZED)));
        }
        BuckyGen::stop(BuckyQ);
        filltime = ftime;
        dualtime = dtime;
    };

    std::vector<double> times_dual(Nruns);
    double not_used = 0;

    for (int i = 0; i < (int)Nruns; ++i) {
        fill_and_dualize(not_used, times_dual[i]);
    }

    std::cout << "N, Nf, BatchSize, device_type \n";
    std::cout << "N: " << N << ", Nf: " << Nf << ", BatchSize: " << BatchSize
              << ", device_type: " << device_type << "\n";
    std::cout << "Dual: " << mean(times_dual)/BatchSize << ", "
              << stddev(times_dual)/BatchSize << " ns \n";
    return 0;
}
