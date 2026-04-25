#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include <fullerenes/buckygen-wrapper.hh>
#include <fullerenes/graph.hh>
#include <fullerenes/kernel-headers/all-kernels.hh>

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
    size_t Nwarmup  = argc > 5 ? strtol(argv[5], 0, 0) : 1;
    auto dual_version = argc > 6 ? strtol(argv[6], 0, 0) : 0;
    (void)dual_version;

    size_t BatchSize = std::ceil((real_t)NumNodes / (real_t)N);

    SyclQueue Q(device_type);

    batch::Batch<RSR> src_dual (Nf, BatchSize, /*dmax*/6);
    batch::Batch<RSR> dst_cubic(N,  BatchSize, /*dmax*/3);
    batch::BatchState st(BatchSize);

    SyclVector<std::array<node_t,6>> faces_cubic_buf(BatchSize * Nf);
    SyclVector<std::array<node_t,3>> faces_dual_buf (BatchSize * N);
    std::span<std::array<node_t,6>> faces_cubic_span(faces_cubic_buf.data(), faces_cubic_buf.size());
    std::span<std::array<node_t,3>> faces_dual_span (faces_dual_buf.data(),  faces_dual_buf.size());

    Graph G(N);
    auto fill = [&]() {
        auto spans = src_dual.view_capacity().spans();
        auto& src_adj = std::get<0>(spans);
        auto& src_deg = std::get<1>(spans);
        st.clear();

        BuckyGen::buckygen_queue BuckyQ = BuckyGen::start(N, 0, 0);
        for (size_t ii = 0; ii < BatchSize; ++ii) {
            auto more = BuckyGen::next_fullerene(BuckyQ, G);
            if (!more) break;
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
            st.push_back(uint64_t(ii),
                         StatusFlag(int(StatusEnum::DUAL_INITIALIZED)));
        }
        BuckyGen::stop(BuckyQ);
    };

    std::vector<double> times_generate(Nruns + Nwarmup);
    std::vector<double> times_memcpy  (Nruns + Nwarmup);
    std::vector<double> times_dual    (Nruns);

    DualizeFunctor<real_t, node_t> dualize;

    auto src_view = src_dual.view_capacity();
    auto dst_view = dst_cubic.view_capacity();

    for (int i = 0; i < (int)(Nruns + Nwarmup); ++i) {
        auto start = std::chrono::steady_clock::now();
        fill();
        auto T0 = std::chrono::steady_clock::now();
        times_generate[i] = std::chrono::duration<double, std::nano>(T0 - start).count();
        auto T1 = std::chrono::steady_clock::now();
        times_memcpy[i]   = std::chrono::duration<double, std::nano>(T1 - T0).count();
        dualize.compute(Q, src_view, dst_view, st.view(), faces_cubic_span, faces_dual_span).wait();
        if (i >= (int)Nwarmup) {
            auto T2 = std::chrono::steady_clock::now();
            times_dual[i - Nwarmup] = std::chrono::duration<double, std::nano>(T2 - T1).count();
        }
    }

    std::cout << "N, Nf, BatchSize, device_type" << std::endl;
    std::cout << "N: " << N << ", Nf: " << Nf << ", BatchSize: " << BatchSize
              << ", device_type: " << device_type << "\n";
    std::cout << "Memory copy: " << mean(times_memcpy)/BatchSize << ", "
              << stddev(times_memcpy)/BatchSize << " ns \n";
    std::cout << "Dual: " << mean(times_dual)/BatchSize << ", "
              << stddev(times_dual)/BatchSize << " ns \n";
    return 0;
}
