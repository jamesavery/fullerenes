#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include <fullerenes/buckygen-wrapper.hh>
#include <fullerenes/graph.hh>
#include <fullerenes/isomerdb.hh>
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

vector<size_t> loadbalanced_chunks(size_t N_chunks, size_t n_chunks, size_t my_node_idx = 0, int /*seed*/ = 42)
{
    vector<size_t> all_chunks(N_chunks), my_chunks(n_chunks);
    for (size_t i = 0; i < N_chunks; i++) all_chunks[i] = i;
    auto rng = std::default_random_engine(42);
    shuffle(all_chunks.begin(), all_chunks.end(), rng);
    for (size_t i = 0; i < n_chunks; i++) my_chunks[i] = all_chunks[my_node_idx * n_chunks + i];
    return my_chunks;
}


int main(int argc, char** argv) {
    using real_t = float;
    using node_t = uint16_t;
    using RSR    = Spanify::RSRAdjacencyView<node_t>;

    size_t N  = argc > 1 ? strtol(argv[1], 0, 0) : 200;
    size_t NisomersInIsomerspace = IsomerDB::number_isomers(N);
    auto   BatchSize = std::min((int)NisomersInIsomerspace, (int)1280 * 10);

    size_t Nf = N / 2 + 2;

    SyclQueue Q("gpu");

    double times_generate = 0., times_memcpy = 0., times_dual = 0., times_tutte = 0.;
    double times_project = 0., times_opt = 0., times_hessian = 0.;
    double times_spectral_ends = 0., times_spectral = 0., times_spectral_vectors = 0.;

    // ---- Batch containers ----------------------------------------------
    batch::Batch<RSR> src_dual (Nf, BatchSize, /*dmax*/6);
    batch::Batch<RSR> dst_cubic(N,  BatchSize, /*dmax*/3);
    batch::BatchState st(BatchSize);

    SyclVector<std::array<node_t,6>> faces_cubic_buf(BatchSize * Nf);
    SyclVector<std::array<node_t,3>> faces_dual_buf (BatchSize * N);
    SyclVector<std::array<real_t,2>> layout2d       (BatchSize * N);
    SyclVector<std::array<real_t,3>> xyz            (BatchSize * N);

    SyclVector<real_t>    hessian_buffer         (BatchSize * N * 90);
    SyclVector<node_t>    cols_buffer            (BatchSize * N * 90);
    SyclVector<real_t>    spectral_ends_buffer   (BatchSize * 2);
    SyclVector<real_t>    spectral_buffer        (BatchSize * N * 3);
    SyclVector<real_t>    spectral_vectors_buffer(BatchSize * N * 3 * N * 3);

    const size_t n_lanczos = 40;
    SyclVector<real_t>    off_diagonal(BatchSize * n_lanczos);
    SyclVector<real_t>    qmat        (BatchSize * n_lanczos * n_lanczos);
    SyclVector<real_t>    lanczos_buf (BatchSize * n_lanczos * N * 3);
    SyclVector<real_t>    diag_buf    (BatchSize * n_lanczos);
    SyclVector<node_t>    ends_idx    (BatchSize * 2);

    // ---- Kernel functors ------------------------------------------------
    DualizeFunctor<real_t, node_t>                                      dualize_V1;
    TutteFunctor<real_t, node_t>                                        tutte_layout;
    SphericalProjectionFunctor<real_t, node_t>                          spherical_projection;
    ForcefieldOptimizeFunctor<PEDERSEN, real_t, node_t>                 forcefield_optimize;
    HessianFunctor<PEDERSEN, real_t, node_t>                            compute_hessians;
    EigenFunctor<EigensolveMode::ENDS, real_t, node_t>                  eigensolve_ends;
    EigenFunctor<EigensolveMode::FULL_SPECTRUM, real_t, node_t>         eigensolve_full;

    // ---- Fill src_dual (and state) from BuckyGen -----------------------
    {
        auto src_spans = src_dual.view_capacity().spans();
        auto& src_adj = std::get<0>(src_spans);
        auto& src_deg = std::get<1>(src_spans);
        auto BQ = BuckyGen::start(N, false, false);
        Graph G(Nf, std::vector<int>(6, -1));
        int generated = 0;
        for (int i = 0; i < BatchSize; ++i) {
            if (!BuckyGen::next_fullerene(BQ, G)) break;
            ++generated;
            for (int u = 0; u < (int)Nf; ++u) {
                for (int k = 0; k < 6; ++k)
                    src_adj[i*Nf*6 + u*6 + k] = node_t(G.neighbours[u*G.dmax + k]);
                src_deg[i*Nf + u] = uint8_t(G.deg[u]);
            }
            st.push_back(uint64_t(i),
                         StatusFlag(int(StatusEnum::DUAL_INITIALIZED)));
        }
        // Duplicate cyclically to fill the batch.
        for (int i = generated; i < BatchSize; ++i) {
            int src = i % std::max(1, generated);
            std::copy(src_adj.begin() + src*Nf*6,
                      src_adj.begin() + (src+1)*Nf*6,
                      src_adj.begin() + i*Nf*6);
            std::copy(src_deg.begin() + src*Nf,
                      src_deg.begin() + (src+1)*Nf,
                      src_deg.begin() + i*Nf);
            st.push_back(uint64_t(i),
                         StatusFlag(int(StatusEnum::DUAL_INITIALIZED)));
        }
        BuckyGen::stop(BQ);
    }

    Span<std::array<real_t,3>> xyz_span(xyz.data(), xyz.size());
    Span<std::array<real_t,2>> layout_span(layout2d.data(), layout2d.size());
    Span<std::array<node_t,6>> faces_cubic_span(faces_cubic_buf.data(), faces_cubic_buf.size());
    Span<std::array<node_t,3>> faces_dual_span (faces_dual_buf .data(), faces_dual_buf .size());

    auto src_view = src_dual .view_capacity();
    auto dst_view = dst_cubic.view_capacity();

    const auto Nruns = 10;
    for (size_t i = 0; i < Nruns; i++) {
        auto T1 = std::chrono::steady_clock::now(); times_generate += std::chrono::duration<double, std::nano>(T1 - T1).count();
        auto T2 = std::chrono::steady_clock::now(); times_generate += std::chrono::duration<double, std::nano>(T2 - T1).count();
        auto T3 = std::chrono::steady_clock::now(); times_memcpy += std::chrono::duration<double, std::nano>(T3 - T2).count();
        dualize_V1.compute(Q, src_view, dst_view, st.view(),
                           faces_cubic_span, faces_dual_span).wait();
        auto T4 = std::chrono::steady_clock::now(); times_dual += std::chrono::duration<double, std::nano>(T4 - T3).count();
        tutte_layout.compute(Q, dst_view, layout_span, st.view()).wait();
        auto T5 = std::chrono::steady_clock::now(); times_tutte += std::chrono::duration<double, std::nano>(T5 - T4).count();
        spherical_projection.compute(Q, dst_view, layout_span, xyz_span, st.view()).wait();
        auto T6 = std::chrono::steady_clock::now(); times_project += std::chrono::duration<double, std::nano>(T6 - T5).count();
        forcefield_optimize.compute(Q, dst_view, xyz_span, st.view(), 4*N, 4*N).wait();
        auto T7 = std::chrono::steady_clock::now(); times_opt += std::chrono::duration<double, std::nano>(T7 - T6).count();
        compute_hessians.compute(Q, dst_view, xyz_span, st.view(),
                                 hessian_buffer, cols_buffer).wait();
        auto T8 = std::chrono::steady_clock::now(); times_hessian += std::chrono::duration<double, std::nano>(T8 - T7).count();
        eigensolve_ends.compute(Q, xyz_span, (int)N, BatchSize,
            hessian_buffer, cols_buffer, n_lanczos,
            spectral_ends_buffer, /*eigenvectors*/Span<real_t>(),
            off_diagonal, qmat, lanczos_buf, diag_buf, ends_idx).wait();
        auto T9 = std::chrono::steady_clock::now(); times_spectral_ends += std::chrono::duration<double, std::nano>(T9 - T8).count();
        eigensolve_full.compute(Q, xyz_span, (int)N, BatchSize,
            hessian_buffer, cols_buffer, n_lanczos,
            spectral_buffer, /*eigenvectors*/Span<real_t>(),
            off_diagonal, qmat, lanczos_buf, diag_buf, ends_idx).wait();
        auto T10 = std::chrono::steady_clock::now(); times_spectral += std::chrono::duration<double, std::nano>(T10 - T9).count();
        auto T11 = std::chrono::steady_clock::now(); times_spectral_vectors += std::chrono::duration<double, std::nano>(T11 - T10).count();
    }

    std::cout << "N, Nf, BatchSize, JOBID, NTASKS, TASK_ID, FILL_ME_UP_SCOTTY, MEMCPY, DUAL, TUTTE, PROJECT, OPT\n"
              << N << ", " << Nf << ", " << (BatchSize*Nruns) << ", "
              << times_generate/(BatchSize*Nruns) << ", "
              << times_memcpy/(BatchSize*Nruns) << ", "
              << times_dual/(BatchSize*Nruns) << ", "
              << times_tutte/(BatchSize*Nruns) << ", "
              << times_project/(BatchSize*Nruns) << ", "
              << times_opt/(BatchSize*Nruns) << ", "
              << times_hessian/(BatchSize*Nruns) << ", "
              << times_spectral_ends/(BatchSize*Nruns) << ", "
              << times_spectral/(BatchSize*Nruns) << ", "
              << times_spectral_vectors/(BatchSize*Nruns) << std::endl;
    return 0;
}
