#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
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

    size_t N                 = argc > 1 ? stoi(argv[1]) : 200;
    size_t N_TASKS_MAX       = argc > 2 ? stoi(argv[2]) : 8;
    size_t workers_per_task  = argc > 3 ? stoi(argv[3]) : 3;
    size_t chunks_per_worker = argc > 4 ? stoi(argv[4]) : 1;

    size_t NisomersInIsomerspace = IsomerDB::number_isomers(N);
    const char *N_TASKS_str    = getenv("N_TASKS");
    const char *MY_TASK_ID_str = getenv("MY_TASK_ID");
    size_t N_TASKS    = N_TASKS_str    ? std::stoi(N_TASKS_str)    : 1;
    size_t MY_TASK_ID = MY_TASK_ID_str ? std::stoi(MY_TASK_ID_str) : 0;

    auto IsomerPerNodeEstimate = NisomersInIsomerspace / N_TASKS;
    auto BatchSize = std::min<int>(IsomerPerNodeEstimate, (1<<17));

    size_t Nf = N/2 + 2;
    std::string filename = "output/full_pipeline_" + std::string(getenv("SLURM_JOB_ID"))
        + "_" + std::to_string(MY_TASK_ID) + "_" + std::to_string(N_TASKS) + ".csv";
    std::ofstream myfile(filename);

    SyclQueue Q("gpu");

    size_t N_chunks = N_TASKS_MAX * workers_per_task * chunks_per_worker;
    size_t n_chunks = N_chunks / N_TASKS;
    auto my_chunks = loadbalanced_chunks(N_chunks, n_chunks, MY_TASK_ID);
    BuckyGen::buckyherd_queue BuckyQ(N, N_chunks, workers_per_task,
                                     false, false, my_chunks);

    size_t isomers_in_queue = 0;
    Graph G(N);
    bool more = true;

    double times_generate = 0., times_memcpy = 0., times_dual = 0.;
    double times_tutte = 0., times_project = 0., times_opt = 0.;

    // ---- Batch containers ----------------------------------------------
    batch::Batch<RSR> src_dual (Nf, BatchSize, /*dmax*/6);
    batch::Batch<RSR> dst_cubic(N,  BatchSize, /*dmax*/3);
    batch::BatchState st(BatchSize);

    SyclVector<std::array<node_t,6>> faces_cubic_buf(BatchSize * Nf);
    SyclVector<std::array<node_t,3>> faces_dual_buf (BatchSize * N);
    SyclVector<std::array<real_t,2>> layout2d       (BatchSize * N);
    SyclVector<std::array<real_t,3>> xyz            (BatchSize * N);

    DualizeFunctor<real_t, node_t>                       dualize_V1;
    TutteFunctor<real_t, node_t>                         tutte_layout;
    SphericalProjectionFunctor<real_t, node_t>           spherical_projection;
    ForcefieldOptimizeFunctor<PEDERSEN, real_t, node_t>  forcefield_optimize;

    std::span<std::array<real_t,3>> xyz_span(xyz.data(), xyz.size());
    std::span<std::array<real_t,2>> layout_span(layout2d.data(), layout2d.size());
    std::span<std::array<node_t,6>> faces_cubic_span(faces_cubic_buf.data(), faces_cubic_buf.size());
    std::span<std::array<node_t,3>> faces_dual_span (faces_dual_buf .data(), faces_dual_buf .size());

    auto generate_and_fill = [&]() {
        auto spans = src_dual.view_capacity().spans();
        auto& src_adj = std::get<0>(spans);
        auto& src_deg = std::get<1>(spans);
        st.clear();
        int isomer_idx = 0;
        while (more && isomer_idx < BatchSize) {
            more &= BuckyQ.next_fullerene(G);
            if (!more) break;
            for (size_t j = 0; j < Nf; ++j) {
                const auto deg_j = G.degree(j);
                for (size_t k = 0; k < deg_j; ++k)
                    src_adj[isomer_idx*Nf*6 + j*6 + k] = node_t(G.nbrs(j)[k]);
                if (deg_j == 5) {
                    src_adj[isomer_idx*Nf*6 + j*6 + 5] = std::numeric_limits<node_t>::max();
                    src_deg[isomer_idx*Nf + j] = 5;
                } else {
                    src_deg[isomer_idx*Nf + j] = 6;
                }
            }
            st.push_back(uint64_t(isomers_in_queue),
                         StatusFlag(int(StatusEnum::DUAL_INITIALIZED)));
            ++isomer_idx;
            ++isomers_in_queue;
        }
    };

    auto src_view = src_dual.view_capacity();
    auto dst_view = dst_cubic.view_capacity();

    while (more) {
        auto T1 = std::chrono::steady_clock::now();
        generate_and_fill();
        auto T2 = std::chrono::steady_clock::now(); times_generate = std::chrono::duration<double, std::nano>(T2 - T1).count();
        auto T3 = std::chrono::steady_clock::now(); times_memcpy   = std::chrono::duration<double, std::nano>(T3 - T2).count();
        dualize_V1.compute(Q, src_view, dst_view, st.view(), faces_cubic_span, faces_dual_span).wait();
        auto T4 = std::chrono::steady_clock::now(); times_dual     = std::chrono::duration<double, std::nano>(T4 - T3).count();
        tutte_layout.compute(Q, dst_view, layout_span, st.view()).wait();
        auto T5 = std::chrono::steady_clock::now(); times_tutte    = std::chrono::duration<double, std::nano>(T5 - T4).count();
        spherical_projection.compute(Q, dst_view, layout_span, xyz_span, st.view()).wait();
        auto T6 = std::chrono::steady_clock::now(); times_project  = std::chrono::duration<double, std::nano>(T6 - T5).count();
        forcefield_optimize.compute(Q, dst_view, xyz_span, st.view(), 5*N, 5*N).wait();
        auto T7 = std::chrono::steady_clock::now(); times_opt      = std::chrono::duration<double, std::nano>(T7 - T6).count();
    }

    myfile << "N, Nf, BatchSize, JOBID, NTASKS, TASK_ID, FILL_ME_UP_SCOTTY, MEMCPY, DUAL, TUTTE, PROJECT, OPT\n"
           << N << ", " << Nf << ", " << isomers_in_queue << ", "
           << getenv("SLURM_JOB_ID") << ", " << N_TASKS << ", " << MY_TASK_ID << ", "
           << times_generate/std::max<size_t>(1,isomers_in_queue) << ", "
           << times_memcpy  /std::max<size_t>(1,isomers_in_queue) << ", "
           << times_dual    /std::max<size_t>(1,isomers_in_queue) << ", "
           << times_tutte   /std::max<size_t>(1,isomers_in_queue) << ", "
           << times_project /std::max<size_t>(1,isomers_in_queue) << ", "
           << times_opt     /std::max<size_t>(1,isomers_in_queue) << "\n";
    return 0;
}
