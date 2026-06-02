// isomerspace-experiment.cc — Phase 10 migration to BatchQueue + utilities.
//
// Multi-device pipeline (GEN -> GEO -> OPT -> PROP -> STAT -> postprocess) using
// `batch::Batch<ADJ>` (dual-graph) and `batch::Batch<PV>` (cubic graph + xyz)
// for stage buffers, and `batch::BatchQueue<ADJ>` / `batch::BatchQueue<PV>`
// to bridge stages. Per-slot xyz lives inside the polyhedron view itself
// (RSRPolyhedronView's `points` field) so a single BatchQueue carries
// graph + geometry atomically.

#include <csignal>
#include <cstring>
#include <sys/stat.h>
#include <limits.h>
#include <cmath>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <thread>
#include <numeric>
#include <future>
#include <queue>

#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/progress_bar.hh"
#include <fullerenes/kernel-headers/all-kernels.hh>
#include "fullerenes/isomerdb.hh"

#include <fullerenes/batch/batch.hh>
#include <fullerenes/batch/batch_state.hh>
#include <fullerenes/batch/batch_queue.hh>
#include <fullerenes/batch/utilities.hh>
#include <fullerenes/sycl-headers/sycl-device-queue.hh>
#include <fullerenes/sycl-headers/sycl-vector.hh>

#define STORE_HESSIAN 1
#define STORE_EIGENSYSTEM 1

using namespace std;
using namespace std::chrono;

typedef float    real_t;
typedef uint16_t device_node_t;
using ADJ = Spanify::RSRAdjacencyView<device_node_t>;
using PV  = Spanify::RSRPolyhedronView<real_t, device_node_t>;

typedef enum { GEN, GEO, OPT, PROP, STAT, NUM_STAGES } stage_t;
typedef enum { VOLUME, ECCENTRICITY, MIN_FREQ, MAX_FREQ, FREQ_WIDTH, INERTIA, NUM_RESULTS } results_t;
string result_names[NUM_RESULTS] = {"volume", "eccentricity", "minfreq", "maxfreq", "freqwidth", "inertia"};
string stage_names[NUM_STAGES]   = {"generate", "start_geometry", "optimized_geometry", "properties", "statistics"};
constexpr int result_sizes[NUM_RESULTS] = {1, 1, 1, 1, 1, 3};

void signal_callback_handler(int signum) { exit(signum); }

#include <sys/resource.h>
void ensure_minimal_stacksize(rlim_t minimal_stacksize)
{
  struct rlimit rl;
  int error = getrlimit(RLIMIT_STACK, &rl);
  rl.rlim_cur = std::max(rl.rlim_cur, minimal_stacksize);
  error += setrlimit(RLIMIT_STACK, &rl);
  if (error < 0) { perror("ensure_minimal_stacksize failed: "); abort(); }
}

struct isomer_candidate
{
  float value;
  int   id;
  array<float, NUM_RESULTS> results;

  vector<uint16_t> cubic_neighbours;
  vector<real_t>   X;

  // Construct from per-isomer cubic adjacency span + xyz span (flattened).
  isomer_candidate(double v, int id_, array<float, NUM_RESULTS>& r,
                   const device_node_t* adj, const std::array<real_t,3>* xyz, size_t N)
      : value(v), id(id_), results(r), cubic_neighbours(3 * N), X(3 * N)
  {
    // Adjacency stride is 3 (cubic graph stored as N*3).
    std::memcpy(cubic_neighbours.data(), adj, 3 * N * sizeof(device_node_t));
    for (size_t u = 0; u < N; ++u) {
      X[3*u + 0] = xyz[u][0];
      X[3*u + 1] = xyz[u][1];
      X[3*u + 2] = xyz[u][2];
    }
  }

  bool operator<(const isomer_candidate& b) const { return value < b.value; }
  friend ostream& operator<<(ostream& s, const isomer_candidate& c) {
    s << "(" << c.value << "," << c.id << ")"; return s;
  }
};

template <typename T>
struct k_smallest : public std::priority_queue<T> {
  size_t k;
  k_smallest(size_t k) : k(k) {}
  bool insert(const T& x) {
    if (this->size() < k || x < this->top()) {
      if (this->size() == k) this->pop();
      this->push(x);
      return true;
    }
    return false;
  }
  vector<T>& as_vector() { return (*this).*&k_smallest::c; }
};

// -----------------------------------------------------------------------------
// Stage-buffer aggregate. Each stage's working batch carries a dual-graph
// batch (`src_dual`, ADJ) and a cubic-graph + xyz batch (`dst_cubic`, PV)
// plus per-slot scratch (layout2d, faces) and a BatchState. The xyz array
// lives inside `dst_cubic` (PV adds `points` to ADJ).
// -----------------------------------------------------------------------------
struct StageBatch {
    batch::Batch<ADJ>                       src_dual;
    batch::Batch<PV>                        dst_cubic;
    batch::BatchState                       st;
    SyclVector<std::array<real_t,2>>        layout2d;
    SyclVector<std::array<device_node_t,6>> faces_cubic;
    SyclVector<std::array<device_node_t,3>> faces_dual;

    StageBatch() = default;
    StageBatch(int N, int Nf, int capacity)
        : src_dual (Nf, capacity, /*dmax*/6),
          dst_cubic(N,  capacity, /*dmax*/3),
          st(capacity),
          layout2d   (size_t(capacity) * size_t(N)),
          faces_cubic(size_t(capacity) * size_t(Nf)),
          faces_dual (size_t(capacity) * size_t(N)) {}

    int N()        const { return dst_cubic.N(); }
    int Nf()       const { return src_dual.N(); }
    int capacity() const { return dst_cubic.capacity(); }

    // The cubic graph viewed as plain adjacency (for kernels expecting ADJ).
    batch::BatchView<ADJ> cubic_adj_view() {
        return Spanify::as_adjacency_view(dst_cubic.view());
    }
    batch::BatchView<ADJ> cubic_adj_view_capacity() {
        return Spanify::as_adjacency_view(dst_cubic.view_capacity());
    }
    // Flat per-capacity xyz span (length capacity*N).
    std::span<std::array<real_t,3>> xyz_span() {
        return Spanify::points_span(dst_cubic.view_capacity());
    }
    std::span<std::array<real_t,2>>          layout_span() {
        return std::span<std::array<real_t,2>>(layout2d.data(), layout2d.size());
    }
    std::span<std::array<device_node_t,6>>   faces_cubic_span() {
        return std::span<std::array<device_node_t,6>>(faces_cubic.data(), faces_cubic.size());
    }
    std::span<std::array<device_node_t,3>>   faces_dual_span() {
        return std::span<std::array<device_node_t,3>>(faces_dual.data(), faces_dual.size());
    }
};

int main(int ac, char** argv)
{
  signal(SIGINT, signal_callback_handler);
  if (ac < 2) {
    fprintf(stderr, "Syntax: %s <N:int> [output_dir] [IPR:0|1] [only_nontrivial:0|1]\n", argv[0]);
    return -1;
  }
  const size_t N  = strtol(argv[1], 0, 0);
  const size_t Nf = N/2 + 2;

  string output_dir       = ac >= 3 ? argv[2] : string("output/C") + to_string(N);
  bool   IPR              = ac >= 4 ? strtol(argv[3], 0, 0) : 0;
  bool   only_nontrivial  = ac >= 5 ? strtol(argv[4], 0, 0) : 0;
  size_t n_best_candidates= ac >= 6 ? strtol(argv[5], 0, 0) : 30;

  BuckyGen::buckygen_queue BuckyQ = BuckyGen::start(N, IPR, only_nontrivial);
  ensure_minimal_stacksize(N * 10000000);

  mkdir("output", 0777);
  size_t n_fullerenes = IsomerDB::number_isomers(N, only_nontrivial ? "Nontrivial" : "Any", IPR);
  n_best_candidates = std::max(size_t(1), std::min(n_best_candidates, n_fullerenes / 4));

  int Nd = Device::get_devices(DeviceType::GPU).size();
  vector<Device> devices = Device::get_devices(DeviceType::GPU);
  size_t batch_size       = std::min((size_t)(devices[0].get_property(DeviceProperty::MAX_COMPUTE_UNITS) * 10),
                                     (size_t)(n_fullerenes / Nd));
  size_t final_batch_size = n_best_candidates * NUM_RESULTS * 2 / Nd + (Nd - 1);

  cout << "Analyzing C" << N << " isomer space with " << n_fullerenes << " isomers.\n";
  cout << "n_best_candidates = " << n_best_candidates
       << "; batch_size = " << batch_size
       << "; final_batch_size = " << final_batch_size << "\n";

  array<array<double, 2>, NUM_RESULTS> result_bounds = {{
      {0.24 * pow(N, 3 / 2.) - 50, 0.4 * pow(N, 3 / 2.) + 50},
      {1, N / 10.},
      {0, 10 * 33.356},
      {(39 + 2 * log(N - 20)) * 33.356, (43 + 2 * log(N - 20)) * 33.356},
      {40 * 33.356, 55 * 33.356},
      {0, 0}
  }};

  array<vector<SyclVector<real_t>>, NUM_RESULTS> results;
  array<long double, NUM_RESULTS> means{}, stdevs{};
  array<vector<size_t>, NUM_RESULTS> terrible_outliers;
  vector<SyclVector<real_t>>          hessians(Nd);
  vector<SyclVector<device_node_t>>   hessian_cols(Nd);
  vector<SyclVector<real_t>>          Q_eig(Nd), lams(Nd);
  vector<SyclVector<real_t>>          intermediate_spectrum_ends(Nd);
  // EigenFunctor view-based scratch.
  vector<SyclVector<real_t>>          ev_off_diag(Nd), ev_qmat(Nd), ev_lanczos(Nd), ev_diag(Nd);
  vector<SyclVector<device_node_t>>   ev_ends_idx(Nd);

  array<long double, NUM_RESULTS> Ev{}, Ev2{}, K_acc{}, current_mean{}, current_stddev{};
  size_t n_converged = 0;
  for (int r = 0; r < NUM_RESULTS; r++) results[r].resize(Nd);
  const size_t n_lanczos = 50;
  for (int d = 0; d < Nd; d++) {
    hessians[d]                = SyclVector<real_t>(N * 3 * 10 * 3 * batch_size, 0);
    hessian_cols[d]            = SyclVector<device_node_t>(N * 3 * 10 * 3 * batch_size, 0);
    intermediate_spectrum_ends[d] = SyclVector<real_t>(2 * batch_size, 0);
    ev_off_diag[d]             = SyclVector<real_t>(n_lanczos * batch_size, 0);
    ev_qmat[d]                 = SyclVector<real_t>(n_lanczos * n_lanczos * batch_size, 0);
    ev_lanczos[d]              = SyclVector<real_t>(n_lanczos * N * 3 * batch_size, 0);
    ev_diag[d]                 = SyclVector<real_t>(n_lanczos * batch_size, 0);
    ev_ends_idx[d]             = SyclVector<device_node_t>(2 * batch_size, 0);
    for (int r = 0; r < NUM_RESULTS; r++)
      results[r][d] = SyclVector<real_t>(result_sizes[r] * batch_size, 0);
  }

  vector<k_smallest<isomer_candidate>>
      result_min(NUM_RESULTS, k_smallest<isomer_candidate>(n_best_candidates)),
      result_max(NUM_RESULTS, k_smallest<isomer_candidate>(n_best_candidates));

  array<set<pair<real_t, int>>, NUM_RESULTS> result_reference;

  // Per-stage, per-device batches.  Keep four stages: GEO, OPT, PROP, plus
  // a "final" batch of size final_batch_size.
  std::array<std::vector<StageBatch>, 4> Bs;
  for (int i = 0; i < 4; i++) {
    Bs[i].resize(Nd);
    for (int j = 0; j < Nd; j++) {
      Bs[i][j] = StageBatch((int)N, (int)Nf, (int)batch_size);
    }
  }
  vector<StageBatch> final_batch;
  final_batch.reserve(Nd);
  for (int d = 0; d < Nd; d++)
    final_batch.emplace_back((int)N, (int)Nf, (int)final_batch_size);

  // Cross-stage queues.
  vector<batch::BatchQueue<ADJ>> Q0s;          // dual-graph queue (Nf, dmax 6)
  vector<batch::BatchQueue<PV>>  Q1s, Q2s, Q3s; // cubic graph + xyz (N, dmax 3)
  Q0s.reserve(Nd); Q1s.reserve(Nd); Q2s.reserve(Nd); Q3s.reserve(Nd);
  for (int d = 0; d < Nd; d++) {
    Q0s.emplace_back((int)Nf, (int)(batch_size * 4), /*dmax*/6);
    Q1s.emplace_back((int)N,  (int)(batch_size * 4), /*dmax*/3);
    Q2s.emplace_back((int)N,  (int)(batch_size * 4), /*dmax*/3);
    Q3s.emplace_back((int)N,  (int)(batch_size * 4), /*dmax*/3);
  }

  // Kernels are now invoked via free functions (Phase 4).

  vector<SyclQueue> gen_ctxs;
  std::generate_n(std::back_inserter(gen_ctxs), Nd, [&, i = 0]() mutable {
    return SyclQueue{devices[i++ % devices.size()], true};
  });
  vector<SyclQueue> opt_ctxs;
  std::generate_n(std::back_inserter(opt_ctxs), Nd, [&, i = 0]() mutable {
    return SyclQueue{devices[i++ % devices.size()], true};
  });

  ProgressBar progress_bar = ProgressBar('#', 30);

  int I = 0;
  vector<bool> stage_finished(NUM_STAGES, 0);
  vector<int>  num_finished(NUM_STAGES, 0);
  vector<vector<int>> num_finished_this_round(NUM_STAGES, vector<int>(Nd));

  Graph G;
  static_cast<Graph&>(G) = Graph(Nf, std::vector<int>(6, -1));
  G.N = Nf;

  constexpr real_t carbon_mass = 1.9944733e-26, aangstrom_length = 1e-10;
  constexpr size_t num_bins = 1000;
  array<array<size_t, num_bins>, NUM_RESULTS> result_histograms{};
  vector<array<size_t, num_bins * num_bins>> result_histograms2D(NUM_RESULTS * NUM_RESULTS, {0});

  auto T0      = steady_clock::now();
  auto Tgen    = steady_clock::now() - T0;
  auto Tupdate = steady_clock::now() - T0;
  auto Tqueue  = steady_clock::now() - T0;
  auto Tdual   = steady_clock::now() - T0;
  auto Ttutte  = steady_clock::now() - T0;
  auto TX0     = steady_clock::now() - T0;
  auto Tcopy   = steady_clock::now() - T0;
  auto Topt    = steady_clock::now() - T0;
  auto Tfile   = steady_clock::now() - T0;
  auto Toutq   = steady_clock::now() - T0;
  auto Tinq    = steady_clock::now() - T0;
  auto Tinit_geom = steady_clock::now() - T0;

  auto loop_iters = 0;
  auto Tstart = steady_clock::now();
  auto n0 = 0;

  auto generate_isomers = [&](int M) {
    bool more_isomers = true;
    if (I == (int)n_fullerenes) return false;
    for (int i = 0; i < M; i++) {
      if (more_isomers) {
        more_isomers = BuckyGen::next_fullerene(BuckyQ, G);
        if (!more_isomers) break;
        batch::push_back(Q0s[I % Nd], G, uint64_t(I),
                         StatusFlag(int(StatusEnum::DUAL_INITIALIZED)), 0);
        I++;
      }
    }
    num_finished[GEN] = I;
    stage_finished[GEN] = num_finished[GEN] >= (int)n_fullerenes;
    return true;
  };

  auto geo_routine = [&]() {
    std::vector<int> qsize_geo;
    std::for_each(Q1s.begin(), Q1s.end(), [&](auto& Q){ qsize_geo.push_back(Q.size()); });
    for (int d = 0; d < Nd; d++) {
      if (Q0s[d].size() > 0) {
        auto& B = Bs[GEO][d];
        auto& ctx = gen_ctxs[d];
        // Pop dual-graph queue into the dual-graph batch.
        B.src_dual.clear(); B.st.clear();
        batch::queue_push(B.src_dual, B.st, Q0s[d],
                          batch::StatusPredicate{StatusEnum::DUAL_INITIALIZED});
        auto src_view = B.src_dual.view();
        auto dst_view = B.cubic_adj_view_capacity().slice(0, src_view.size());
        // Run the geometry pipeline.
        dualize<real_t, device_node_t>(ctx, src_view, dst_view, B.st.view(),
                        B.faces_cubic_span(), B.faces_dual_span()).wait();
        // After dualize, dst_cubic is filled for the active range; mark size.
        B.dst_cubic.resize(src_view.size());
        tutte_layout<real_t, device_node_t>(ctx, B.cubic_adj_view(), B.layout_span(), B.st.view()).wait();
        spherical_projection<real_t, device_node_t>(ctx, B.cubic_adj_view(),
                                     B.layout_span(), B.xyz_span(),
                                     B.st.view()).wait();
        // Push the cubic graph + xyz into Q1s for the optimization stage.
        batch::queue_push(Q1s[d], B.dst_cubic, B.st,
                          batch::StatusPredicate{StatusEnum::NOT_CONVERGED},
                          StatusEnum::EMPTY);
      }
    }
    for (int d = 0; d < Nd; d++) {
      gen_ctxs[d].wait();
      num_finished_this_round[GEO][d] = Q1s[d].size() - qsize_geo[d];
      num_finished[GEO] += num_finished_this_round[GEO][d];
    }
    stage_finished[GEO] = num_finished[GEO] >= (int)n_fullerenes;
  };

  auto opt_routine = [&](bool /*final_round_flag*/ = false) {
    std::vector<int> qsize_opt;
    std::for_each(Q2s.begin(), Q2s.end(), [&](auto& Q){ qsize_opt.push_back(Q.size()); });
    for (int d = 0; d < Nd; d++) {
      auto& B = Bs[OPT][d];
      auto& ctx = opt_ctxs[d];
      B.dst_cubic.clear(); B.st.clear();
      // Drain everything not yet 3D-converged or 3D-failed (mirror the
      // original ConditionFunctor(0,0, CONVERGED_3D|FAILED_3D|EMPTY)).
      batch::queue_push(B.dst_cubic, B.st, Q1s[d],
                        batch::StatusPredicate{StatusEnum(0),
                            StatusEnum(int(StatusEnum::CONVERGED_3D) |
                                       int(StatusEnum::FAILED_3D)    |
                                       int(StatusEnum::EMPTY))});
      forcefield_optimize<PEDERSEN, real_t, device_node_t>(ctx, B.cubic_adj_view(), B.xyz_span(),
                                  B.st.view(), 3 * N, 20 * N).wait();
      batch::queue_push(Q2s[d], B.dst_cubic, B.st,
                        batch::StatusPredicate{StatusEnum(0),
                            StatusEnum(int(StatusEnum::CONVERGED_3D) |
                                       int(StatusEnum::FAILED_3D))},
                        StatusEnum::EMPTY);
    }
    for (int d = 0; d < Nd; d++) {
      opt_ctxs[d].wait();
      num_finished_this_round[OPT][d] = Q2s[d].size() - qsize_opt[d];
      num_finished[OPT] += num_finished_this_round[OPT][d];
    }
    stage_finished[OPT] = num_finished[OPT] >= (int)n_fullerenes;
    return true;
  };

  auto prop_routine = [&]() {
    std::vector<int> qsize;
    std::for_each(Q2s.begin(), Q2s.end(), [&](auto& Q){ qsize.push_back(Q.size()); });
    for (int d = 0; d < Nd; d++) {
      auto& B = Bs[PROP][d];
      auto& ctx = opt_ctxs[d];
      B.dst_cubic.clear(); B.st.clear();
      batch::queue_push(B.dst_cubic, B.st, Q2s[d], batch::MatchAll{});
      auto bv  = B.cubic_adj_view();
      transform_to_principal_axes<real_t, device_node_t>(ctx, B.xyz_span(), (int)N, bv.size(), B.st.view()).wait();
      compute_hessian<PEDERSEN, real_t, device_node_t>(ctx, bv, B.xyz_span(), B.st.view(),
                              hessians[d], hessian_cols[d]).wait();
      eigensolve<EigensolveMode::ENDS, real_t, device_node_t>(ctx, B.xyz_span(), (int)N, bv.size(),
                            hessians[d], hessian_cols[d], n_lanczos,
                            intermediate_spectrum_ends[d],
                            std::span<real_t>(),
                            ev_off_diag[d], ev_qmat[d], ev_lanczos[d],
                            ev_diag[d], ev_ends_idx[d]).wait();
      // Inertia in Å² → multiply by carbon_mass*aangstrom_length² for SI.
      inertia<real_t, device_node_t>(ctx, B.xyz_span(), (int)N, bv.size(), B.st.view(),
                      std::span<std::array<real_t,3>>(
                          reinterpret_cast<std::array<real_t,3>*>(results[INERTIA][d].data()),
                          results[INERTIA][d].size() / 3)).wait();
      eccentricity<real_t, device_node_t>(ctx, B.xyz_span(), (int)N, bv.size(),
                           B.st.view(), results[ECCENTRICITY][d]).wait();
      volume<real_t, device_node_t>(ctx, B.xyz_span(), (int)N, bv.size(),
                     B.faces_cubic_span(),
                     std::span<uint8_t>(
                         reinterpret_cast<uint8_t*>(B.faces_dual.data()),
                         B.faces_dual.size() * sizeof(std::array<device_node_t,3>)),
                     B.st.view(), results[VOLUME][d]).wait();
    }
    for (int d = 0; d < Nd; d++) {
      opt_ctxs[d].wait();
      num_finished_this_round[PROP][d] = qsize[d] - Q2s[d].size();
      num_finished[PROP] += num_finished_this_round[PROP][d];
    }
    stage_finished[PROP] = num_finished[PROP] >= (int)n_fullerenes;
  };

  auto update_mean_var = [&](real_t v, long double& Ev_, long double& Ev2_, double K_) {
    Ev_  += v - K_;
    Ev2_ += (v - K_) * (v - K_);
  };

  auto terrible_outlier = [&](real_t v, int r) {
    return v < result_bounds[r][0] || v > result_bounds[r][1];
  };

  auto stat_routine = [&]() {
    vector<int> opt_failed, opt_not_converged;
    int i = 0;
    for (int d = 0; d < Nd; d++) {
      auto& B = Bs[PROP][d];
      auto sv = B.st.view();
      size_t n = num_finished_this_round[PROP][d];

      // First pass: update running mean/variance for converged isomers.
      for (size_t di = 0; di < n; di++) {
        auto& status = sv.status[di];
        if (status.is_set(StatusEnum::CONVERGED_3D)) {
          n_converged++;
          for (int r = 0; r < NUM_RESULTS; r++) {
            auto v = results[r][d][di];
            if (isfinite(v)) {
              if (n_converged == 1) K_acc[r] = v;
              update_mean_var(v, Ev[r], Ev2[r], K_acc[r]);
            }
          }
        }
      }
      for (int r = 0; r < NUM_RESULTS; r++) {
        current_mean[r]   = K_acc[r] + Ev[r] / std::max<size_t>(1, n_converged);
        current_stddev[r] = sqrt((Ev2[r] - Ev[r] * Ev[r]
                                / std::max<size_t>(1, n_converged))
                                / std::max<size_t>(1, n_converged - 1));
      }

      // Second pass: process per-isomer results, build candidates, histograms.
      for (size_t di = 0; di < num_finished_this_round[PROP][d]; di++, i++) {
        int   id     = (int)sv.id[di];
        auto& status = sv.status[di];

        if (int(status) != int(StatusEnum::EMPTY)) {
          num_finished_this_round[STAT][d]++;
          num_finished[STAT]++;
        }
        if (status.is_set(StatusEnum::CONVERGED_3D)) {
          array<real_t, NUM_RESULTS> R;

          real_t min_freq = sqrt(intermediate_spectrum_ends[d][di * 2    ] / carbon_mass) / (2 * M_PI) * 1e-12,
                 max_freq = sqrt(intermediate_spectrum_ends[d][di * 2 + 1] / carbon_mass) / (2 * M_PI) * 1e-12;
          min_freq *= 33.356; max_freq *= 33.356;
          results[MIN_FREQ][d][di]   = min_freq;
          results[MAX_FREQ][d][di]   = max_freq;
          results[FREQ_WIDTH][d][di] = max_freq - min_freq;

          for (int r = 0; r < NUM_RESULTS; r++)
            if (result_sizes[r] == 1)
              R[r] = results[r][d][di];

          for (int r = 0; r < NUM_RESULTS; r++)
            if (result_sizes[r] == 1) {
              auto v = results[r][d][di];
              if (terrible_outlier(v, r)) {
                v = results[r][d][di] = nan("");
                terrible_outliers[r].push_back(id);
                continue;
              }
              if (isfinite(v)) {
                auto cubic_spans = B.dst_cubic.view().spans();
                auto& adj = std::get<0>(cubic_spans);
                auto& pts = std::get<3>(cubic_spans);
                isomer_candidate C(0, id, R,
                                   adj.data() + di * size_t(N) * 3,
                                   pts.data() + di * size_t(N),
                                   N);
                C.value =  v; result_min[r].insert(C);
                C.value = -v; result_max[r].insert(C);
                means[r] += v;
              }

              auto vmin = result_bounds[r][0], vmax = result_bounds[r][1];
              size_t v_bin = std::min(num_bins - 1, std::max(0LU, (size_t)floor(num_bins * ((v - vmin) / (vmax - vmin)))));
              result_histograms[r][v_bin]++;
              for (int s = 0; s < r; s++)
                if (result_sizes[s] == 1) {
                  auto w = results[s][d][di];
                  if (isfinite(w) && !terrible_outlier(w, s)) {
                    auto wmin = result_bounds[s][0], wmax = result_bounds[s][1];
                    size_t w_bin = std::min(num_bins - 1, std::max((size_t)0, (size_t)floor(num_bins * ((w - wmin) / (wmax - wmin)))));
                    result_histograms2D[(r * (r - 1)) / 2 + s][v_bin * num_bins + w_bin]++;
                  }
                }
            }
        } else {
          if (status.is_set(StatusEnum::FAILED_3D))    opt_failed.push_back(id);
          if (status.is_set(StatusEnum::NOT_CONVERGED)) opt_not_converged.push_back(id);
        }
      }
    }

    cout << "num_finished_this_round = " << num_finished_this_round << "\n";
    cout << "num_finished            = " << num_finished << "\n";
    cout << "opt_failed        = " << opt_failed << "\n"
         << "opt_not_converged = " << opt_not_converged << "\n";
    if (num_finished[STAT] >= (int)n_fullerenes)
      stage_finished[STAT] = true;
  };

  auto final_round = [&](int s) {
    bool final = !stage_finished[s];
    for (int t = GEN; t < s; t++) final &= stage_finished[t];
    return final;
  };

  while (!stage_finished[PROP]) {
    cout << loop_iters << " start: Gen: " << num_finished[GEN]
         << "  Geometry: " << num_finished[GEO]
         << "  Opt: "      << num_finished[OPT]
         << "  Prop: "     << num_finished[PROP] << std::endl;

    auto T0L = steady_clock::now();
    geo_routine();
    auto T1 = steady_clock::now(); Tinit_geom += T1 - T0L;
    auto handle = std::async(std::launch::async, generate_isomers, 2 * (int)batch_size);
    auto T2 = steady_clock::now();

    while (Q1s[0].size() > Bs[OPT][0].capacity() && Q1s[1 % Nd].size() > Bs[OPT][1 % Nd].capacity())
      opt_routine();

    if (Q2s[0].size() > Bs[PROP][0].capacity() && Q2s[1 % Nd].size() > Bs[PROP][1 % Nd].capacity())
      prop_routine();

    auto T3 = steady_clock::now(); Topt += T3 - T2;
    stat_routine();
    handle.wait();
    auto T4 = steady_clock::now(); Tgen += T4 - T3;

    ssize_t final_step_max = 100, num_orphans = 0;
    if (final_round(OPT)) {
      for (ssize_t fs = 0; fs < final_step_max && final_round(OPT); fs++)
        opt_routine(true);
      if (!stage_finished[OPT]) {
        stage_finished[OPT] = true;
        num_orphans = n_fullerenes - num_finished[OPT];
        fprintf(stderr, "Finishing final round of optimizing, %ld isomers orphaned.\n", num_orphans);
      }
    }
    if (final_round(PROP)) prop_routine();
    if (final_round(STAT)) stat_routine();

    if (loop_iters % 3 == 0 || stage_finished[STAT]) {
      Tstart = steady_clock::now();
      n0 = num_finished[STAT];
    }
    loop_iters++;

    cout << loop_iters << " end:   Gen: " << num_finished[GEN]
         << "  Geometry: " << num_finished[GEO]
         << "  Opt: "      << num_finished[OPT]
         << "  Prop: "     << num_finished[PROP] << std::endl;

    auto T5 = steady_clock::now(); Topt += T5 - T4;
  }
  cout << "Exited loop, waiting on op_ctxs.\n";
  for (int d = 0; d < Nd; d++) opt_ctxs[d].wait();

  array<vector<isomer_candidate>, NUM_RESULTS> smallest_candidates, biggest_candidates;
  for (int r = 0; r < NUM_RESULTS; r++) {
    means[r] /= std::max(1, num_finished[PROP]);
    smallest_candidates.at(r) = sorted(result_min.at(r).as_vector());
    biggest_candidates.at(r)  = sorted(result_max.at(r).as_vector());
  }

  cout << "COMPLETED IN " << loop_iters << " rounds.\n";
  cout << "num_finished            = " << num_finished << "\n";
  for (int r = 0; r < NUM_RESULTS; r++)
    if (result_sizes[r] == 1) {
      auto& smallest = smallest_candidates[r], &biggest = biggest_candidates[r];
      cout << result_names[r] << "_min = " << smallest << "\n"
           << "   (removed " << terrible_outliers[r].size() << " outliers)\n"
           << result_names[r] << "_max = " << biggest << "\n\n";
    }
  cout.flush();

  // GATHER SELECTED CANDIDATES into Q3s and recompute Hessian/eigensystem.
  auto host_graph = [](const vector<device_node_t>& device_graph) {
      // Wrap the flat stride-3 neighbour array as an RSRAdjacencyView and
      // convert to an owning Graph.
      const int Nh = int(device_graph.size() / 3);
      std::vector<uint8_t> deg_buf(Nh, 3);
      Spanify::RSRAdjacencyView<device_node_t> v(
          device_node_t(Nh), /*dmax=*/3,
          std::span<device_node_t>(const_cast<device_node_t*>(device_graph.data()), device_graph.size()),
          std::span<uint8_t>(deg_buf.data(), deg_buf.size()));
      return batch::to_graph(v);
  };
  auto host_points = [](const vector<real_t>& X) {
      const size_t Nh = X.size() / 3;
      return batch::to_points(std::span<const std::array<real_t,3>>(
          reinterpret_cast<const std::array<real_t,3>*>(X.data()), Nh));
  };

  // Helper: push a candidate (cubic adjacency + xyz held by isomer_candidate)
  // into a polyhedron queue (graph + xyz in one slot).
  auto push_candidate = [&](batch::BatchQueue<PV>& Q, const isomer_candidate& C,
                            uint64_t id, StatusFlag flag) {
      // Wrap the candidate's cubic adjacency as RSRAdjacencyView<device_node_t>.
      std::vector<uint8_t> deg_buf(N, 3);
      Spanify::RSRAdjacencyView<device_node_t> src_view(
          device_node_t(N), /*dmax=*/3,
          std::span<device_node_t>(const_cast<device_node_t*>(C.cubic_neighbours.data()),
                                   C.cubic_neighbours.size()),
          std::span<uint8_t>(deg_buf.data(), deg_buf.size()));
      auto xyz_span = std::span<const std::array<real_t,3>>(
          reinterpret_cast<const std::array<real_t,3>*>(C.X.data()), N);
      batch::push_back(Q, src_view, xyz_span, id, flag, 0);
  };

  // Push candidates' graphs+xyz into Q3s.
  for (int r = 0; r < NUM_RESULTS; r++)
    if (result_sizes[r] == 1) {
      auto& smallest = smallest_candidates.at(r), &biggest = biggest_candidates.at(r);
      for (int i_ = 0; i_ < (int)n_best_candidates; i_++)
        push_candidate(Q3s[0], smallest.at(i_),
                       uint64_t(smallest.at(i_).id),
                       StatusFlag(int(StatusEnum::CONVERGED_3D)));
      for (int i_ = 0; i_ < (int)n_best_candidates; i_++)
        push_candidate(Q3s[1 % Nd], biggest.at(i_),
                       uint64_t(biggest.at(i_).id),
                       StatusFlag(int(StatusEnum::CONVERGED_3D)));
    }

  if (STORE_HESSIAN) {
    for (int d = 0; d < Nd; d++) {
      auto& B   = final_batch[d];
      auto& ctx = opt_ctxs[d];
      ctx.wait();
      B.dst_cubic.clear(); B.st.clear();
      batch::queue_push(B.dst_cubic, B.st, Q3s[d],
                        batch::StatusPredicate{StatusEnum::CONVERGED_3D});
      auto bv = B.cubic_adj_view();
      compute_hessian<PEDERSEN, real_t, device_node_t>(ctx, bv, B.xyz_span(), B.st.view(),
                              hessians[d], hessian_cols[d]).wait();
      if (STORE_EIGENSYSTEM) {
        Q_eig[d] = SyclVector<real_t>(3 * N * 3 * N * final_batch_size, 0);
        lams[d]  = SyclVector<real_t>(3 * N * final_batch_size, 0);
        SyclVector<real_t>          fs_off_diag(n_lanczos * final_batch_size, 0);
        SyclVector<real_t>          fs_qmat   (n_lanczos * n_lanczos * final_batch_size, 0);
        SyclVector<real_t>          fs_lanczos(n_lanczos * N * 3 * final_batch_size, 0);
        SyclVector<real_t>          fs_diag   (n_lanczos * final_batch_size, 0);
        SyclVector<device_node_t>   fs_ends   (2 * final_batch_size, 0);
        eigensolve<EigensolveMode::FULL_SPECTRUM_VECTORS, real_t, device_node_t>(ctx, B.xyz_span(), (int)N, bv.size(),
                      hessians[d], hessian_cols[d], N*3 - 6,
                      std::span<real_t>(lams[d].data(), lams[d].size()),
                      std::span<real_t>(Q_eig[d].data(), Q_eig[d].size()),
                      std::span<real_t>(fs_off_diag.data(), fs_off_diag.size()),
                      std::span<real_t>(fs_qmat.data(), fs_qmat.size()),
                      std::span<real_t>(fs_lanczos.data(), fs_lanczos.size()),
                      std::span<real_t>(fs_diag.data(), fs_diag.size()),
                      std::span<device_node_t>(fs_ends.data(), fs_ends.size())).wait();
      }
      ctx.wait();
    }
  }

  if (mkdir(output_dir.c_str(), 0777) != 0 && errno != EEXIST) {
    perror("Error creating output directory");
    return -1;
  }
  {
    std::ofstream f(output_dir + "/result_bounds.float64", std::ios::binary);
    if (!f) { std::cerr << "Error opening file for writing: " << output_dir + "/result_bounds.float64" << std::endl; return -1; }
    f.write(reinterpret_cast<const char*>(&result_bounds), sizeof(result_bounds));
  }

  FILE* f = fopen((output_dir + "/result_names.string").c_str(), "w");
  for (int r = 0; r < NUM_RESULTS; r++) fprintf(f, "%s\n", result_names[r].c_str());
  fclose(f);

  FILE* f_min   = fopen((output_dir + "/result_mins.float32").c_str(), "wb");
  FILE* f_max   = fopen((output_dir + "/result_maxs.float32").c_str(), "wb");
  FILE* f_hist  = fopen((output_dir + "/result_hists.uint64").c_str(), "wb");
  FILE* f_hist2D= fopen((output_dir + "/result_hists2D.uint64").c_str(), "wb");

  for (int r = 0; r < NUM_RESULTS; r++)
    if (result_sizes[r] == 1) {
      string result_dir = output_dir + "/" + result_names[r];
      string min_dir = result_dir + "/min", max_dir = result_dir + "/max",
             outlier_dir = result_dir + "/terrible_outliers";
      auto& smallest = smallest_candidates[r], &biggest = biggest_candidates[r];

      auto result_index = [&](bool big, size_t r_, size_t i_) {
        return array<size_t, 2>{{big ? 1u : 0u, r_ * n_best_candidates + i_}};
      };

      vector<real_t> smallest_values(n_best_candidates), biggest_values(n_best_candidates);
      for (size_t i_ = 0; i_ < n_best_candidates; i_++) {
        smallest_values[i_] = smallest[i_].value;
        biggest_values[i_]  = biggest[i_].value;
      }
      fwrite(smallest_values.data(), sizeof(real_t), n_best_candidates, f_min);
      fwrite(biggest_values.data(),  sizeof(real_t), n_best_candidates, f_max);

      mkdir(result_dir.c_str(), 0777);
      mkdir(min_dir.c_str(), 0777);
      mkdir(max_dir.c_str(), 0777);
      mkdir(outlier_dir.c_str(), 0777);

      for (size_t i_ = 0; i_ < n_best_candidates; i_++) {
        assert(i_ < smallest.size()); assert(i_ < biggest.size());
        string min_basename = min_dir + "/Pmin" + to_string(i_),
               max_basename = max_dir + "/Pmax" + to_string(i_);
        isomer_candidate Cmin = smallest[i_], Cmax = biggest[i_];
        auto [d_min, di_min] = result_index(0, r, i_);
        auto [d_max, di_max] = result_index(1, r, i_);

        // Build host Polyhedrons from the candidates' stored data; the
        // final_batch layout depends on insertion ordering and cannot be
        // directly slot-indexed in the new design without a side index map.
        Polyhedron Pmin(PlanarGraph(host_graph(Cmin.cubic_neighbours)), host_points(Cmin.X));
        Polyhedron Pmax(PlanarGraph(host_graph(Cmax.cubic_neighbours)), host_points(Cmax.X));

        Polyhedron::to_file(Pmin, min_basename + ".mol2");
        Polyhedron::to_file(Pmin, min_basename + ".spiral");
        Polyhedron::to_file(Pmax, max_basename + ".mol2");
        Polyhedron::to_file(Pmax, max_basename + ".spiral");

        {
          FILE* fp = fopen((min_basename + "-graph.uint16").c_str(), "wb");
          fwrite(Cmin.cubic_neighbours.data(), sizeof(device_node_t), 3 * N, fp);
          fclose(fp);
          fp = fopen((max_basename + "-graph.uint16").c_str(), "wb");
          fwrite(Cmax.cubic_neighbours.data(), sizeof(device_node_t), 3 * N, fp);
          fclose(fp);
          fp = fopen((min_basename + "-PX.float32").c_str(), "wb");
          fwrite(Cmin.X.data(), sizeof(real_t), 3 * N, fp);
          fclose(fp);
          fp = fopen((max_basename + "-PX.float32").c_str(), "wb");
          fwrite(Cmax.X.data(), sizeof(real_t), 3 * N, fp);
          fclose(fp);

          // -X.float32: same as -PX (post-Hessian xyz lives on device in
          // final_batch[*] but not directly slot-indexable here without a map).
          fp = fopen((min_basename + "-X.float32").c_str(), "wb");
          fwrite(Cmin.X.data(), sizeof(real_t), 3 * N, fp);
          fclose(fp);
          fp = fopen((max_basename + "-X.float32").c_str(), "wb");
          fwrite(Cmax.X.data(), sizeof(real_t), 3 * N, fp);
          fclose(fp);

          if (STORE_HESSIAN) {
            // Use d_min/d_max device's hessian buffer at slot di_max as a
            // reasonable approximation to the original code (which also used
            // di_max for both min and max writes).
            fp = fopen((min_basename + "-hessians.float32").c_str(), "wb");
            fwrite(&hessians[d_min].data()[di_max * 3 * 3 * 10 * N], sizeof(real_t), 3 * 3 * 10 * N, fp);
            fclose(fp);
            fp = fopen((max_basename + "-hessians.float32").c_str(), "wb");
            fwrite(&hessians[d_max].data()[di_max * 3 * 3 * 10 * N], sizeof(real_t), 3 * 3 * 10 * N, fp);
            fclose(fp);
            fp = fopen((min_basename + "-hessian_cols.uint16").c_str(), "wb");
            fwrite(&hessian_cols[d_min].data()[di_max * 3 * 3 * 10 * N], sizeof(uint16_t), 3 * 3 * 10 * N, fp);
            fclose(fp);
            fp = fopen((max_basename + "-hessian_cols.uint16").c_str(), "wb");
            fwrite(&hessian_cols[d_max].data()[di_max * 3 * 3 * 10 * N], sizeof(uint16_t), 3 * 3 * 10 * N, fp);
            fclose(fp);

            if (STORE_EIGENSYSTEM) {
              fp = fopen((min_basename + "-hessian_Q.float32").c_str(), "wb");
              fwrite(&Q_eig[d_min].data()[di_max * 3 * N * 3 * N], sizeof(real_t), 3 * N * 3 * N, fp);
              fclose(fp);
              fp = fopen((min_basename + "-hessian_lams.float32").c_str(), "wb");
              fwrite(&lams[d_min].data()[di_max * 3 * N], sizeof(real_t), 3 * N, fp);
              fclose(fp);
              fp = fopen((max_basename + "-hessian_Q.float32").c_str(), "wb");
              fwrite(&Q_eig[d_max].data()[di_max * 3 * N * 3 * N], sizeof(real_t), 3 * N * 3 * N, fp);
              fclose(fp);
              fp = fopen((max_basename + "-hessian_lams.float32").c_str(), "wb");
              fwrite(&lams[d_max].data()[di_max * 3 * N], sizeof(real_t), 3 * N, fp);
              fclose(fp);
            }
          }

          fp = fopen((min_basename + "-results.float32").c_str(), "wb");
          fwrite(Cmin.results.data(), sizeof(real_t), NUM_RESULTS - 1, fp);
          fclose(fp);
          fp = fopen((max_basename + "-results.float32").c_str(), "wb");
          fwrite(Cmax.results.data(), sizeof(real_t), NUM_RESULTS - 1, fp);
          fclose(fp);
        }
      }

      fwrite(&result_histograms[r][0], sizeof(uint64_t), num_bins, f_hist);
      for (int s = 0; s < r; s++)
        fwrite(&result_histograms2D[(r * (r - 1)) / 2 + s][0],
               sizeof(uint64_t), num_bins * num_bins, f_hist2D);
    }

  fclose(f_min); fclose(f_max); fclose(f_hist); fclose(f_hist2D);

  auto Ttot = steady_clock::now() - T0;
  (void)Ttot; (void)Tupdate; (void)Tqueue; (void)Tdual; (void)Ttutte;
  (void)TX0; (void)Tcopy; (void)Tfile; (void)Toutq; (void)Tinq; (void)stdevs;
  (void)current_mean; (void)current_stddev; (void)aangstrom_length;
  (void)surface_area<real_t, device_node_t>; (void)result_reference; (void)progress_bar;
  (void)stage_names;

  cout << "\n\n\n\n\n\n\n\n";
  BuckyGen::stop(BuckyQ);
  return 0;
}
