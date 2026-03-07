// Benchmark: ExtPathOptimized quality across M evenly-spaced isomers of C_N.
// Uses buckygen enumeration with strided sampling (every total/M-th isomer).
// Writes results incrementally: JSONL after each isomer, full JSON periodically.
// Supports resumption from partial output.
//
// Build: cmake --build build --target bench_epopt
// Run:   ./build/benchmarks/bench_epopt 60 1812 benchmarks/bench_C60.json

#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/deltahedron.hh"
#include "fullerenes/buckinverse.hh"
#include "fullerenes/isomerdb.hh"
#include "fullerenes/stats.hh"
#include <chrono>
#include <cmath>
#include <fstream>
#include <cstdio>
#include <algorithm>
#include <string>

using namespace std;
using namespace std::chrono;
using namespace buckinverse;

struct IsomerStats {
    int idx;
    double ms;
    int iters;
    bool converged;
    double edge_cv;
    double edge_relerr_max;  // max |L_i - L_mean| / L_mean
    double h_min;
    int n_concave;
    double ang_min, ang_max, ang_std;
    double ang_relerr_max;   // max |angle - 60| / 60
    double gmax_L;
    int seed;
    int n_steps;
};

static IsomerStats compute_stats(const Deltahedron& D, int idx, int seed, int n_steps) {
    IsomerStats s{};
    s.idx = idx;
    s.seed = seed;
    s.n_steps = n_steps;

    // Edge lengths
    vector<double> edge_lens;
    edge_lens.reserve(D.N * 3);
    for (int u = 0; u < D.N; u++)
        for (int v : D.neighbours[u])
            if (v > u) edge_lens.push_back((D.points[u] - D.points[v]).norm());
    s.edge_cv = cv_twopass(edge_lens);
    double L_mean = accumulate(edge_lens.begin(), edge_lens.end(), 0.0) / edge_lens.size();
    s.edge_relerr_max = 0;
    for (double l : edge_lens)
        s.edge_relerr_max = max(s.edge_relerr_max, fabs(l - L_mean) / L_mean);

    // Convexity
    s.h_min = INFINITY;
    s.n_concave = 0;
    for (int v = 0; v < D.N; v++) {
        int d = (int)D.neighbours[v].size();
        if (d > 6) continue;
        coord3d centroid(0,0,0);
        for (int i = 0; i < d; i++) centroid += D.points[D.neighbours[v][i]];
        centroid /= (double)d;
        coord3d n_fan(0,0,0);
        for (int i = 0; i < d; i++) {
            coord3d e1 = D.points[D.neighbours[v][i]] - D.points[v];
            coord3d e2 = D.points[D.neighbours[v][(i+1)%d]] - D.points[v];
            n_fan += e1.cross(e2);
        }
        double n_len = n_fan.norm();
        if (n_len < 1e-15) continue;
        double h = (D.points[v] - centroid).dot(n_fan / n_len);
        if (h < s.h_min) s.h_min = h;
        if (h < -1e-6) s.n_concave++;
    }

    // Triangle angles
    vector<double> angles;
    s.ang_min = 180; s.ang_max = 0;
    for (const auto& tri : D.triangles)
        for (int c = 0; c < 3; c++) {
            coord3d va = D.points[tri[(c+1)%3]] - D.points[tri[c]];
            coord3d vb = D.points[tri[(c+2)%3]] - D.points[tri[c]];
            double ang = coord3d::angle(va, vb) * 180.0 / M_PI;
            angles.push_back(ang);
            s.ang_min = min(s.ang_min, ang); s.ang_max = max(s.ang_max, ang);
        }
    s.ang_std = stddev_twopass(angles);
    s.ang_relerr_max = max(60.0 - s.ang_min, s.ang_max - 60.0) / 60.0;

    return s;
}

// Write one JSONL line for a single isomer result
static void append_jsonl(FILE* f, const IsomerStats& s) {
    fprintf(f, "{\"idx\":%d,\"ms\":%.2f,\"iters\":%d,\"converged\":%s,"
            "\"edge_cv\":%.6e,\"edge_relerr_max\":%.6e,"
            "\"h_min\":%.6f,\"n_concave\":%d,"
            "\"ang_min\":%.3f,\"ang_max\":%.3f,\"ang_std\":%.6e,\"ang_relerr_max\":%.6e,"
            "\"gmax_L\":%.6e,\"seed\":%d,\"n_steps\":%d}\n",
            s.idx, s.ms, s.iters, s.converged ? "true" : "false",
            s.edge_cv, s.edge_relerr_max,
            s.h_min, s.n_concave,
            s.ang_min, s.ang_max, s.ang_std, s.ang_relerr_max,
            s.gmax_L, s.seed, s.n_steps);
    fflush(f);
}

// Read existing JSONL file to recover partial results for resumption
static vector<IsomerStats> read_jsonl(const string& path) {
    vector<IsomerStats> results;
    FILE* f = fopen(path.c_str(), "r");
    if (!f) return results;

    char line[4096];
    while (fgets(line, sizeof(line), f)) {
        IsomerStats s{};
        char conv_str[8] = {};
        int n = sscanf(line,
            "{\"idx\":%d,\"ms\":%lf,\"iters\":%d,\"converged\":%7[^,],"
            "\"edge_cv\":%lf,\"edge_relerr_max\":%lf,"
            "\"h_min\":%lf,\"n_concave\":%d,"
            "\"ang_min\":%lf,\"ang_max\":%lf,\"ang_std\":%lf,\"ang_relerr_max\":%lf,"
            "\"gmax_L\":%lf,\"seed\":%d,\"n_steps\":%d}",
            &s.idx, &s.ms, &s.iters, conv_str,
            &s.edge_cv, &s.edge_relerr_max,
            &s.h_min, &s.n_concave,
            &s.ang_min, &s.ang_max, &s.ang_std, &s.ang_relerr_max,
            &s.gmax_L, &s.seed, &s.n_steps);
        if (n >= 13) {
            s.converged = (conv_str[0] == 't');
            results.push_back(s);
        }
    }
    fclose(f);
    return results;
}

static void write_json(const char* path, int N, int total_isomers, int stride,
                       const vector<IsomerStats>& stats, double elapsed_ms) {
    FILE* f = fopen(path, "w");
    if (!f) { fprintf(stderr, "Cannot open %s\n", path); return; }

    int dual_N = N/2 + 2;

    // Summary statistics
    double iters_sum = 0;
    double ecv_max = 0, ermax_max = 0;
    double hmin_min = INFINITY;
    int concave_total = 0, concave_isomers = 0;
    double amin_min = 180, amax_max = 0;
    double astd_max = 0, armax_max = 0;
    double gmax_max = 0;
    int converged_count = 0;
    int n_bad_edge = 0, n_bad_ang = 0;

    for (const auto& s : stats) {
        iters_sum += s.iters;
        ecv_max = max(ecv_max, s.edge_cv);
        ermax_max = max(ermax_max, s.edge_relerr_max);
        if (s.h_min < hmin_min) hmin_min = s.h_min;
        if (s.n_concave > 0) concave_isomers++;
        concave_total += s.n_concave;
        amin_min = min(amin_min, s.ang_min);
        amax_max = max(amax_max, s.ang_max);
        astd_max = max(astd_max, s.ang_std);
        armax_max = max(armax_max, s.ang_relerr_max);
        gmax_max = max(gmax_max, s.gmax_L);
        if (s.converged) converged_count++;
        if (s.edge_relerr_max > 0.01) n_bad_edge++;
        if (s.ang_relerr_max > 0.01) n_bad_ang++;
    }

    int actual_M = (int)stats.size();
    auto ms_mean  = actual_M > 0 ? accumulate(stats.begin(), stats.end(), 0.0, [](double a, const IsomerStats& s){ return a + s.ms; }) / actual_M : 0.0;
    auto ecv_mean = actual_M > 0 ? accumulate(stats.begin(), stats.end(), 0.0, [](double a, const IsomerStats& s){ return a + s.edge_cv; }) / actual_M : 0.0;
    auto astd_mean= actual_M > 0 ? accumulate(stats.begin(), stats.end(), 0.0, [](double a, const IsomerStats& s){ return a + s.ang_std; }) / actual_M : 0.0;
    auto gmax_mean= actual_M > 0 ? accumulate(stats.begin(), stats.end(), 0.0, [](double a, const IsomerStats& s){ return a + s.gmax_L; }) / actual_M : 0.0;
    double ms_std   = stddev_twopass(actual_M, [&](int i){ return stats[i].ms; });
    double ecv_std  = stddev_twopass(actual_M, [&](int i){ return stats[i].edge_cv; });
    double astd_std = stddev_twopass(actual_M, [&](int i){ return stats[i].ang_std; });
    double gmax_std = stddev_twopass(actual_M, [&](int i){ return stats[i].gmax_L; });
    double iters_per_vertex = actual_M > 0 ? iters_sum / ((double)dual_N * actual_M) : 0;

    fprintf(f, "{\n");
    fprintf(f, "  \"N\": %d,\n", N);
    fprintf(f, "  \"M\": %d,\n", actual_M);
    fprintf(f, "  \"total_isomers\": %d,\n", total_isomers);
    fprintf(f, "  \"stride\": %d,\n", stride);
    fprintf(f, "  \"dual_N\": %d,\n", dual_N);
    fprintf(f, "  \"total_ms\": %.1f,\n", elapsed_ms);
    fprintf(f, "  \"summary\": {\n");
    fprintf(f, "    \"ms_mean\": %.3f, \"ms_std\": %.3f,\n", ms_mean, ms_std);
    fprintf(f, "    \"iters_per_vertex\": %.1f,\n", iters_per_vertex);
    fprintf(f, "    \"converged_frac\": %.4f,\n", actual_M > 0 ? (double)converged_count / actual_M : 0.0);
    fprintf(f, "    \"edge_cv_mean\": %.6e, \"edge_cv_std\": %.6e, \"edge_cv_max\": %.6e,\n",
            ecv_mean, ecv_std, ecv_max);
    fprintf(f, "    \"edge_relerr_max\": %.6e,\n", ermax_max);
    fprintf(f, "    \"h_min_worst\": %.6f,\n", hmin_min);
    fprintf(f, "    \"concave_isomers\": %d, \"concave_vertices_total\": %d,\n",
            concave_isomers, concave_total);
    fprintf(f, "    \"ang_min_worst\": %.2f, \"ang_max_worst\": %.2f,\n", amin_min, amax_max);
    fprintf(f, "    \"ang_std_mean\": %.6e, \"ang_std_std\": %.6e, \"ang_std_max\": %.6e,\n",
            astd_mean, astd_std, astd_max);
    fprintf(f, "    \"ang_relerr_max\": %.6e,\n", armax_max);
    fprintf(f, "    \"gmax_L_mean\": %.6e, \"gmax_L_std\": %.6e, \"gmax_L_max\": %.6e,\n",
            gmax_mean, gmax_std, gmax_max);
    fprintf(f, "    \"n_bad_edge\": %d, \"n_bad_ang\": %d, \"n_not_converged\": %d\n",
            n_bad_edge, n_bad_ang, actual_M - converged_count);
    fprintf(f, "  },\n");

    // Per-isomer data
    fprintf(f, "  \"per_isomer\": [\n");
    for (int i = 0; i < actual_M; i++) {
        const auto& s = stats[i];
        fprintf(f, "    {\"idx\":%d,\"ms\":%.2f,\"iters\":%d,\"converged\":%s,"
                "\"edge_cv\":%.6e,\"edge_relerr_max\":%.6e,"
                "\"h_min\":%.6f,\"n_concave\":%d,"
                "\"ang_min\":%.3f,\"ang_max\":%.3f,\"ang_std\":%.6e,\"ang_relerr_max\":%.6e,"
                "\"gmax_L\":%.6e,\"seed\":%d,\"n_steps\":%d}%s\n",
                s.idx, s.ms, s.iters, s.converged ? "true" : "false",
                s.edge_cv, s.edge_relerr_max,
                s.h_min, s.n_concave,
                s.ang_min, s.ang_max, s.ang_std, s.ang_relerr_max,
                s.gmax_L, s.seed, s.n_steps,
                (i + 1 < actual_M) ? "," : "");
    }
    fprintf(f, "  ]\n");
    fprintf(f, "}\n");
    fclose(f);
}

int main(int argc, char** argv) {
    if (argc < 4) {
        fprintf(stderr, "Usage: %s N M output.json [total_isomers]\n", argv[0]);
        fprintf(stderr, "\n");
        fprintf(stderr, "  N              - fullerene size (number of carbon atoms)\n");
        fprintf(stderr, "  M              - number of isomers to sample\n");
        fprintf(stderr, "  output.json    - output file (also writes output.jsonl incrementally)\n");
        fprintf(stderr, "  total_isomers  - total isomer count (optional, looked up from IsomerDB)\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "Samples every total/M-th isomer from the buckygen enumeration.\n");
        fprintf(stderr, "Results are written incrementally to output.jsonl; if interrupted,\n");
        fprintf(stderr, "re-running the same command resumes from where it left off.\n");
        return 1;
    }
    int N = atoi(argv[1]);
    int M = atoi(argv[2]);
    const char* json_path = argv[3];

    // Determine total isomer count
    int total_isomers = 0;
    if (argc >= 5) {
        total_isomers = atoi(argv[4]);
    } else {
        int64_t db_count = IsomerDB::number_isomers(N, "Any", false);
        if (db_count > 0) {
            total_isomers = (int)db_count;
        } else {
            fprintf(stderr, "C%d not in IsomerDB. Provide total_isomers as 4th argument.\n", N);
            return 1;
        }
    }

    if (M > total_isomers) M = total_isomers;
    int stride = max(1, total_isomers / M);
    // Actual number of samples: ceil(total_isomers / stride)
    int planned_M = (total_isomers + stride - 1) / stride;

    fprintf(stderr, "C%d: %d total isomers, sampling %d (stride %d)\n",
            N, total_isomers, planned_M, stride);

    // JSONL path for incremental output
    string jsonl_path = string(json_path);
    // Replace .json with .jsonl, or append .jsonl
    if (jsonl_path.size() >= 5 && jsonl_path.substr(jsonl_path.size()-5) == ".json")
        jsonl_path = jsonl_path.substr(0, jsonl_path.size()-5) + ".jsonl";
    else
        jsonl_path += "l";

    // Check for existing partial results (resumption)
    vector<IsomerStats> stats = read_jsonl(jsonl_path);
    int n_done = (int)stats.size();
    if (n_done > 0) {
        fprintf(stderr, "Resuming: %d samples already completed in %s\n",
                n_done, jsonl_path.c_str());
        if (n_done >= planned_M) {
            fprintf(stderr, "All %d samples already done. Writing final JSON.\n", planned_M);
            double elapsed_ms = 0;
            for (const auto& s : stats) elapsed_ms += s.ms;
            write_json(json_path, N, total_isomers, stride, stats, elapsed_ms);
            fprintf(stderr, "Written to %s\n", json_path);
            return 0;
        }
    }

    // Open JSONL for appending
    FILE* jsonl_f = fopen(jsonl_path.c_str(), "a");
    if (!jsonl_f) {
        fprintf(stderr, "Cannot open %s for writing\n", jsonl_path.c_str());
        return 1;
    }

    // The first sample we need is at enumeration index n_done * stride.
    // We must enumerate through buckygen to reach it.
    int next_sample_enum = n_done * stride;  // 0-based enumeration index of next sample
    int sample_count = n_done;               // samples completed so far

    fprintf(stderr, "Starting buckygen enumeration (skipping to index %d)...\n", next_sample_enum);
    fflush(stderr);

    auto t_start = steady_clock::now();
    auto Q = BuckyGen::start(N, false);
    Graph G;
    int enum_idx = 0;

    while (BuckyGen::next_fullerene(Q, G)) {
        if (enum_idx == next_sample_enum) {
            // Process this isomer
            ReducibleDual rd(G);
            ExtensionPath ep = rd.reduceToExtensionPath(20);
            if (ep.seed == SeedType::NotASeed) {
                fprintf(stderr, "  isomer %d: failed to reduce, skipping\n", enum_idx);
            } else {
                auto t0 = steady_clock::now();
                Deltahedron D = Deltahedron::fromExtensionPathOptimized(ep);
                auto t1 = steady_clock::now();

                IsomerStats s = compute_stats(D, enum_idx, (int)ep.seed, (int)ep.steps.size());
                s.ms = duration<double,milli>(t1 - t0).count();
                s.iters = D.iterations_used;
                s.converged = (D.final_gmax_L < 1e-5);
                s.gmax_L = D.final_gmax_L;
                stats.push_back(s);

                // Write incrementally to JSONL
                append_jsonl(jsonl_f, s);
            }

            sample_count++;

            // Progress report
            if (sample_count % 100 == 0 || sample_count == planned_M) {
                double elapsed = duration<double>(steady_clock::now() - t_start).count();
                int remaining = planned_M - sample_count;
                double ms_per = elapsed * 1000.0 / (sample_count - n_done);
                double eta_s = remaining * ms_per / 1000.0;
                fprintf(stderr, "  %d/%d samples (enum %d/%d) %.0f ms/iso ETA %.0fs\n",
                        sample_count, planned_M, enum_idx+1, total_isomers,
                        ms_per, eta_s);
            }

            // Write full JSON periodically (every 500 samples)
            if (sample_count % 500 == 0) {
                double elapsed_ms = duration<double,milli>(steady_clock::now() - t_start).count();
                // Add time from resumed samples
                for (int i = 0; i < n_done; i++) elapsed_ms += stats[i].ms;
                write_json(json_path, N, total_isomers, stride, stats, elapsed_ms);
            }

            // Advance to next sample
            next_sample_enum += stride;
            if (next_sample_enum >= total_isomers) break;
        }

        enum_idx++;
        if (enum_idx % 1000000 == 0 && enum_idx < next_sample_enum)
            fprintf(stderr, "  enumerating %dM / %dM (next sample at %d)...\n",
                    enum_idx / 1000000, total_isomers / 1000000, next_sample_enum);
    }
    BuckyGen::stop(Q);
    fclose(jsonl_f);

    // Write final JSON
    double elapsed_ms = duration<double,milli>(steady_clock::now() - t_start).count();
    for (int i = 0; i < n_done; i++) elapsed_ms += stats[i].ms;
    write_json(json_path, N, total_isomers, stride, stats, elapsed_ms);
    fprintf(stderr, "Done: %d samples in %.1f s. Written to %s\n",
            (int)stats.size(), elapsed_ms / 1000.0, json_path);
    return 0;
}
