// programs/eisenstein_paint/validate
//
// Sweep all general (non-IPR) fullerene duals C20..Cmax via buckygen,
// run the Eisenstein-paint pipeline on each, dump per-N stats + per-
// failure records (incl. raw triangulation .tri dumps for offline
// reproduction).
//
// Usage: eisenstein_paint_validate [-r start_N] [--no-alex] [max_N] [results_dir]
//
//   -r N        resume from size N (skip sizes < N, append to summary.txt)
//   --no-alex   only test sorted_dual + dual_idt (sort + iDT compute +
//               chartability check); ~20x faster, useful for filtering
//               inputs whose raw iDT is non-chartable.
//   defaults: max_N=200, results_dir=./validate_results
//
// Per-N _fails.txt records: idx, stage_name, why.
// Per-N _fails.tri dumps the failing triangulations (header + neighbour
// lists) so downstream tools can reproduce without re-enumerating.

#include "fullerenes/eisenstein_paint_geometry.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/buckygen-wrapper.hh"
#if defined(_OPENMP)
#include <omp.h>
#endif
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <numeric>
#include <string>
#include <vector>

namespace fs  = std::filesystem;
namespace ep  = eisenstein_paint;
using clk     = std::chrono::high_resolution_clock;

static constexpr int BATCH = 2048;

static void dump_triangulation(FILE* f, const Triangulation& T, long long idx,
                               int N_carbon, const char* stage)
{
    std::fprintf(f, "N %d idx %lld fail %s  dual_N %d\n",
                 N_carbon, idx, stage, T.N);
    for (int u = 0; u < T.N; ++u) {
        auto row = T[u];
        std::fprintf(f, "%d", (int)row.size());
        for (int v : row) std::fprintf(f, " %d", v);
        std::fprintf(f, "\n");
    }
    std::fprintf(f, "\n");
    std::fflush(f);
}

// One isomer outcome under the chosen mode.  In --no-alex we only run
// sorted_dual + dual_idt; the Status-style fields are filled to mirror
// the same shape as the full-pipeline path.
struct Outcome {
    ep::Stage   stage;
    std::string why;
    double      us;
};

static Outcome run_one(const Triangulation& T, bool with_alex) {
    Outcome o;
    auto t0 = clk::now();
    if (with_alex) {
        ep::DualGeometry R = ep::dual_geometry(T);
        o.stage = R.status.stage;
        o.why   = std::move(R.status.why);
    } else {
        try {
            (void)ep::dual_idt(ep::sorted_dual(T));
            o.stage = ep::Stage::OK;
        } catch (const ep::StageError& e) {
            o.stage = e.stage;
            o.why   = e.what();
        } catch (const std::exception& e) {
            o.stage = ep::Stage::UNEXPECTED;
            o.why   = std::string("unexpected: ") + e.what();
        }
    }
    auto t1 = clk::now();
    o.us = std::chrono::duration<double, std::micro>(t1 - t0).count();
    // Strip newlines from why so per-failure records stay one-line.
    for (char& c : o.why) if (c == '\n' || c == '\r') c = ' ';
    return o;
}

int main(int argc, char** argv) {
    int max_N = 200;
    int start_N = 20;
    std::string dir = "validate_results";
    bool with_alex = true;

    int ai = 1;
    while (ai < argc && argv[ai][0] == '-') {
        if (std::strcmp(argv[ai], "-r") == 0 && ai + 1 < argc) {
            start_N = std::atoi(argv[++ai]);
        } else if (std::strcmp(argv[ai], "--no-alex") == 0) {
            with_alex = false;
        } else {
            std::fprintf(stderr, "unknown flag: %s\n", argv[ai]);
            return 1;
        }
        ++ai;
    }
    if (ai < argc) max_N = std::atoi(argv[ai++]);
    if (ai < argc) dir   = argv[ai++];
    fs::create_directories(dir);

    const char* mode_s = with_alex ? "with_alexandrov" : "no_alexandrov";
    const char* fmode  = (start_N > 20) ? "a" : "w";
    FILE* master = std::fopen((dir + "/summary.txt").c_str(), fmode);
    if (start_N <= 20) {
        std::fprintf(master, "# mode: %s\n", mode_s);
        std::fprintf(master,
            "%-6s %14s %10s %10s %10s %10s %10s %10s "
            "%10s %12s %12s %12s %10s\n",
            "N", "isomers", "ok",
            "fail_idt", "fail_alex", "fail_nsim",
            "fail_embed", "fail_cov", "fail_other",
            "med_us", "max_us", "mean_us", "wall_s");
        std::fflush(master);
    }

    for (int N = start_N; N <= max_N; N += 2) {
        if (N == 22) continue;

        auto wall0 = clk::now();
        auto Q = BuckyGen::start(N, false, false);

        long long n_ok = 0;
        long long fail_idt = 0, fail_alex = 0, fail_nsim = 0;
        long long fail_embed = 0, fail_cov = 0, fail_other = 0;
        std::vector<double> us;

        struct FailRec { long long idx; const char* stage_s; std::string why; };
        std::vector<FailRec> fails;
        FILE* graph_dump = nullptr;
        char gpath[256];
        std::snprintf(gpath, sizeof gpath, "%s/C%03d_fails.tri", dir.c_str(), N);

        std::vector<Triangulation> Tbuf(BATCH);
        long long total_idx = 0;

        while (true) {
            int n_read = 0;
            for (int i = 0; i < BATCH; ++i) {
                if (!BuckyGen::next_fullerene(Q, Tbuf[i])) break;
                ++n_read;
            }
            if (n_read == 0) break;

            std::vector<Outcome> outs(n_read);
            #pragma omp parallel for schedule(dynamic, 16)
            for (int i = 0; i < n_read; ++i)
                outs[i] = run_one(Tbuf[i], with_alex);

            for (int i = 0; i < n_read; ++i) {
                if (outs[i].stage == ep::Stage::OK) {
                    ++n_ok;
                    us.push_back(outs[i].us);
                    continue;
                }
                const long long idx = total_idx + i;
                switch (outs[i].stage) {
                    case ep::Stage::IDT_COMPUTE:    ++fail_idt;   break;
                    case ep::Stage::ALEXANDROV:     ++fail_alex;  break;
                    case ep::Stage::NON_SIMPLICIAL: ++fail_nsim;  break;
                    case ep::Stage::EMBED:          ++fail_embed; break;
                    case ep::Stage::COVERAGE:       ++fail_cov;   break;
                    case ep::Stage::INTERPOLATE:    ++fail_embed; break;
                    case ep::Stage::UNEXPECTED:     ++fail_other; break;
                    case ep::Stage::OK:             break;
                }
                const char* sname = ep::stage_name(outs[i].stage);
                fails.push_back({idx, sname, outs[i].why});
                if (!graph_dump) graph_dump = std::fopen(gpath, "w");
                dump_triangulation(graph_dump, Tbuf[i], idx, N, sname);
            }
            total_idx += n_read;

            if (total_idx % 100000 == 0) {
                std::fprintf(stderr,
                    "  C%d: %lld tested, ok=%lld idt=%lld alex=%lld "
                    "nsim=%lld embed=%lld cov=%lld\n",
                    N, total_idx, n_ok,
                    fail_idt, fail_alex, fail_nsim, fail_embed, fail_cov);
            }
        }
        BuckyGen::stop(Q);
        if (graph_dump) std::fclose(graph_dump);

        auto wall1 = clk::now();
        const double wall_s = std::chrono::duration<double>(wall1 - wall0).count();

        auto stat_med = [](std::vector<double>& v) {
            if (v.empty()) return 0.0;
            std::sort(v.begin(), v.end());
            return v[v.size() / 2];
        };
        auto stat_max = [](const std::vector<double>& v) {
            return v.empty() ? 0.0 : *std::max_element(v.begin(), v.end());
        };
        auto stat_mean = [](const std::vector<double>& v) {
            return v.empty() ? 0.0
                : std::accumulate(v.begin(), v.end(), 0.0) / v.size();
        };
        const double max_us  = stat_max (us);
        const double mean_us = stat_mean(us);
        const double med_us  = stat_med (us);

        const long long n_tested = total_idx;

        // Per-N summary file.
        char path[256];
        std::snprintf(path, sizeof path, "%s/C%03d.txt", dir.c_str(), N);
        FILE* fpn = std::fopen(path, "w");
        std::fprintf(fpn, "N           %d\n",   N);
        std::fprintf(fpn, "mode        %s\n",   mode_s);
        std::fprintf(fpn, "isomers     %lld\n", n_tested);
        std::fprintf(fpn, "ok          %lld\n", n_ok);
        std::fprintf(fpn, "fail_idt    %lld\n", fail_idt);
        std::fprintf(fpn, "fail_alex   %lld\n", fail_alex);
        std::fprintf(fpn, "fail_nsim   %lld\n", fail_nsim);
        std::fprintf(fpn, "fail_embed  %lld\n", fail_embed);
        std::fprintf(fpn, "fail_cov    %lld\n", fail_cov);
        std::fprintf(fpn, "fail_other  %lld\n", fail_other);
        std::fprintf(fpn, "med_us      %.3f\n", med_us);
        std::fprintf(fpn, "max_us      %.3f\n", max_us);
        std::fprintf(fpn, "mean_us     %.3f\n", mean_us);
        std::fprintf(fpn, "wall_s      %.3f\n", wall_s);
        #if defined(_OPENMP)
        std::fprintf(fpn, "threads     %d\n",   omp_get_max_threads());
        #endif
        std::fclose(fpn);

        if (!fails.empty()) {
            char fpath[256];
            std::snprintf(fpath, sizeof fpath, "%s/C%03d_fails.txt",
                          dir.c_str(), N);
            FILE* ff = std::fopen(fpath, "w");
            for (const auto& r : fails)
                std::fprintf(ff, "%lld %s %s\n",
                             r.idx, r.stage_s, r.why.c_str());
            std::fclose(ff);
        }

        std::fprintf(master,
            "C%-5d %14lld %10lld %10lld %10lld %10lld %10lld %10lld "
            "%10lld %12.1f %12.1f %12.1f %10.1f\n",
            N, n_tested, n_ok,
            fail_idt, fail_alex, fail_nsim,
            fail_embed, fail_cov, fail_other,
            med_us, max_us, mean_us, wall_s);
        std::fflush(master);

        std::fprintf(stderr,
            "C%d: %lld isomers, ok=%lld idt=%lld alex=%lld nsim=%lld "
            "embed=%lld cov=%lld, %.1fus med, %.1fs wall\n",
            N, n_tested, n_ok,
            fail_idt, fail_alex, fail_nsim, fail_embed, fail_cov,
            med_us, wall_s);
    }

    std::fclose(master);
    return 0;
}
