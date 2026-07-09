// Replays saved extension paths and reports per-step quality at each sub-phase.
//
// Output (stderr): per-isomer summary (angle relerr, edge cv, h_min, concavity, time).
// Output (stdout): CSV with per-step data at placed/reflected/patched/relaxed phases.
//
// Usage: replay_perstep [--limit N] [--mol2-step S]... paths1.json [paths2.json ...]
//   --limit N: only process first N isomers per file (default: all)
//   --mol2-step S: dump mol2 + per-vertex detail at step S (repeatable, -1 = auto-detect worst)

#include "fullerenes/buckinverse.hh"
#include "fullerenes/deltahedron.hh"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/stats.hh"
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <chrono>
#include <set>

using namespace std;
using namespace buckinverse;

// Minimal JSON parser
struct JsonVal {
    enum Type { Null, Num, Str, Arr, Obj } type = Null;
    double num = 0;
    string str;
    // Only vector<JsonVal> may hold the incomplete JsonVal; a composite such as
    // vector<pair<string,JsonVal>> is UB. Obj keeps its keys alongside in `keys`.
    vector<JsonVal> arr;
    vector<string> keys;
    const JsonVal& operator[](const string& key) const {
        for (size_t i = 0; i < keys.size(); i++) if (keys[i] == key) return arr[i];
        static JsonVal null_val; return null_val;
    }
    const JsonVal& operator[](size_t i) const { return arr[i]; }
    int as_int() const { return (int)num; }
    double as_double() const { return num; }
    const string& as_str() const { return str; }
    size_t size() const { return type == Arr || type == Obj ? arr.size() : 0; }
};

static void skip_ws(const char*& p) { while (*p == ' ' || *p == '\n' || *p == '\r' || *p == '\t') p++; }

static JsonVal parse_json(const char*& p) {
    skip_ws(p);
    JsonVal v;
    if (*p == '"') {
        v.type = JsonVal::Str; p++;
        while (*p && *p != '"') {
            if (*p == '\\') { p++; v.str += *p++; }
            else v.str += *p++;
        }
        if (*p == '"') p++;
        return v;
    }
    if (*p == '[') {
        v.type = JsonVal::Arr; p++; skip_ws(p);
        if (*p == ']') { p++; return v; }
        while (true) {
            v.arr.push_back(parse_json(p)); skip_ws(p);
            if (*p == ',') { p++; continue; }
            if (*p == ']') { p++; break; }
            break;
        }
        return v;
    }
    if (*p == '{') {
        v.type = JsonVal::Obj; p++; skip_ws(p);
        if (*p == '}') { p++; return v; }
        while (true) {
            auto key = parse_json(p); skip_ws(p);
            if (*p == ':') p++;
            auto val = parse_json(p);
            v.keys.push_back(key.str); v.arr.push_back(val); skip_ws(p);
            if (*p == ',') { p++; skip_ws(p); continue; }
            if (*p == '}') { p++; break; }
            break;
        }
        return v;
    }
    v.type = JsonVal::Num;
    char* end; v.num = strtod(p, &end); p = end;
    return v;
}

static ExtensionPath json_to_ep(const JsonVal& bp) {
    ExtensionPath ep;
    string seed = bp["seed"].as_str();
    ep.seed = seed == "C20" ? SeedType::C20 : seed == "C28" ? SeedType::C28 : SeedType::C30;
    ep.full_N = bp["full_N"].as_int();
    auto& ss = bp["seed_state"];
    for (size_t i = 0; i < ss.size(); i++) {
        ExtensionPath::SeedVertex sv;
        sv.id = ss[i]["id"].as_int();
        auto& nbr = ss[i]["nbr"];
        for (int j = 0; j < 6; j++) sv.nbr[j] = nbr[j].as_int();
        sv.active = ss[i]["active"].as_int();
        ep.seed_state.push_back(sv);
    }
    auto& steps = bp["steps"];
    for (size_t i = 0; i < steps.size(); i++) {
        ExtensionStep st;
        string kind_str = steps[i]["kind"].as_str();
        if (kind_str[0] == 'L') { st.kind.type = ExpKind::L_type; st.kind.i = stoi(kind_str.substr(1)); }
        else if (kind_str[0] == 'B') {
            st.kind.type = ExpKind::B_type;
            auto paren = kind_str.find('('); auto comma = kind_str.find(','); auto rparen = kind_str.find(')');
            st.kind.i = stoi(kind_str.substr(paren+1, comma-paren-1));
            st.kind.j = stoi(kind_str.substr(comma+1, rparen-comma-1));
        } else { st.kind.type = ExpKind::F_type; }
        st.dir = steps[i]["dir"].as_str() == "R" ? Dir::DRight : Dir::DLeft;
        auto& strip = steps[i]["strip"];
        for (size_t j = 0; j < strip.size(); j++) st.strip.push_back(strip[j].as_int());
        auto& path = steps[i]["path"];
        for (size_t j = 0; j < path.size(); j++) st.path.push_back(path[j].as_int());
        auto& tp = steps[i]["tp"];
        for (size_t j = 0; j < tp.size(); j++) st.tp.push_back(tp[j].as_int());
        ep.steps.push_back(st);
    }
    return ep;
}

static string step_kind_str(const ExpKind& k) {
    if (k.type == ExpKind::L_type) return "L" + to_string(k.i);
    if (k.type == ExpKind::B_type) return "B(" + to_string(k.i) + "," + to_string(k.j) + ")";
    return "F";
}

static void write_mol2(const Deltahedron& D, const string& path) {
    Polyhedron P(static_cast<const PlanarGraph&>(D), vector<coord3d>(D.points.begin(), D.points.end()));
    Polyhedron::to_file(P, path);
}

// Print per-vertex detail for concave and worst-h vertices.
// strip_compact: set of compact vertex IDs that are strip vertices.
// path_compact/tp_compact: similar for path/tp.
static void print_vertex_detail(const Deltahedron& D, int N, int bucky_idx, int step,
                                const char* phase,
                                const set<int>& strip_compact,
                                const set<int>& path_compact,
                                const set<int>& tp_compact) {
    // Compute h for all vertices
    struct VH { int id; double h; };
    vector<VH> all_vh;
    for (int v = 0; v < D.N; v++) {
        int d = D.degree(v);
        if (d > 6) continue;
        coord3d cen(0,0,0);
        for (int i = 0; i < d; i++) cen += D.points[D[v][i]];
        cen /= (double)d;
        coord3d nf(0,0,0);
        for (int i = 0; i < d; i++) {
            coord3d e1 = D.points[D[v][i]] - D.points[v];
            coord3d e2 = D.points[D[v][(i+1)%d]] - D.points[v];
            nf += e1.cross(e2);
        }
        double nl = nf.norm();
        if (nl < 1e-15) continue;
        double h = (D.points[v] - cen).dot(nf / nl);
        all_vh.push_back({v, h});
    }
    sort(all_vh.begin(), all_vh.end(), [](const VH& a, const VH& b){ return a.h < b.h; });

    auto classify = [&](int v) -> const char* {
        if (strip_compact.count(v)) return "STRIP";
        if (path_compact.count(v))  return "PATH";
        if (tp_compact.count(v))    return "TP";
        return "other";
    };

    // Print concave vertices
    int n_concave = 0;
    for (auto& vh : all_vh) if (vh.h < 0) n_concave++;

    fprintf(stderr, "  === C%d #%d step %d [%s]: %d concave, worst h = %+.6f ===\n",
            N, bucky_idx, step, phase, n_concave, all_vh.empty() ? 0.0 : all_vh[0].h);

    // Print all concave + worst 5
    int printed = 0;
    for (auto& vh : all_vh) {
        if (vh.h >= 0 && printed >= 5) break;
        int v = vh.id;
        int d = D.degree(v);
        // Edge lengths to neighbors
        fprintf(stderr, "    v%-3d deg%d %-5s h=%+.6f  edges:", v, d, classify(v), vh.h);
        for (int i = 0; i < d; i++) {
            int nb = D[v][i];
            double elen = (D.points[v] - D.points[nb]).norm();
            fprintf(stderr, " %.4f(%s)", elen, classify(nb));
        }
        fprintf(stderr, "\n");
        printed++;
    }
}

// Quality snapshot at a sub-phase.
struct PhaseSnap {
    double ang_re;      // max angle relative error (global)
    double ang_local;   // max angle relative error at strip-adjacent faces only
    double edge_cv;     // edge length CV
    double h_min;       // minimum convexity height
    int n_concave;      // number of concave vertices
    int iters;          // optimizer iterations used (meaningful for patched/cg/final)
};

// Compute the face angle relative error |theta - pi/3| / (pi/3) for a single triangle.
static double face_angle_relerr(const coord3d& a, const coord3d& b, const coord3d& c) {
    const double pi3 = M_PI / 3.0;
    double max_re = 0;
    // Three angles: at a, at b, at c
    coord3d ab = b - a, ac = c - a;
    double cos_a = ab.dot(ac) / (ab.norm() * ac.norm());
    cos_a = max(-1.0, min(1.0, cos_a));
    max_re = max(max_re, fabs(acos(cos_a) - pi3) / pi3);

    coord3d ba = a - b, bc = c - b;
    double cos_b = ba.dot(bc) / (ba.norm() * bc.norm());
    cos_b = max(-1.0, min(1.0, cos_b));
    max_re = max(max_re, fabs(acos(cos_b) - pi3) / pi3);

    coord3d ca = a - c, cb = b - c;
    double cos_c = ca.dot(cb) / (ca.norm() * cb.norm());
    cos_c = max(-1.0, min(1.0, cos_c));
    max_re = max(max_re, fabs(acos(cos_c) - pi3) / pi3);
    return max_re;
}

// Compute max angle relerr of faces touching any vertex in strip_set.
static double strip_local_angle_relerr(const Deltahedron& D, const set<int>& strip_set) {
    if (strip_set.empty()) return 0;
    double max_re = 0;
    for (int u = 0; u < D.N; u++) {
        int d = D.degree(u);
        for (int i = 0; i < d; i++) {
            int v = D[u][i];
            int w = D[u][(i + 1) % d];
            if (v < u || w < u) continue; // avoid counting each face 3x
            // Check if any vertex of this face is in strip_set
            if (strip_set.count(u) || strip_set.count(v) || strip_set.count(w)) {
                double re = face_angle_relerr(D.points[u], D.points[v], D.points[w]);
                if (re > max_re) max_re = re;
            }
        }
    }
    return max_re;
}

static PhaseSnap snap_quality(const Deltahedron& D, const set<int>& strip_compact = {}) {
    PhaseSnap s{};
    s.ang_re = D.max_angle_relerr();
    s.ang_local = strip_compact.empty() ? s.ang_re : strip_local_angle_relerr(D, strip_compact);
    s.n_concave = D.count_concave();

    // Edge CV
    vector<double> lens;
    for (int u = 0; u < D.N; u++)
        for (int v : D[u])
            if (v > u) lens.push_back((D.points[u] - D.points[v]).norm());
    s.edge_cv = cv_twopass(lens);

    // h_min
    s.h_min = 1e30;
    for (int v = 0; v < D.N; v++) {
        int d = D.degree(v);
        if (d > 6) continue;
        coord3d cen(0,0,0);
        for (int i = 0; i < d; i++) cen += D.points[D[v][i]];
        cen /= (double)d;
        coord3d nf(0,0,0);
        for (int i = 0; i < d; i++) {
            coord3d e1 = D.points[D[v][i]] - D.points[v];
            coord3d e2 = D.points[D[v][(i+1)%d]] - D.points[v];
            nf += e1.cross(e2);
        }
        double nl = nf.norm();
        if (nl < 1e-15) continue;
        double h = (D.points[v] - cen).dot(nf / nl);
        if (h < s.h_min) s.h_min = h;
    }
    s.iters = D.iterations_used;
    return s;
}

// Per-step data for one method run.
struct StepRecord {
    int step;
    string kind;
    int Nv;
    int n_strip;        // number of strip vertices in this step
    PhaseSnap placed, reflected, patched, relaxed;
    bool complete = false;
};

int main(int argc, char** argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s [--limit N] [--mol2-step S]... paths1.json [paths2.json ...]\n", argv[0]);
        return 1;
    }

    int limit = 0;  // 0 = all
    set<int> mol2_steps;  // steps to dump mol2 for (-1 = auto-detect worst)
    vector<string> files;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--limit") && i+1 < argc) {
            limit = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--mol2-step") && i+1 < argc) {
            mol2_steps.insert(atoi(argv[++i]));
        } else {
            files.push_back(argv[i]);
        }
    }

    // CSV header
    printf("method,size,idx,step,kind,Nv,n_strip,"
           "ang_placed,ang_reflected,ang_patched,ang_relaxed,"
           "local_placed,local_reflected,local_patched,local_relaxed,"
           "ecv_placed,ecv_reflected,ecv_patched,ecv_relaxed,"
           "hmin_placed,hmin_reflected,hmin_patched,hmin_relaxed,"
           "conc_placed,conc_reflected,conc_patched,conc_relaxed,"
           "patch_iters,relax_iters\n");

    for (auto& fname : files) {
        ifstream f(fname);
        string content((istreambuf_iterator<char>(f)), istreambuf_iterator<char>());
        f.close();
        const char* p = content.c_str();
        JsonVal root = parse_json(p);
        int n_isomers = (int)root.size();
        if (limit > 0 && n_isomers > limit) n_isomers = limit;
        fprintf(stderr, "Processing %s (%d isomers)\n", fname.c_str(), n_isomers);

        for (int ii = 0; ii < n_isomers; ii++) {
            auto& entry = root[ii];
            int N = entry["N"].as_int();
            int bucky_idx = entry["buckygen_idx"].as_int();
            ExtensionPath ep = json_to_ep(entry["bucky_path"]);

            vector<StepRecord> records;
            int cur_step = -1;

            // Track worst step for auto-detect mode (-1)
            int worst_step = -1;
            double worst_relax_iters = 0;

            auto diag = [&](int step, const char* phase, const Deltahedron& D) {
                string ph(phase);

                if (ph == "seed") return;

                int si = step - 1;
                string kind = (si >= 0 && si < (int)ep.steps.size())
                              ? step_kind_str(ep.steps[si].kind) : "?";

                // Build strip_compact: the strip vertex IDs from ep.steps[si]
                // mapped to compact numbering. Strip vertices have the highest
                // full-graph IDs (they were appended last). In extractCompact,
                // active vertices are enumerated in order 0..full_N-1, so strip
                // vertices get the highest compact IDs = D.N - n_strip .. D.N - 1.
                set<int> strip_compact;
                set<int> path_compact, tp_compact;
                int n_strip = 0;
                if (si >= 0 && si < (int)ep.steps.size()) {
                    n_strip = (int)ep.steps[si].strip.size();
                    for (int k = 0; k < n_strip; k++)
                        strip_compact.insert(D.N - n_strip + k);

                    // Build path/tp compact IDs by scanning the full-graph IDs
                    // through the compact mapping. The compact ID of full-graph
                    // vertex u is its rank among active vertices (sorted by ID).
                    // Since all vertices 0..D.N-1 in the compact graph are active
                    // and ordered, compact(u) = u for vertices below the strip range.
                    // But we need the actual mapping. For simplicity, use the
                    // full-graph IDs directly and check: vertices in ep.steps[si].path
                    // and .tp that are < D.N - n_strip are their own compact IDs
                    // (since the compact extraction preserves order of pre-existing vertices).
                    for (int pid : ep.steps[si].path)
                        if (pid < D.N) path_compact.insert(pid);
                    for (int tid : ep.steps[si].tp)
                        if (tid < D.N) tp_compact.insert(tid);
                }

                PhaseSnap snap = snap_quality(D, strip_compact);
                snap.iters = D.iterations_used;

                // Check if this step should get mol2 + detail
                bool dump_this = mol2_steps.count(si) > 0;

                if (dump_this) {
                    char mol2_path[256];
                    snprintf(mol2_path, sizeof(mol2_path), "/tmp/C%d_%d_step%d_%s.mol2",
                             N, bucky_idx, si, phase);
                    write_mol2(D, mol2_path);
                    print_vertex_detail(D, N, bucky_idx, si, phase,
                                        strip_compact, path_compact, tp_compact);
                }

                if (ph == "placed") {
                    StepRecord r;
                    r.step = si;
                    r.kind = kind;
                    r.Nv = D.N;
                    r.n_strip = n_strip;
                    r.placed = snap;
                    records.push_back(r);
                    cur_step = (int)records.size() - 1;
                } else if (ph == "reflected" && cur_step >= 0) {
                    records[cur_step].reflected = snap;
                } else if (ph == "patched" && cur_step >= 0) {
                    records[cur_step].patched = snap;
                } else if (ph == "relaxed" && cur_step >= 0) {
                    records[cur_step].relaxed = snap;
                    records[cur_step].complete = true;
                    // Track worst step by relaxation iterations
                    if (snap.iters > worst_relax_iters) {
                        worst_relax_iters = snap.iters;
                        worst_step = si;
                    }
                }
            };

            auto t0 = chrono::steady_clock::now();
            Deltahedron D = Deltahedron::fromExtensionPathOptimized(
                ep, nullptr, diag,
                OptMethod::LBFGS, 1e-4, 1e-11, 0, 0, 0,
                OptMethod::STEIHAUG, 1e-10);
            double ms = chrono::duration<double, milli>(chrono::steady_clock::now() - t0).count();

            PhaseSnap final_snap = snap_quality(D);

            fprintf(stderr, "  C%d #%d: ang=%.2e ecv=%.2e h_min=%+.4f conc=%d  %.0fms\n",
                    N, bucky_idx,
                    final_snap.ang_re, final_snap.edge_cv,
                    final_snap.h_min, final_snap.n_concave, ms);

            // Auto-detect: report worst step so user can re-run with --mol2-step N
            if (mol2_steps.count(-1) && worst_step >= 0) {
                // Find the record for worst_step
                for (auto& r : records) {
                    if (r.step == worst_step && r.complete) {
                        fprintf(stderr, "  WORST STEP: %d (%s) Nv=%d relax_iters=%d "
                                "hmin: placed=%+.4f reflected=%+.4f patched=%+.4f relaxed=%+.4f\n",
                                worst_step, r.kind.c_str(), r.Nv, r.relaxed.iters,
                                r.placed.h_min, r.reflected.h_min,
                                r.patched.h_min, r.relaxed.h_min);
                        fprintf(stderr, "  Re-run with: --mol2-step %d\n", worst_step);
                        break;
                    }
                }
            }

            // Output CSV rows for each step
            for (auto& r : records) {
                if (!r.complete) continue;
                printf("epopt,%d,%d,%d,%s,%d,%d,"
                       "%.6e,%.6e,%.6e,%.6e,"
                       "%.6e,%.6e,%.6e,%.6e,"
                       "%.6e,%.6e,%.6e,%.6e,"
                       "%.6e,%.6e,%.6e,%.6e,"
                       "%d,%d,%d,%d,"
                       "%d,%d\n",
                       N, bucky_idx, r.step, r.kind.c_str(), r.Nv, r.n_strip,
                       r.placed.ang_re, r.reflected.ang_re, r.patched.ang_re, r.relaxed.ang_re,
                       r.placed.ang_local, r.reflected.ang_local, r.patched.ang_local, r.relaxed.ang_local,
                       r.placed.edge_cv, r.reflected.edge_cv, r.patched.edge_cv, r.relaxed.edge_cv,
                       r.placed.h_min, r.reflected.h_min, r.patched.h_min, r.relaxed.h_min,
                       r.placed.n_concave, r.reflected.n_concave, r.patched.n_concave, r.relaxed.n_concave,
                       r.patched.iters, r.relaxed.iters);
            }
            fflush(stdout);
        }
    }
    return 0;
}
