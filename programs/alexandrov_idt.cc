// alexandrov_idt -- the Alexandrov polytope of one fullerene isomer,
// reconstructed from its canonical generalized-spiral name and written as a
// .ply mesh.
//
// Usage: alexandrov_idt "<canonical name>" [cubic|dual] [output.ply]
//
//   canonical name  e.g. "C60-[GS:1,7,9,11,13,15,18,20,22,24,26,32]-fullerene"
//                   (the human-readable C<N>- prefix is optional)
//   cubic | dual    which metric to embed (default: cubic).  These are
//                   reserved keywords; the two trailing positionals (metric,
//                   output path) may appear in either order.
//   output.ply      defaults to the sanitized name + "-<metric>.ply" in the
//                   current directory
//
// Two metrics, one output pathway (positions + T̄(0) polygonal cells + status):
//
//   cubic  -- AlexandrovIDTCubic on the CUBIC polyhedral metric: every cubic
//             face is a flat regular unit polygon, curvature at the 20..60
//             pentagon-incident cubic vertices.  Closed forms: C20 -> regular
//             dodecahedron (12 pentagons), C60-Ih -> truncated icosahedron
//             (32 faces).  Reports the flat-face census of T̄(0) against the
//             cubic graph's face lattice (pentagons/hexagons realized flat).
//
//   dual   -- AlexandrovSolver on the 12-cone equilateral dual metric: one
//             unit equilateral triangle per cubic vertex, the hexagon-removed
//             iDT leaving the 12 π/3 cones.  Reports the cell count and
//             whether T̄(0) is fully simplicial (n_cells == 2·12 − 4 == 20);
//             the pentagon/hexagon census does not apply.  Closed form: C20 ->
//             regular icosahedron.
//
// The mesh is the polytope with its TRUE 2-faces: vertices are the cone
// points, and each PLY face is one polygonal cell of T̄(0) -- flat regions
// appear as single polygons, not as an arbitrary triangulation (PLY supports
// n-gon faces; any mesh viewer renders them).  Cells are CCW seen from
// outside: the tesselation walk follows the DCEL's CCW orientation and the
// reconstruction enforces positive signed volume.
//
// The PLY is written directly rather than through Polyhedron::to_file: that
// path needs a Polyhedron (a PlanarGraph + points), and the polytope graph
// object does not exist here -- constructing one from the cells would mean
// deriving oriented neighbour rings after the fact, which the orientation
// invariant forbids.  This tool serializes exactly what solve_polytope
// returns: positions + tesselation.
//
// On a failed validation the mesh is still written when positions exist
// (failed cases must be inspectable); the status is printed and reflected in
// the exit code: 0 = OK, 1 = solved-but-invalid or no positions,
// 2 = bad usage / unparsable name.

#include "fullerenes/delaunay_alexandrov.hh"
#include "fullerenes/delaunay.hh"
#include "fullerenes/spiral.hh"
#include "fullerenes/triangulation.hh"

#include <cctype>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;

namespace {

enum class Metric { CUBIC, DUAL };
const char* metric_name(Metric m) { return m == Metric::DUAL ? "dual" : "cubic"; }

// Strip the human-readable "C<N>-" prefix if present; the bracket form is
// what spiral_nomenclature parses (parent CLAUDE.md, canonical names).
string strip_carbon_prefix(const string& name)
{
    if (name.size() > 2 && name[0] == 'C' && isdigit((unsigned char)name[1])) {
        size_t i = 1;
        while (i < name.size() && isdigit((unsigned char)name[i])) i++;
        if (i + 1 < name.size() && name[i] == '-' && name[i + 1] == '[')
            return name.substr(i + 1);
    }
    return name;
}

// Filesystem-safe default output stem: alnum kept, every other run of
// characters collapsed to one underscore.
string sanitize(const string& name)
{
    string s;
    bool gap = false;
    for (char c : name) {
        if (isalnum((unsigned char)c)) {
            if (gap && !s.empty()) s += '_';
            gap = false;
            s += c;
        } else gap = true;
    }
    return s;
}

void usage(FILE* out, const char* prog)
{
    fprintf(out,
        "usage: %s \"<canonical name>\" [cubic|dual] [output.ply]\n"
        "  metric selector defaults to cubic; the two trailing positionals\n"
        "  (metric, output path) may appear in either order.\n"
        "  e.g. %s \"C60-[GS:1,7,9,11,13,15,18,20,22,24,26,32]-fullerene\"\n"
        "       %s \"C20-[GS:1,2,3,4,5,6,7,8,9,10,11,12]-fullerene\" dual\n",
        prog, prog, prog);
}

// Serialize positions + T̄(0) polygonal cells as ASCII PLY.  `census_comment`
// is the mode-specific status line placed as the second PLY comment.
bool write_ply(const string& path, Metric metric, const string& name,
               const AlexandrovSolver::AlexandrovPolytope& P,
               const string& census_comment)
{
    FILE* f = fopen(path.c_str(), "w");
    if (!f) return false;
    fprintf(f, "ply\nformat ascii 1.0\n");
    fprintf(f, "comment %s Alexandrov polytope of %s\n", metric_name(metric), name.c_str());
    fprintf(f, "comment %s\n", census_comment.c_str());
    fprintf(f, "element vertex %zu\n", P.positions.size());
    fprintf(f, "property double x\nproperty double y\nproperty double z\n");
    fprintf(f, "element face %d\n", P.tesselation.n_cells());
    fprintf(f, "property list uchar int vertex_indices\n");
    fprintf(f, "end_header\n");
    for (const coord3d& p : P.positions)
        fprintf(f, "%.17g %.17g %.17g\n", p[0], p[1], p[2]);
    for (const auto& cell : P.tesselation.cells) {
        fprintf(f, "%zu", cell.size());
        for (const auto& [cone, len2] : cell) fprintf(f, " %d", cone);
        fprintf(f, "\n");
    }
    bool ok = (ferror(f) == 0);
    fclose(f);
    return ok;
}

// The B-I solver stats line, identical for both metrics.
void print_solver_stats(const AlexandrovSolver& S,
                        const AlexandrovSolver::AlexandrovPolytope& P)
{
    printf("cones %d  status %s  steps %d  flips %d  newton %d  kappa %.2e\n",
           S.D.nv, AlexandrovSolver::status_str(P.status),
           S.stats_steps, S.stats_flips, S.stats_newton_total, S.stats_final_kappa);
}

}  // namespace

int main(int argc, char** argv)
{
    // Parse: one required name, an optional "cubic"/"dual" keyword, an
    // optional output path.  Metric keyword and output path are order-free.
    string given, out_path;
    Metric metric = Metric::CUBIC;
    bool have_name = false, have_metric = false, have_out = false;

    for (int i = 1; i < argc; i++) {
        const string a = argv[i];
        if (a == "-h" || a == "--help") { usage(stdout, argv[0]); return 0; }
        if (a == "cubic" || a == "dual") {
            if (have_metric) { fprintf(stderr, "metric given twice\n\n");
                               usage(stderr, argv[0]); return 2; }
            metric = (a == "dual") ? Metric::DUAL : Metric::CUBIC;
            have_metric = true;
        } else if (!have_name) {
            given = a; have_name = true;
        } else if (!have_out) {
            out_path = a; have_out = true;
        } else {
            fprintf(stderr, "unexpected extra argument: %s\n\n", a.c_str());
            usage(stderr, argv[0]);
            return 2;
        }
    }
    if (!have_name) { usage(stderr, argv[0]); return 2; }

    Triangulation T;
    try {
        spiral_nomenclature sn(strip_carbon_prefix(given));
        T = Triangulation(sn);
    } catch (const exception& e) {
        fprintf(stderr, "cannot reconstruct isomer from name: %s\n", e.what());
        return 2;
    }
    const int N = 2 * T.N - 4;
    const string name = "C" + to_string(N) + "-" + strip_carbon_prefix(given);
    const string path = have_out ? out_path
                                  : sanitize(name) + "-" + metric_name(metric) + ".ply";

    // Solve in the requested metric.  Both metrics share the output pathway
    // (positions + T̄(0) cells + status); they differ only in the solver they
    // drive and the mode-specific census reported.  `Sp` points at whichever
    // solver ran, for the shared stats line + PLY writing.
    AlexandrovIDTCubic AC;           // storage for the cubic path
    AlexandrovSolver   S_dual;       // storage for the dual path
    AlexandrovSolver::AlexandrovPolytope P;
    const AlexandrovSolver* Sp = nullptr;
    string summary_line;             // stdout T̄(0) summary
    string census_comment;           // PLY status comment
    char buf[256];

    if (metric == Metric::CUBIC) {
        P = AC.solve_polytope(T);
        Sp = &AC.solver;
        const auto census = AC.flat_face_census(T, P.tesselation);
        snprintf(buf, sizeof buf,
                 "status %s  cells %d  pent_flat %d/12  hex_flat %d/%d",
                 AlexandrovSolver::status_str(P.status),
                 census.n_cells, census.pent_flat, census.hex_flat, census.n_hex);
        census_comment = buf;
        snprintf(buf, sizeof buf,
                 "Tbar(0): %d cells, pent_flat %d/12, hex_flat %d/%d%s",
                 census.n_cells, census.pent_flat, census.hex_flat, census.n_hex,
                 census.face_lattice() ? "  (face lattice: every face flat)" : "");
        summary_line = buf;
    } else {
        // 12-cone equilateral dual metric: remove the flat (hexagon) vertices,
        // leaving the 12 π/3 cones, then run the same B-I solver.
        DelaunayTriangulation D = DelaunayTriangulation::compute(T);
        S_dual.D = std::move(D);
        P = S_dual.solve_polytope();
        Sp = &S_dual;
        const int n_cells = P.tesselation.n_cells();
        const int n_simplex = 2 * 12 - 4;   // 20: a fully-triangulated icosahedron
        const bool simplicial = (n_cells == n_simplex);
        snprintf(buf, sizeof buf, "status %s  cells %d/%d  simplicial %s",
                 AlexandrovSolver::status_str(P.status), n_cells, n_simplex,
                 simplicial ? "yes" : "no");
        census_comment = buf;
        snprintf(buf, sizeof buf, "Tbar(0): %d cells (simplicial target %d)%s",
                 n_cells, n_simplex,
                 simplicial ? "  (fully simplicial)" : "  (has flat 2-faces)");
        summary_line = buf;
    }

    printf("%s  [%s metric]\n", name.c_str(), metric_name(metric));
    print_solver_stats(*Sp, P);
    printf("%s\n", summary_line.c_str());

    if (P.positions.empty()) {
        fprintf(stderr, "no positions (status %s) -- nothing to write\n",
                AlexandrovSolver::status_str(P.status));
        return 1;
    }
    if (!write_ply(path, metric, name, P, census_comment)) {
        fprintf(stderr, "FAILED to write %s\n", path.c_str());
        return 1;
    }
    printf("wrote %s  (%zu vertices, %d faces)\n",
           path.c_str(), P.positions.size(), P.tesselation.n_cells());
    return P.ok() ? 0 : 1;
}
