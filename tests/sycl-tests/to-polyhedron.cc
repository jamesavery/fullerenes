// Polyhedron <-> batch round-trip via the new low-level batch API.
//
// Replaces the legacy FullereneBatch::push_back(Polyhedron) /
// (Polyhedron)batch[0] cast test. The new API stores adjacency in
// batch::Batch<RSR> and xyz in a parallel SyclVector; we round-trip both
// and re-build a host Polyhedron to verify volume_divergence() is preserved.

#include <iostream>
#include <fullerenes/buckygen-wrapper.hh>
#include <fullerenes/graph.hh>
#include <fullerenes/triangulation.hh>
#include <fullerenes/polyhedron.hh>
#include <fullerenes/batch/batch.hh>
#include <fullerenes/dense_graph.hh>
#include <fullerenes/sycl-headers/sycl-vector.hh>

int main() {
    constexpr int N    = 20;
    using K = uint16_t;
    using RSR = Spanify::RSRAdjacencyView<K>;

    // Build a Polyhedron host-side.
    Graph G(N/2 + 2, GRAPH_DMAX);
    BuckyGen::buckygen_queue BQ = BuckyGen::start(N, false, false);
    BuckyGen::next_fullerene(BQ, G);
    BuckyGen::stop(BQ);
    FullereneDual D(G);
    PlanarGraph PG = D.dual_graph();
    Polyhedron P(PG);
    P.set_points(P.zero_order_geometry());
    std::cout.setstate(std::ios_base::failbit);
    P.optimize();
    std::cout.clear();

    // Push adjacency into a 1-isomer Batch<RSR>(N, 1, dmax=3).
    batch::Batch<RSR> b(N, 1, /*dmax*/3);
    SyclVector<std::array<float,3>> xyz(N);
    {
        auto vc  = b.view_capacity();
        auto& adj = std::get<0>(vc.spans());
        auto& deg = std::get<1>(vc.spans());
        for (int u = 0; u < N; ++u) {
            for (int k = 0; k < 3; ++k)
                adj[u*3 + k] = K(P.nbrs(u)[k]);
            deg[u] = 3;
            xyz[u] = { float(P.points[u][0]),
                       float(P.points[u][1]),
                       float(P.points[u][2]) };
        }
        b.resize(1);
    }

    // Read back.
    Graph G_out(N, true);
    std::vector<coord3d> pts_out(N);
    {
        auto v = b.view()[0];
        auto& adj = std::get<0>(v.to_tuple());
        for (int u = 0; u < N; ++u) {
            G_out.assign_row(u, { adj[u*3], adj[u*3+1], adj[u*3+2] });
            pts_out[u] = { xyz[u][0], xyz[u][1], xyz[u][2] };
        }
    }
    Polyhedron Pout(PlanarGraph(G_out), pts_out);

    std::cout.precision(10);
    std::cout << "in  volume_divergence = " << P.volume_divergence()    << std::endl;
    std::cout << "out volume_divergence = " << Pout.volume_divergence() << std::endl;
    if (std::abs(P.volume_divergence() - Pout.volume_divergence()) > 1e-4) {
        std::cerr << "Volume mismatch after round-trip!" << std::endl;
        return 1;
    }
    return 0;
}
