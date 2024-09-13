#include <iostream>
#include <fullerenes/sycl-headers/sycl-fullerene-structs.hh>
#include <fullerenes/buckygen-wrapper.hh>
#include <fullerenes/graph.hh>
#include <fullerenes/polyhedron.hh>

int main() {
    FullereneDual G(20);
    BuckyGen::buckygen_queue BQ = BuckyGen::start(20, false, false);
    BuckyGen::next_fullerene(BQ, G);
    G.update();
    PlanarGraph PG = G.dual_graph();
    PG.layout2d = PG.tutte_layout();
    Polyhedron P = Polyhedron(PG);
    P.points = P.zero_order_geometry();
    std::cout.setstate(std::ios_base::failbit);
    P.optimize();
    std::cout.clear();
    FullereneBatch<float, uint16_t> batch(20, 1);
    batch.push_back(P);
    Polyhedron Pout = (Polyhedron)batch[0];
    //assert(std::abs(Pout.volume_divergence() - P.volume_divergence()) < 1e-6);
    std::cout.precision(10);
    std::cout << P.volume_divergence() << std::endl;
    std::cout << Pout.volume_divergence() << std::endl;

}