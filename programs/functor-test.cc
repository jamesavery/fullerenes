#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <fullerenes/buckygen-wrapper.hh>
#include <fullerenes/fullerenegraph.hh>
#include <fullerenes/graph.hh>
#include <fullerenes/kernel-headers/dualize-functor.hh>
#include <fullerenes/kernel-headers/forcefield-optimize-functor.hh>
#include <fullerenes/kernel-headers/spherical-projection-functor.hh>
#include <fullerenes/kernel-headers/tutte-functor.hh>
#include <fullerenes/polyhedron.hh>
#include <fullerenes/sycl-headers/sycl-fullerene-structs.hh>
#include <fullerenes/sycl-headers/fill.hh>

namespace {

using real_t = double;
using index_t = std::uint16_t;

constexpr int kBatchSize = 2;

std::string output_name(const std::string& prefix, int n, int index)
{
    return prefix + "-C" + std::to_string(n) + "-iso" + std::to_string(index) + ".mol2";
}


} // namespace

int main(int argc, char** argv)
{
    try {
        const int n = (argc > 1) ? std::stoi(argv[1]) : 60;
        const std::string prefix = (argc > 2) ? argv[2] : "functor-test";

        SyclQueue queue = SyclQueue();
        FullereneBatch<real_t, index_t> batch(n, kBatchSize);
        fill(batch);

        auto batch_view = static_cast<FullereneBatchView<real_t, index_t>>(batch);

        DualizeFunctor<real_t, index_t> dual_functor = DualizeFunctor<real_t, index_t>();
        TutteFunctor<real_t, index_t> tutte_functor = TutteFunctor<real_t, index_t>();
        SphericalProjectionFunctor<real_t, index_t> spherical_project_functor = SphericalProjectionFunctor<real_t, index_t>();
        ForcefieldOptimizeFunctor<ForcefieldType::PEDERSEN, real_t, index_t> optimize_functor = ForcefieldOptimizeFunctor<ForcefieldType::PEDERSEN, real_t, index_t>();

        dual_functor(queue, batch_view, LaunchPolicy::SYNC);
        tutte_functor(queue, batch_view, LaunchPolicy::SYNC);
        spherical_project_functor(queue, batch_view, LaunchPolicy::SYNC);
        optimize_functor(queue, batch_view, LaunchPolicy::SYNC, 5*n, 5*n);

        for (int i = 0; i < batch.size(); ++i) {
            Polyhedron polyhedron = static_cast<Polyhedron>(batch[i]);
            const std::string path = output_name(prefix, n, i);
            if (!Polyhedron::to_file(polyhedron, path)) {
                throw std::runtime_error("failed to write mol2 file: " + path);
            }
            std::cout << "wrote " << path << '\n';
        }

        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "fatal: " << ex.what() << '\n';
        return 1;
    }
}