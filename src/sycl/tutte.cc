#include <sycl/sycl.hpp>
#include "numeric"
#include <vector>
#include <tuple>
#include <iterator>
#include <type_traits>
#include <algorithm>
#include <fullerenes/sycl-headers/execution-compat.hh>
#include <numeric>
#include <fullerenes/kernel-headers/tutte.hh>
#include "kernel.hh"
#include "queue-impl.hh"
#include "forcefield-includes.hh"

//Global memory arrays needed for the general tutte layout kernel:
//  - XY double buffer: The position of each vertex
//  - Fixed array: Whether the vertex is fixed or not
//  - Reduction array: The maximum change in position of each vertex


// ---------------------------------------------------------------------------
// View-based batch implementation (Phase 7)
// ---------------------------------------------------------------------------
template <typename T, typename K>
static SyclEvent tutte_view_batch_impl(
    SyclQueue& Q,
    batch::BatchView<Spanify::RSRAdjacencyView<K>> graph,
    std::span<std::array<T,2>> layout,
    batch::BatchStateView state)
{
    TEMPLATE_TYPEDEFS(T,K);
    auto [adj_flat, deg_flat, twin_flat] = graph.spans();
    (void)deg_flat; (void)twin_flat;

    // Reinterpret flat K* as array<K,3>* (dmax==3 for cubic)
    std::span<std::array<K,3>> A_cubic(
        reinterpret_cast<std::array<K,3>*>(adj_flat.data()),
        adj_flat.size() / 3);
    // coord2d == std::array<T,2> (see FLOAT_TYPEDEFS macro)
    std::span<coord2d> layout_cd(
        reinterpret_cast<coord2d*>(layout.data()),
        layout.size());

    auto statuses   = state.status;
    const int N        = graph.N();
    const int capacity = graph.size();
    const auto max_iter = (size_t)N * 50;

    SyclEventImpl tutte_done = Q->submit([&](handler& h) {
        local_accessor<bool, 1>    smem(N, h);
        local_accessor<coord2d, 1> xys_smem(N, h);
        local_accessor<coord2d, 1> newxys_smem(N, h);

        h.parallel_for(sycl::nd_range(sycl::range(N*capacity), sycl::range(N)),
        [=](nd_item<1> nditem) {
            const auto cta        = nditem.get_group();
            const auto a          = nditem.get_local_linear_id();
            const auto isomer_idx = nditem.get_group_linear_id();

            if (!statuses[isomer_idx].is_set(StatusEnum::FULLERENEGRAPH_PREPARED)) return;

            const auto isomer_adj = A_cubic.subspan(isomer_idx * N, N);
            auto xys_acc          = layout_cd.subspan(isomer_idx * N, N);

            DeviceCubicGraph FG(isomer_adj);
            node3 ns = FG[a];
            xys_smem[a] = {T(0),T(0)};

            std::array<node_t, 6> outer_face;
            uint8_t Nface = FG.get_face_oriented(0, FG[0][0], outer_face);

            smem[a] = false;
            sycl::group_barrier(cta);
            if (a < Nface) smem[outer_face[a]] = true;
            sycl::group_barrier(cta);
            bool fixed = smem[a];

            real_t phase = real_t(2) * real_t(M_PI) / Nface;
            if (a < Nface) xys_smem[outer_face[a]] = {sycl::sin(a*phase), sycl::cos(a*phase)};
            sycl::group_barrier(cta);

            bool   converged  = false;
            real_t max_change = real_t(0);
            if (fixed) newxys_smem[a] = xys_smem[a];

            for (size_t i = 0; i < max_iter && !converged; i++) {
                max_change = real_t(0);
                sycl::group_barrier(cta);
                coord2d neighbour_sum = {T(0),T(0)};
                for (uint8_t j = 0; j < 3; j++) neighbour_sum += xys_smem[ns[j]];

                if (!fixed) newxys_smem[a] = xys_smem[a]*real_t(0.15) + (neighbour_sum/real_t(3.))*real_t(0.85);
                real_t neighbour_dist = real_t(0);
                for (uint8_t j = 0; j < 3; j++) neighbour_dist += norm(xys_smem[a] - xys_smem[d_get(ns,j)]) / real_t(3);

                sycl::group_barrier(cta);
                real_t relative_change = real_t(0);
                if (neighbour_dist > real_t(0) && !fixed)
                    relative_change = norm(xys_smem[a] - newxys_smem[a]) / neighbour_dist;

                real_t iteration_max = sycl::reduce_over_group(cta, relative_change, sycl::maximum<real_t>());
                if (iteration_max > max_change) max_change = iteration_max;
                converged = max_change <= 10 * std::numeric_limits<real_t>::epsilon();

                xys_smem[a] = newxys_smem[a];
            }
            sycl::group_barrier(cta);
            xys_acc[a] = xys_smem[a];
            if (a == 0 && converged) statuses[isomer_idx] |= StatusEnum::CONVERGED_2D;
        });
    });
    return SyclEvent(std::move(tutte_done));
}

template <typename T, typename K>
SyclEvent tutte_layout(SyclQueue& Q,
                       batch::BatchView<Spanify::RSRAdjacencyView<K>> graph,
                       std::span<std::array<T,2>> layout,
                       batch::BatchStateView state,
                       Workspace /*ws*/) {
    return tutte_view_batch_impl<T,K>(Q, graph, layout, state);
}

template SyclEvent tutte_layout<float, uint16_t>(SyclQueue&, batch::BatchView<Spanify::RSRAdjacencyView<uint16_t>>, std::span<std::array<float,2>>, batch::BatchStateView, Workspace);
template SyclEvent tutte_layout<float, uint32_t>(SyclQueue&, batch::BatchView<Spanify::RSRAdjacencyView<uint32_t>>, std::span<std::array<float,2>>, batch::BatchStateView, Workspace);
template SyclEvent tutte_layout<double,uint16_t>(SyclQueue&, batch::BatchView<Spanify::RSRAdjacencyView<uint16_t>>, std::span<std::array<double,2>>, batch::BatchStateView, Workspace);
template SyclEvent tutte_layout<double,uint32_t>(SyclQueue&, batch::BatchView<Spanify::RSRAdjacencyView<uint32_t>>, std::span<std::array<double,2>>, batch::BatchStateView, Workspace);

// Phase 2: see dualize_buffer_size — returns 0 (local_accessor for scratch).
template <typename T, typename K>
size_t tutte_layout_buffer_size(const Device&, int, int) { return 0; }
template size_t tutte_layout_buffer_size<float, uint16_t>(const Device&, int, int);
template size_t tutte_layout_buffer_size<float, uint32_t>(const Device&, int, int);
template size_t tutte_layout_buffer_size<double,uint16_t>(const Device&, int, int);
template size_t tutte_layout_buffer_size<double,uint32_t>(const Device&, int, int);