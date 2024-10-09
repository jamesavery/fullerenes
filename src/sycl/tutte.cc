#include <CL/sycl.hpp>
#include "numeric"
#include <vector>
#include <tuple>
#include <iterator>
#include <type_traits>
#include <fullerenes/sycl-headers/tutte-kernel.hh>
#include "kernel.cc"
#include "primitives.cc"
#include "queue-impl.cc"
#include "forcefield-includes.cc"

//Global memory arrays needed for the general tutte layout kernel:
//  - XY double buffer: The position of each vertex
//  - Fixed array: Whether the vertex is fixed or not
//  - Reduction array: The maximum change in position of each vertex

template <typename T, typename K> struct TutteKernel_1 {};

template <typename T, typename K>
SyclEvent tutte_isomer_impl( SyclQueue& Q, 
                        Span<std::array<K,3>> cubic_neighbours,
                        Span<std::array<T,2>> xys,
                        Span<std::array<T,2>> newxys,
                        Span<bool> fixed,
                        Span<T> max_change,
                        int N)
{   
    LaunchConfig<TutteKernel_1<T,K>> config1(Q);

    primitives::fill(Q, fixed, false);
    primitives::fill(Q, xys, std::array<T,2>{0,0});
    primitives::fill(Q, newxys, std::array<T,2>{0,0});

    Q->submit([&](handler& h) { h.parallel_for<TutteKernel_1<T,K>> ( config1.isomer_nd_range(6) , [=](nd_item<1> nditem){
            auto tid = nditem.get_global_linear_id();
            DeviceCubicGraph FG(cubic_neighbours);
            std::array<K, 6> outer_face;
            uint8_t Nface = FG.get_face_oriented(0, FG[0][0], outer_face);
            if(tid < Nface){ fixed[outer_face[tid]] = true;}
            T phase = M_PI*2/Nface;
            if(tid < Nface) xys[outer_face[tid]] = {sycl::sin(tid*phase), sycl::cos(tid*phase)};
    });});

    primitives::copy(Q, xys.subspan(0,N), newxys);

    for (size_t i = 0; i < N*10; i++){
        primitives::transform(Q, cubic_neighbours.template as_span<std::array<K,3>>(), xys, newxys, [=](auto& ns, auto& xy){
            auto idx = &xy - &xys[0];
            std::array<T,2> neighbour_sum = {0,0};
            for (int j = 0; j < 3; j++) neighbour_sum += xys[ns[j]];
            //Weighted neighbour coordinate sum based on the distance between the neighbours
            /* T max_dist = 0;
            K max_idx = 0;
            for (int j = 0; j < 3; j++) {
                T dist = norm(xys[ns[j]] - xy);
                if (dist > max_dist) {max_dist = dist; max_idx = ns[j];}
            }
            for (int j = 0; j < 3; j++) {
                T dist = norm(xys[ns[j]] - xy);
                neighbour_sum += xys[ns[j]]*(dist/max_dist);
            } */
            if(!fixed[idx]) return xy*T(0.15) + neighbour_sum*T(0.85)/T(3);
            return xy;
        });

        if (i%N == 0){
            Q.wait();
            T max_diff = primitives::transform_reduce(Q, xys, newxys, T(0), Max{}, [&](auto& old_xys, auto& new_xys){return norm(old_xys - new_xys);});
            if (max_diff <= 10*std::numeric_limits<T>::epsilon()) break;
        }
        primitives::copy(Q, newxys.subspan(0,N), xys);
    }
    return Q.get_event();
}


template<typename T, typename K>
SyclEvent tutte_batch_impl(SyclQueue& Q, FullereneBatchView<T,K> batch){
    TEMPLATE_TYPEDEFS(T,K);
    auto statuses = batch.m_.flags_;
    SyclEventImpl tutte_done = Q->submit([&](handler& h) {
        const auto N        = batch.N_;
        const auto capacity = batch.capacity();
        const auto max_iter = N * 50;

        local_accessor<bool, 1>     smem(N, h);        // TODO: More transparent name?
        local_accessor<coord2d, 1>  xys_smem(N, h);    
        local_accessor<coord2d, 1>  newxys_smem(N, h);
        h.parallel_for(sycl::nd_range(sycl::range(N*capacity), sycl::range(N)), [=](nd_item<1> nditem) {
 
            const auto   cta               = nditem.get_group();
            const auto   a                 = nditem.get_local_linear_id(); // Atom  vertex index in graph
            const auto   isomer_idx        = nditem.get_group_linear_id(); // Isomer graph index in batch
            auto fullerene                 = batch[isomer_idx];

            const auto& isomer_neighbours  = fullerene.d_.A_cubic_;
            auto xys_acc                  = fullerene.d_.X_cubic_.template as_span<coord2d>();
	    
            if (statuses[isomer_idx] & StatusFlag::FULLERENEGRAPH_PREPARED){

            DeviceCubicGraph FG(isomer_neighbours); 

            node3 ns       = FG[a];
            xys_smem[a]    = {0,0};

            std::array<node_t, 6> outer_face;
            node_t outer_face_vertex   = 0;
            uint8_t Nface = FG.get_face_oriented(0,FG[0][0], outer_face);    
            
            smem[a] =  false; 
            sycl::group_barrier(cta);
            if(a < Nface){
                outer_face_vertex = outer_face[a];
                smem[outer_face_vertex] =  true; 
            }
            sycl::group_barrier(cta);;
            bool fixed = smem[a];

	    real_t phase = 2*real_t(M_PI)/Nface;
            if(a < Nface) xys_smem[outer_face_vertex] = {sycl::sin(a*phase),sycl::cos(a*phase)};

	    sycl::group_barrier(cta);
	    
            bool   converged       = false;
            real_t max_change       = 0;

            if(fixed) newxys_smem[a] = xys_smem[a];
	    
            for (size_t i = 0; i < max_iter && !converged; i++)
            {   
                max_change = 0;
                sycl::group_barrier(cta);
                coord2d neighbour_sum   = {0,0};
                for (uint8_t j = 0; j < 3; j++) neighbour_sum += xys_smem[ns[j]];

                // Calculate the new position of the point
                if(!fixed) newxys_smem[a] = xys_smem[a]*real_t(0.15) + (neighbour_sum/real_t(3.))*real_t(0.85);
                real_t neighbour_dist = 0;

                // Calculate the distance between neighbours
                for (uint8_t j = 0; j < 3; j++) neighbour_dist += norm(xys_smem[a] - xys_smem[d_get(ns,j)])/real_t(3);
                
                sycl::group_barrier(cta);
                real_t relative_change = 0;

                // Calculate the relative change
                if (neighbour_dist > 0 && !fixed){ 
                    relative_change = norm(xys_smem[a] - newxys_smem[a])/neighbour_dist;
                }

                // Reduce the relative change to find the maximum change
                real_t iteration_max = sycl::reduce_over_group(cta, relative_change, sycl::maximum<real_t>());
                if (iteration_max > max_change) max_change = iteration_max;

                converged = max_change <= 10*std::numeric_limits<real_t>::epsilon();

                // Update the position of the point
                xys_smem[a] = newxys_smem[a];
            }
            sycl::group_barrier(cta);
            xys_acc[a]  =  xys_smem[a];
            }
        });
    });
    return SyclEvent(std::move(tutte_done));
}

/* template <typename T, typename K>
void TutteFunctor<T,K>::operator() (SyclQueue& Q, Fullerene<T,K> fullerene, LaunchPolicy policy){
    execute_with_policy(Q, (Policy)policy, [&]{
        tutte_isomer_impl(  Q, 
                            fullerene.d_.A_cubic_, 
                            fullerene.d_.X_cubic_.template as_span<std::array<T,2>>(), 
                            newxys_[{Q,0}].template as_span<std::array<T,2>>(),
                            fixed_[{Q,0}],
                            max_change_[{Q,0}],
                            fullerene.N_);
        }
    );
    fullerene.m_.flags_ |= StatusFlag::CONVERGED_2D;
}

template <typename T, typename K>
void TutteFunctor<T,K>::operator() (SyclQueue& Q, FullereneBatchView<T,K> batch, LaunchPolicy policy){
    execute_with_policy(Q, (Policy)policy, [&]{
        tutte_batch_impl(Q, batch);
        }
    );
} */

template <typename T, typename K>
SyclEvent TutteFunctor<T,K>::compute(SyclQueue& Q, Fullerene<T,K> fullerene, 
                                    Span<std::array<T,2>> newxys,
                                    Span<bool> fixed,
                                    Span<T> max_change){
    if (! (fullerene.m_.flags_.get() & (int)StatusFlag::FULLERENEGRAPH_PREPARED)) return SyclEvent();
    auto ret_val = tutte_isomer_impl(  Q, 
                        fullerene.d_.A_cubic_, 
                        fullerene.d_.X_cubic_.template as_span<std::array<T,2>>(), 
                        newxys,
                        fixed,
                        max_change,
                        fullerene.N_);
    fullerene.m_.flags_ |= StatusFlag::CONVERGED_2D;
    return ret_val;
}

template <typename T, typename K>
SyclEvent TutteFunctor<T,K>::compute(SyclQueue& Q, FullereneBatchView<T,K> batch){
    return tutte_batch_impl(Q, batch);
}

template struct TutteFunctor<float,uint16_t>;
template struct TutteFunctor<float,uint32_t>;
//template struct TutteFunctor<double,uint16_t>;
//template struct TutteFunctor<double,uint32_t>;