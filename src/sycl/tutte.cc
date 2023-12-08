#include <CL/sycl.hpp>
#include "numeric"
#include <vector>
#include <tuple>
#include <iterator>
#include <type_traits>
#include <fullerenes/sycl-isomer-batch.hh>
#include "forcefield-includes.cc"


template<typename T, typename K>
void tutte_layout(sycl::queue& Q, IsomerBatch<T,K>& batch, const LaunchPolicy policy){
    TEMPLATE_TYPEDEFS(T,K);

    
    if(policy == LaunchPolicy::SYNC) Q.wait();
    Q.submit([&](handler& h) {
        auto N = batch.N();
        auto Nf = batch.Nf();
        auto capacity = batch.capacity();
        auto max_iter = N * 10;
        accessor xys_acc(batch.xys, h, write_only); 
        accessor cubic_neighbours(batch.cubic_neighbours, h, read_only);
        accessor statuses(batch.statuses, h, read_only);

        local_accessor<bool, 1>     smem(N, h);
        local_accessor<coord2d, 1>  xys_smem(N, h);
        local_accessor<coord2d, 1>  newxys_smem(N, h);

        h.parallel_for<class tutte>(sycl::nd_range(sycl::range(N*capacity), sycl::range(N)), [=](nd_item<1> nditem) {

            auto cta = nditem.get_group();
            auto a = nditem.get_local_linear_id();
            auto isomer_idx = nditem.get_group_linear_id();

            if (statuses[isomer_idx] != IsomerStatus::EMPTY){
            size_t offset = isomer_idx * N;

            DeviceCubicGraph FG(cubic_neighbours, offset*3); 

            node3 ns       = {cubic_neighbours[(a + offset)*3], cubic_neighbours[(a + offset)*3 + 1], cubic_neighbours[(a + offset)*3 + 2]};
            xys_smem[a]    = {0,0};

            node_t outer_face[6];
            node_t outer_face_vertex   = 0;
            uint8_t Nface = FG.get_face_oriented(0,FG[0], outer_face);    
            
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

                converged = max_change <= 100*std::numeric_limits<real_t>::epsilon();

                // Update the position of the point
                xys_smem[a] = newxys_smem[a];
            }
            sycl::group_barrier(cta);
            xys_acc[isomer_idx*N + a]  =  xys_smem[a];
            }
        });
    });
    if(policy == LaunchPolicy::SYNC) Q.wait();
}

template void tutte_layout<float,uint16_t>(sycl::queue& Q, IsomerBatch<float,uint16_t>& batch, const LaunchPolicy policy);
