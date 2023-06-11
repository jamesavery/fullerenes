#include "vector"
#include "fullerenes/gpu/kernels.hh"
#include "omp.h"
typedef std::pair<uint16_t,uint16_t> arc_t;
namespace gpu_kernels {
    namespace isomerspace_dual{
    
template <int MaxDegree>
struct GraphWrapper{
    const uint16_t* neighbours;
    const uint8_t* degrees;

    GraphWrapper(const uint16_t* neighbours, const uint8_t* degrees) : neighbours(neighbours), degrees(degrees) {}

    uint16_t dedge_ix(const uint16_t u, const uint16_t v) const{
        for (uint8_t j = 0; j < degrees[u]; j++){
            if (neighbours[u*MaxDegree + j] == v) return j;
        }
        assert(false);
    }

    uint16_t next(const uint16_t u, const uint16_t v) const{
        uint16_t j = dedge_ix(u,v);
        return neighbours[u*MaxDegree + ((j+1)%degrees[u])];
    }

    uint16_t prev(const uint16_t u, const uint16_t v) const{
        uint16_t j = dedge_ix(u,v);
        return neighbours[u*MaxDegree + ((j-1+degrees[u])%degrees[u])];
    }

    arc_t canon_arc(const uint16_t u, const uint16_t v) const{
        arc_t edge = {u,v};
        uint16_t w = next(u,v);
        if (v < u && v < w) return {v,w};
        if (w < u && w < v) return {w,u};
        return edge;
    }

};
void dualise_3(IsomerBatch& B){
    std::vector<uint16_t> triangle_numbers(6*B.n_faces, UINT16_MAX);
    std::vector<uint16_t> canon_arcs(6*B.n_faces, UINT16_MAX);
    std::vector<uint16_t> n_triangles(B.n_faces, 0); //Number of triangles that each face owns.
    std::vector<uint16_t> scan_array(B.n_faces, 0); //Scan array for prefix sum.
    std::vector<arc_t> triangle_arcs(B.n_atoms);
    #pragma omp parallel
    {
        for (size_t i = 0; i < B.isomer_capacity; ++i){
            GraphWrapper<6> G(B.dual_neighbours + i*B.n_faces*6, B.face_degrees + i*B.n_faces);
            #pragma omp for schedule(auto)
            for (size_t j = 0; j < B.n_faces; ++j){
                n_triangles[j] = 0;
                //#pragma omp simd
                for (size_t k = 0; k < G.degrees[j]; ++k){
                    arc_t carc = G.canon_arc(j, G.neighbours[j*6 + k]);
                    if(carc.first == j){
                        canon_arcs[j*6 + k] = carc.second;
                        n_triangles[j]++;
                    } else {
                        canon_arcs[j*6 + k] = UINT16_MAX;
                    }
                }
            }
            #pragma omp barrier

            uint16_t accumulator = 0;
            #pragma omp simd reduction(inscan,+:accumulator)
            for (size_t j = 0; j < B.n_faces; ++j){
                scan_array[j] = accumulator;
                #pragma omp scan exclusive(accumulator)
                accumulator += n_triangles[j];
            }
            #pragma omp barrier
            #pragma omp for schedule(auto)
            for (size_t j = 0; j < B.n_faces; ++j){
                uint8_t n = 0;
                //#pragma omp simd
                for (size_t k = 0; k < G.degrees[j]; ++k){
                    if (canon_arcs[j*6 + k] != UINT16_MAX){
                        triangle_numbers[j*6 + k] = scan_array[j] + n;
                        n++;
                        triangle_arcs[triangle_numbers[j*6 + k]] = {j,canon_arcs[j*6 + k]};
                    }
                }
            }
            #pragma omp barrier
            #pragma omp for schedule(auto)
            for (size_t j = 0; j < B.n_atoms; j++){
                uint16_t u = triangle_arcs[j].first;
                uint16_t v = triangle_arcs[j].second;
                uint16_t w = G.next(u,v);
                arc_t arc_a = G.canon_arc(v,u);
                arc_t arc_b = G.canon_arc(w,v);
                arc_t arc_c = G.canon_arc(u,w);
                B.cubic_neighbours[i*B.n_atoms*3 + j*3 + 0] = triangle_numbers[arc_a.first*6 + G.dedge_ix(arc_a.first, arc_a.second)];
                B.cubic_neighbours[i*B.n_atoms*3 + j*3 + 1] = triangle_numbers[arc_b.first*6 + G.dedge_ix(arc_b.first, arc_b.second)];
                B.cubic_neighbours[i*B.n_atoms*3 + j*3 + 2] = triangle_numbers[arc_c.first*6 + G.dedge_ix(arc_c.first, arc_c.second)];
            }
        }
    }
}
void dualise_4(IsomerBatch& B){ //CPU Parallelised version of the CUDA kernel.
    #pragma omp parallel 
    {
    std::vector<uint16_t> triangle_numbers(6*B.n_faces, UINT16_MAX);
    std::vector<arc_t> triangle_arcs(B.n_atoms);
    #pragma omp for schedule(auto)
    for (size_t i = 0; i < B.isomer_capacity; ++i){
        GraphWrapper<6> G(B.dual_neighbours + i*B.n_faces*6, B.face_degrees + i*B.n_faces);
        //#pragma omp for schedule(auto)
        uint16_t accumulator = 0;
        for (size_t j = 0; j < B.n_faces; ++j){
            //#pragma omp simd
            for (size_t k = 0; k < G.degrees[j]; ++k){
                arc_t carc = G.canon_arc(j, G.neighbours[j*6 + k]);
                if(carc.first == j){
                    triangle_numbers[j*6 + k] = accumulator;
                    triangle_arcs[accumulator] = {j,carc.second};
                    accumulator++;
                }
            }
        }
        for (int j = B.n_atoms - 1; j > -1  ; j--){
            uint16_t u = triangle_arcs[j].first;
            uint16_t v = triangle_arcs[j].second;
            uint16_t w = G.next(u,v);
            arc_t arc_a = G.canon_arc(v,u);
            arc_t arc_b = G.canon_arc(w,v);
            arc_t arc_c = G.canon_arc(u,w);
            B.cubic_neighbours[i*B.n_atoms*3 + j*3 + 0] = triangle_numbers[arc_a.first*6 + G.dedge_ix(arc_a.first, arc_a.second)];
            B.cubic_neighbours[i*B.n_atoms*3 + j*3 + 1] = triangle_numbers[arc_b.first*6 + G.dedge_ix(arc_b.first, arc_b.second)];
            B.cubic_neighbours[i*B.n_atoms*3 + j*3 + 2] = triangle_numbers[arc_c.first*6 + G.dedge_ix(arc_c.first, arc_c.second)];
        }

    }   
    }
         
} 


}  
}