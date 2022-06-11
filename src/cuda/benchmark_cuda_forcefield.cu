#include "isomerspace_forcefield.cu"
#include "coord3d.cuh"
#include "fullerenes/gpu/isomerspace_forcefield.hh"

#include <unistd.h>
constexpr
int minpow2(int v)
{
  v--;
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  v++;
  return v;
}

int main(){

    /**
     * Args: 
     * 
     * device_real_t X*                    -Pointer to starting geometry 
     * device_node_t cubic_neighbours      -Pointer to neighbour array
     * device_node_t next_on_face          -Pointer to array with cyclical neighbour information
     * device_node_t prev_on_face          -Pointer to array with cyclical (counter clockwise) neighbour information
     * uint8_t face_right           -Pointer to array with face information to the right of the vertex A->B
     * 
     * These arrays are assumed to be contigous containers with information for all fullerenes in the batch. I.e. 
     * device_real_t X[N*3*M] = {{x0_0,y0_0,z0_0, x0_1,y0_1,z0_1 , ..., x0_N,y0_N,z0_N},{x1_0,y1_0,z1_0, x1_1,y1_1,z1_1 , ..., x1_N,y1_N,z1_N}, ..., {x2_0,y2_0,z2_0, x2_1,y2_1,z2_1 , ..., x2_N,y2_N,z2_N}}
     * Dimensions: N x 3 = Fullerene size * 3 coordinates, M = Batch size. 
     * 
     * 
     * Batch size depends on the fullerene size N in the following way:
     * batchSize = ( maxThreadsPerBlock // N ) * numOfMultiProcessors
     * However the API call cudaOccupancyMaxActiveBlocksPerMultiprocessor() should be used.
     * 
    **/
    const size_t N = 256;


    //size_t batch_size = IsomerspaceForcefield::computeBatchSize(N);
    IsomerspaceForcefield kernel = IsomerspaceForcefield(N);
    size_t batch_capacity = kernel.get_batch_capacity();
    printf("Solving %d fullerenes of size: %d \n", (int)batch_capacity, (int)N);

    //TODO find a better way to make a batch of fullerenes and benchmark the optimizer. Or simply use isomerspace-optimize.

    for (size_t i = 0; i < batch_capacity; i++)
    {
      //      int offset = i*N*3;
      //kernel.insert_isomer(synth_X + offset, synth_cubic_neighbours + offset, synth_next_on_face + offset, synth_prev_on_face + offset, synth_face_right + offset, i);
      
    }
    

    //IsomerspaceForcefield::IsomerspaceGraph graph = IsomerspaceForcefield::IsomerspaceGraph(synth_X,synth_cubic_neighbours, synth_next_on_face, synth_prev_on_face, synth_face_right);
    //graph.N = N; graph.batch_size = batch_size; graph.buffer_type = IsomerspaceForcefield::HOST_BUFFER;
    


    //kernel.optimize_batch(N*3);
    //kernel.check_batch(N*3);
    //kernel.optimize_batch(N*1);
    //kernel.check_batch(N*3);
    //kernel.to_file(0);
    //kernel.batch_statistics_to_file();
}
