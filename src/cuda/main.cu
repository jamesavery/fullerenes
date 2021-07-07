#include "kernel_shared.cu"
#include "C60ih.cu"

int main(){

    /**
     * Args: 
     * 
     * real_t X*                    -Pointer to starting geometry 
     * node_t cubic_neighbours      -Pointer to neighbour array
     * node_t next_on_face          -Pointer to array with cyclical neighbour information
     * node_t prev_on_face          -Pointer to array with cyclical (counter clockwise) neighbour information
     * uint8_t face_right           -Pointer to array with face information to the right of the vertex A->B
     * 
     * These arrays are assumed to be contigous containers with information for all fullerenes in the batch. I.e. 
     * real_t X[N*3*M] = {{x0_0,y0_0,z0_0, x0_1,y0_1,z0_1 , ..., x0_N,y0_N,z0_N},{x1_0,y1_0,z1_0, x1_1,y1_1,z1_1 , ..., x1_N,y1_N,z1_N}, ..., {x2_0,y2_0,z2_0, x2_1,y2_1,z2_1 , ..., x2_N,y2_N,z2_N}}
     * Dimensions: N x 3 = Fullerene size * 3 coordinates, M = Batch size. 
     * 
     * 
     * Batch size depends on the fullerene size N in the following way:
     * batchSize = ( maxThreadsPerBlock // N ) * numOfMultiProcessors
     * However the API call cudaOccupancyMaxActiveBlocksPerMultiprocessor() should be used.
     * 
    **/
    size_t N = 60;
    size_t batch_size = computeBatchSize(N);
    printf("Solving %d fullerenes of size: %d \n", batch_size, N);

    /** Generates a synthetic load from a single set of fullerene pointers **/
    real_t* synth_X = reinterpret_cast<real_t*>(synthetic_array<real_t>(N, batch_size, &X[0]));
    node_t* synth_cubic_neighbours = reinterpret_cast<node_t*>(synthetic_array<node_t>(N, batch_size, &cubic_neighbours[0]));
    node_t* synth_next_on_face = reinterpret_cast<node_t*>(synthetic_array<node_t>(N, batch_size, &next_on_face[0]));
    node_t* synth_prev_on_face = reinterpret_cast<node_t*>(synthetic_array<node_t>(N, batch_size, &prev_on_face[0]));
    uint8_t* synth_face_right = reinterpret_cast<uint8_t*>(synthetic_array<uint8_t>(N, batch_size, &face_right[0]));

    callKernelSingleBlockFullerenes(synth_X,synth_cubic_neighbours,synth_next_on_face,synth_prev_on_face,synth_face_right,N,batch_size);
}
