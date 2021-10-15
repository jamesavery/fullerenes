#include "kernel_clean.cu"
#include "coord3d.cu"
#include "C256ih.cu"


typedef uint16_t node_t;
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
    const size_t N = 256;
    const int Block_Size_Pow2 = minpow2((int)N);

    size_t batch_size = IsomerspaceForcefield::computeBatchSize<Block_Size_Pow2>(N);
    
    printf("Solving %d fullerenes of size: %d \n", (int)batch_size, (int)N);

    /** Generates a synthetic load from a single set of fullerene pointers **/
    real_t* synth_X = reinterpret_cast<real_t*>(IsomerspaceForcefield::synthetic_array<real_t>(N, batch_size, &X[0]));
    node_t* synth_cubic_neighbours = reinterpret_cast<node_t*>(IsomerspaceForcefield::synthetic_array<node_t>(N, batch_size, &cubic_neighbours[0]));
    node_t* synth_next_on_face = reinterpret_cast<node_t*>(IsomerspaceForcefield::synthetic_array<node_t>(N, batch_size, &next_on_face[0]));
    node_t* synth_prev_on_face = reinterpret_cast<node_t*>(IsomerspaceForcefield::synthetic_array<node_t>(N, batch_size, &prev_on_face[0]));
    uint8_t* synth_face_right = reinterpret_cast<uint8_t*>(IsomerspaceForcefield::synthetic_array<uint8_t>(N, batch_size, &face_right[0]));
    
    real_t* bonds = new real_t[batch_size*N*3];
    real_t* angles = new real_t[batch_size*N*3];
    real_t* dihedrals = new real_t[batch_size*N*3];
    real_t* bond_0 = new real_t[batch_size*N*3];
    real_t* angle_0 = new real_t[batch_size*N*3];
    real_t* dihedral_0 = new real_t[batch_size*N*3];
    real_t* gradients =  new real_t[batch_size*N*3];
    
    real_t* d_X; real_t* d_X_temp; real_t* d_X2; node_t* d_neighbours; node_t* d_prev_on_face; node_t* d_next_on_face; uint8_t* d_face_right; real_t* d_gdata; real_t* d_bonds; real_t* d_angles; real_t* d_dihedrals; real_t* d_angle_0; real_t* d_bond_0; real_t* d_dihedral_0; real_t* d_gradients;
    IsomerspaceForcefield::DevicePointers d_pointers = IsomerspaceForcefield::DevicePointers(d_X,d_X_temp,d_X2,d_neighbours,d_prev_on_face, d_next_on_face, d_face_right, d_gdata,d_bonds,d_angles,d_dihedrals,d_bond_0,d_angle_0,d_dihedral_0,d_gradients);
    IsomerspaceForcefield::HostPointers h_pointers = IsomerspaceForcefield::HostPointers(synth_X, synth_cubic_neighbours, synth_next_on_face, synth_prev_on_face, synth_face_right);
    IsomerspaceForcefield::AllocateDevicePointers(d_pointers, N, batch_size);
    IsomerspaceForcefield::OptimizeBatch<Block_Size_Pow2>(d_pointers,h_pointers,N,batch_size,N*5);
    IsomerspaceForcefield::CheckBatch<Block_Size_Pow2>(d_pointers, h_pointers, N, batch_size);
    IsomerspaceForcefield::InternalCoordinates<Block_Size_Pow2>(d_pointers,h_pointers,N,batch_size,bonds,angles,dihedrals);
    IsomerspaceForcefield::HarmonicConstants<Block_Size_Pow2>(d_pointers,h_pointers,N,batch_size,bond_0,angle_0,dihedral_0);
    IsomerspaceForcefield::Gradients<Block_Size_Pow2>(d_pointers,h_pointers,N,batch_size,gradients);
    IsomerspaceForcefield::FreePointers(d_pointers);

    //IsomerspaceForcefield::print_array(reinterpret_cast<IsomerspaceForcefield::coord3d*>(bonds),N,0);


}
