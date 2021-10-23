#pragma once
#include "coord3d.cu"
#include "fullerenes/gpu/isomerspace_forcefield.hh"
#include <exception>

/** Container for internal fullerene coordinates, bond lengths, inner angles, outer angles, inner dihedral angles, outer dihedral angles.
/*  Doubles as container for harmonic equillibirium constants corresponding to, bond equillibrium length, inner equillbirium angle  ... etc.
/*  Used diagnostically to test quality of optimized structure or for comparison to other implementations fullerene optimization. **/



void IsomerspaceForcefield::InternalCoordinates::allocate(const size_t N, const size_t batch_size)
{   
    if((!allocated)){
        this->N = N; 
        this->batch_size = batch_size; 
        size_t bytes = sizeof(device_coord3d)*N*batch_size;
        cudaMalloc(&this->bonds, bytes);
        cudaMalloc(&this->angles, bytes);
        cudaMalloc(&this->outer_angles_m, bytes);
        cudaMalloc(&this->outer_angles_p, bytes);
        cudaMalloc(&this->dihedrals, bytes);
        cudaMalloc(&this->outer_dihedrals_a, bytes);
        cudaMalloc(&this->outer_dihedrals_m, bytes);
        cudaMalloc(&this->outer_dihedrals_p, bytes);
        printLastCudaError("Failed to allocate Internal Coordinates on device: ");
        this->allocated = true;
    }
    else{
        printf("\n Warning: Arrays already allocated, allocation is ignored\n");
    }
}
  

void IsomerspaceForcefield::InternalCoordinates::free()
{
    if(allocated){
        cudaFree(this->bonds);
        cudaFree(this->angles);
        cudaFree(this->outer_angles_m);
        cudaFree(this->outer_angles_p);
        cudaFree(this->dihedrals);
        cudaFree(this->outer_dihedrals_a);
        cudaFree(this->outer_dihedrals_m);
        cudaFree(this->outer_dihedrals_p);
        printLastCudaError("Failed to free Coordinates Device :");
        this->allocated = false;
    }
}


void IsomerspaceForcefield::InternalCoordinates::to_file(const IsomerspaceForcefield::InternalCoordinates& coords,size_t fullereneID, std::string ID){
    size_t N = coords.N;
    toBinary("Bonds" + ID + ".bin",coords.bonds,N,fullereneID);
    toBinary("Angles" + ID + ".bin",coords.angles,N,fullereneID);
    toBinary("Outer_Angles_m" + ID + ".bin",coords.outer_angles_m,N,fullereneID);
    toBinary("Outer_Angles_p" + ID + ".bin",coords.outer_angles_p,N,fullereneID);
    toBinary("Dihedrals" + ID + ".bin",coords.dihedrals,N,fullereneID);
    toBinary("Outer_Dihedrals_a" + ID + ".bin",coords.outer_dihedrals_a,N,fullereneID);
    toBinary("Outer_Dihedrals_m" + ID + ".bin",coords.outer_dihedrals_m,N,fullereneID);
    toBinary("Outer_Dihedrals_p" + ID + ".bin",coords.outer_dihedrals_p,N,fullereneID);

}
    



/** Container for cartesian fullerene coordinates and fullerene graph information generated by auxillary methods 
 *  Used to primarily to copy starting geometry and graph from host to device and copy the end-state back to the host.
 * **/
void IsomerspaceForcefield::DeviceGraph::allocate(const size_t N, const size_t batch_size)
{
    if (!allocated)
    {   
        this->N = N;
        this->batch_size = batch_size;
        cudaMalloc(&this->X, sizeof(device_coord3d)*N*batch_size);
        cudaMalloc(&this->neighbours, sizeof(device_node_t)*3*N*batch_size);
        cudaMalloc(&this->next_on_face, sizeof(device_node_t)*3*N*batch_size);
        cudaMalloc(&this->prev_on_face, sizeof(device_node_t)*3*N*batch_size);
        cudaMalloc(&this->face_right, sizeof(uint8_t)*3*N*batch_size);
        printLastCudaError("Failed to allocated device graph :");
        this->allocated = true;
    }
}

void IsomerspaceForcefield::DeviceGraph::allocate_host(const size_t N, const size_t batch_size){
    if (!allocated)
    {   
        this->N = N;
        this->batch_size = batch_size;
        X = new device_real_t[3*N*batch_size];
        neighbours = new device_node_t[3*N*batch_size];
        next_on_face = new device_node_t[3*N*batch_size];
        prev_on_face = new device_node_t[3*N*batch_size];
        face_right = new uint8_t[3*N*batch_size];
        this->allocated = true;
    }
}

void IsomerspaceForcefield::DeviceGraph::free()
{   
    if (allocated)
    {
        cudaFree(this->X);
        cudaFree(this->neighbours);
        cudaFree(this->next_on_face);
        cudaFree(this->prev_on_face);
        cudaFree(this->face_right);
        printLastCudaError("Failed to free device graph :");
        this->allocated = false;
    }
}

void IsomerspaceForcefield::DeviceGraph::free_host()
{   
    if (allocated)
    {
        delete X;
        delete neighbours;
        delete next_on_face;
        delete prev_on_face;
        delete face_right;
        printLastCudaError("Failed to free device graph :");
        this->allocated = false;
    }
}




void IsomerspaceForcefield::DeviceGraph::copy_to_gpu(const IsomerspaceForcefield::DeviceGraph& G){
    cudaMemcpy(X, G.X, sizeof(device_coord3d)*N*G.batch_size , cudaMemcpyHostToDevice);
    cudaMemcpy(neighbours, G.neighbours, sizeof(device_node_t)*3*N*G.batch_size, cudaMemcpyHostToDevice);
    cudaMemcpy(next_on_face, G.next_on_face, sizeof(device_node_t)*3*N*G.batch_size, cudaMemcpyHostToDevice);
    cudaMemcpy(prev_on_face, G.prev_on_face, sizeof(device_node_t)*3*N*G.batch_size, cudaMemcpyHostToDevice);
    cudaMemcpy(face_right, G.face_right, sizeof(uint8_t)*3*N*G.batch_size, cudaMemcpyHostToDevice);
    printLastCudaError("Failed to copy DeviceGraph :");
}



