#include "coord3d.cu"

template <typename T>
void pad_and_append(T* memory, const T* fullerene, const size_t N){
    size_t padded_molecule_size = (32 - (N % 32)) + N;
    if (sizeof(T) != sizeof(coord3d))
    {
        N *= 3;
        padded_molecule_size *= 3;
    }
    
    for (size_t i = 0; i < N; i++)
    {
        memory[i] = fullerene[i];
    }
    
    for (size_t i = N; i < padded_molecule_size; i++)
    {
        memory[i] = (T)0;
    }
}

template <typename T>
T* synthetic_array(const size_t N, const size_t num_molecules, const T* fullerene){
    size_t padded_molecule_size = (32 - (N % 32)) + N;
    if (sizeof(T) != sizeof(coord3d))
    {
        padded_molecule_size *= 3;
    }
    T* storage_array = new T[padded_molecule_size*num_molecules];
    for (size_t i = 0; i < num_molecules; i++)
    {
        pad_and_append(&storage_array[padded_molecule_size*i],fullerene,N);
    }
    return storage_array;
}

