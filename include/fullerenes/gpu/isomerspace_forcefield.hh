#pragma once
#include <inttypes.h>

namespace IsomerspaceForcefield {
  typedef float device_real_t;
  typedef uint16_t device_node_t;


struct DevicePointers {
	device_real_t* X;
    device_real_t* X1;
    device_real_t* X2;
    device_node_t* neighbours; 
    device_node_t* next_on_face;
    device_node_t* prev_on_face;
    uint8_t* face_right;
    device_real_t* gdata;
	device_real_t* bonds;
	device_real_t* angles;
	device_real_t* dihedrals;
	device_real_t* bond_0;
	device_real_t* angle_0;
	device_real_t* dihedral_0;
	device_real_t* gradients;
	
	DevicePointers(device_real_t* X, device_real_t* X_temp, device_real_t* X2, device_node_t* neighbours, device_node_t* next_on_face, device_node_t* prev_on_face, uint8_t* face_right, device_real_t* gdata,
					device_real_t* bonds, device_real_t* angles, device_real_t* dihedrals,	device_real_t* bond_0, device_real_t* angle_0, device_real_t* dihedral_0, device_real_t* gradients) : 
        X(X), X1(X1), X2(X2), neighbours(neighbours), next_on_face(next_on_face), prev_on_face(prev_on_face), face_right(face_right), gdata(gdata), bonds(bonds), angles(angles), dihedrals(dihedrals), bond_0(bond_0), angle_0(angle_0), dihedral_0(dihedral_0), gradients(gradients){}

};

struct HostPointers
{
	device_real_t* h_X;
	device_node_t* h_cubic_neighbours;
	device_node_t* h_next_on_face;
	device_node_t* h_prev_on_face;
	uint8_t* h_face_right;
	device_real_t* bonds;
	device_real_t* angles;
	device_real_t* dihedrals;
	device_real_t* bond_0;
	device_real_t* angle_0;
	device_real_t* dihedral_0;
	device_real_t* gradients;

	HostPointers(device_real_t* h_X,
		       device_node_t* h_cubic_neighbours,
		       device_node_t* h_next_on_face,
		       device_node_t* h_prev_on_face, uint8_t* h_face_right,
			   device_real_t* bonds,
	device_real_t* angles,
	device_real_t* dihedrals,
	device_real_t* bond_0,
	device_real_t* angle_0,
	device_real_t* dihedral_0,
	device_real_t* gradients
		      ) : h_X(h_X), h_cubic_neighbours(h_cubic_neighbours), h_next_on_face(h_next_on_face), h_prev_on_face(h_prev_on_face), h_face_right(h_face_right), bonds(bonds), angles(angles), dihedrals(dihedrals), bond_0(bond_0), angle_0(angle_0), dihedral_0(dihedral_0), gradients(gradients) {}
};


void OptimizeBatch(
				DevicePointers& d_pointers,
			   HostPointers& h_pointers,
		       const size_t N,
		       const size_t batch_size,
			   const size_t MaxIter);

size_t computeBatchSize(size_t N);

void AllocateDevicePointers(DevicePointers& p, size_t N, size_t batch_size);
void FreePointers(DevicePointers& p);
void CheckBatch(DevicePointers& p, const HostPointers& h, const size_t N, const size_t batch_size);
void Gradients(DevicePointers& p, device_real_t* gradients);
void InternalCoordinates(DevicePointers& p, device_real_t* bonds,device_real_t* angles,device_real_t* dihedrals);
void HarmonicConstants(DevicePointers& p, device_real_t* bond_0, device_real_t* angle_0, device_real_t* dihedral_0);

};
