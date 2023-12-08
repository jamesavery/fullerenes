import numpy as np

(R55, R56, R66) = (1.479,1.458, 1.401); # Bond-lengths in Ångstrom

R0_constants = np.array([[R55, R56],
                         [R56, R66]]);

# neighbours:      n x d  - naboer orienteret med uret -> dermed også kanter
# face_size_left:  n x d
# face_size_right: n x d


R0 = R0_constants[face_size_left-5, face_size_right-5]; # n x d

NA = np.newaxis 
Rsqr = np.sum((points[:,NA] - points[neigbours])**2,axis=2) # n x d x 3 -> n x d
R    = np.sqrt(Rsqr)          # Optional: optimize Rsqr instead of R

# Computes displacement vectors along all edges
#  X(nx3):          Cartesian x,y,z coordinates of all vertices
#  neighbours(nxd): Neighbour list representation of edges; d is degree of graph (cubic => d=3)
#
# Result (u,v)-edge lengths Ruv(n x d), displacement vector directions Duv(n x d x 3) (unit vectors)
def edge_displacements(X,neighbours):
    Duv = X[neighbours]-X[:,None,:]
    Ruv = np.sqrt(np.sum(Duv*Duv,axis=2))
    return (Ruv,Duv/Ruv)

# Input:  Ruv(n x d), R0(n x d), Duv(n x d x 3)
# Result: grad(n x 3)
def edge_energy_gradient(Ruv,Duv,R0):
    return np.sum(2*(Ruv-R0)[:,:,None]*Duv, axis=1)

# Input: Duv(n x d x 3)
# Output: cos_angle(n x d)
def corner_cos_angles(Duv):
    Duv0, Duv1 = Duv, np.roll(Duv,shift=1,axis=1);
    return np.sum(Duv0*Duv1,axis=2)

def edge_energy(Ruv,R0):
    error = Ruv-R0;
    return np.sum(error*error)

def corner_cos_energy(cos_angles, cos_angles0):
    error = cos_angles-cos_angles0
    return np.sum(error*error)

# TODO: This is the gradient of cos_angles, not the energy. Fix -> corner_cos_energy_gradient
def corner_cos_gradient(Duv,Ruv,cos_angle):
    Dab = Duv;
    Rab = Ruv[:,:,None];
    Dac = np.roll(Duv,shift=1,axis=1);
    Rac = np.roll(Ruv,shift=1,axis=1)[:,:,None];    

    return cos_angle*(Dab/Rab+Dac/Rac) - Dab/Rac - Dac/Rab;






