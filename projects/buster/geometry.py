import numpy as np
import numpy.linalg as la
NA = np.newaxis

def split_norm(X, axis=-1):
    R = np.sqrt(np.sum(X*X,axis=axis));
    D = X/R[...,NA]
    return R,D

# X: n x 3, neighbours: n x d
# X[neighbours]: (n x d) x 3
def edge_displacements(X,neighbours):
    Xab = X[neighbours]-X[:,NA]       # Displacement vectors Xv-Xu (n x d x 3)
    return split_norm(Xab)


# The three angles around node a
# b   c
#  \_/
#  (a)
#   |
#   d
def a_corner_cos(Dab): #Dab: n x d x 3 are the displacement unit vectors. 
    Dab, Dac = Dab, np.roll(Dab,shift=-1,axis=1); #a = np.array(1,2,3), np.roll(a,shift=1) => (3,1,2)
    return np.sum(Dab*Dac,axis=2)                 # All the angle cosines around a, n x d


def a_corner_cos_gradient(Rab,Dab): 
    Rab = Rab[:,:,NA];
    Dab = Dab;
    
    Dac = np.roll(Dab,shift=-1,axis=1);
    Rac = np.roll(Rab,shift=-1,axis=1);

    cos_angle = np.sum(Dab*Dac,axis=2)[...,NA];
    
    return (Dab*cos_angle - Dac)/Rab + (Dac*cos_angle - Dab)/Rac;


# The three angles around node a's neighbours, either on left or right face
#b-  b+  c- c+
# \ /    \ /
#  b) F0 c)
#   \   /
# F2  a  F1
#     |
#    (d
#    / \
#   d+  d-
def bcd_corner_cos(X,neighbours,next_on_face): 
    bcd      = neighbours
    bcdplus  = next_on_face     # next_on_face(a,b)=b+ | next_on_face(a,c)=c+ | next_on_face(a,d)=d+

    Rba,  Dba  = split_norm(X[:,NA]    - X[bcd])
    Rbbp, Dbbp = split_norm(X[bcdplus] - X[bcd])

    return np.sum(Dba*Dbbp,axis=2)

def bcd_corner_cos_gradient(X, neighbours, next_on_face):
    bcd      = neighbours
    bcdplus  = next_on_face     # next_on_face(a,b)=b+ | next_on_face(a,c)=c+ | next_on_face(a,d)=d+
    
    Rba,  Dba  = split_norm(X[:,NA] - X[bcd])
    Rbbp, Dbbp = split_norm(X[bcdplus] - X[bcd])
    cos_angle  = np.sum(Dba*Dbbp,axis=-1)

    return (Dbbp - Dba*cos_angle[...,NA])/Rba[...,NA]



def harmonic_energy(Q,Q0,from_axis=0):
    error = Q-Q0;
    return np.sum(error*error,axis=tuple(range(from_axis,len(Q.shape))))

def harmonic_energy_gradient(Q,Q0,
                             gradQ,from_axis=0):
    error = Q-Q0;
    return 2*(Q-Q0)[...,NA]*gradQ

def corner_cos_gradient(X,neighbours,next_on_face,prev_on_face):
    Rab,Dab = edge_displacements(X,neighbours)
    
    grad_a        = a_corner_cos_gradient(Rab,Dab)
    grad_bcdplus  = bcd_corner_cos(X,neighbours,next_on_face)
    grad_bcdminus = bcd_corner_cos(X,neighbours,prev_on_face)

    return grad_a + grad_bcdplus + grad_bcdminus;


def edge_energy_gradient(X,neighbours,R0): #also needs constants
    Rab,Dab = edge_displacements(X,neighbours)
    return np.sum(2*(Rab-R0)[...,NA]*Dab, axis=1) 
    
    



