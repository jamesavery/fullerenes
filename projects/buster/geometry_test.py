import unittest 
import numpy as np
import numpy.linalg as la
from numpy import array, pi, cos, sin, sqrt, linspace
from geometry import *

NA = np.newaxis;

############### CIRCLE-TRIANGLE TESTS (CT) ###############
#
# Generate a geometry with easily understood solution.
#            c = (cos(alpha),sin(alpha))
#           /  
#  (0,0) = a----b = (1,0)
# 
# for alpha \in [-2\pi;2\pi[. We look at the effect of moving c around.
#
alphas = linspace(-2*pi+0.0001,2*pi, 1000, endpoint=False)
dalpha = alphas[1]-alphas[0]

CT_neighbours    = array([[2,1],[0,2],[1,0]])
CT_geometries    = array([[[0,0,0],[1,0,0],[cos(a),sin(a),0]] for a in alphas])
CT_geometries_da = array([[[0,0,0],[0,0,0],[-sin(a),cos(a),0]] for a in alphas])
CT_next_on_face  = array([[1,2],[2,0],[0,1]])  # Because it's just one triangle
CT_prev_on_face  = CT_next_on_face             # Because it's just one triangle


################# LENGTHS AND DERIVATIVES ##############
# |ca| = 1,
# |cb|^2 = sin(alpha)^2 + (1-cos(alpha))^2 = sin(alpha)^2 + 1 + cos(alpha)^2 - 2*cos(alpha)
#        = 2-2*cos(alpha)
# with derivatives
# d/dalpha |ca| = 0, d/dalpha |cb| = sin(alpha)/sqrt(2-2*cos(alpha))
exact_lengths    = array([sqrt(2-2*cos(alphas)), np.ones(alphas.shape)]).T
exact_lengths_da = array([sin(alphas)/sqrt(2-2*cos(alphas)), np.zeros(alphas.shape)]).T

# BRUTE-FORCE NUMERICAL SOLUTION
num_lengths    = array([la.norm(CT_geometries[:,2]-CT_geometries[:,i],axis=-1) for i in CT_neighbours[2]]).T
num_lengths_da = np.gradient(num_lengths[:,0],dalpha)

# OUR SOLUTION
lengths    = array([edge_displacements(g,CT_neighbours)[0] for g in CT_geometries])
directions = array([edge_displacements(g,CT_neighbours)[1] for g in CT_geometries])

length_gradients = -directions; # n x d x 3
lengths_da       = np.sum(length_gradients*CT_geometries_da[:,:,NA,:],axis=-1)


################# LENGTH-ENERGY AND DERIVATIVES ##############
R0 = 0.5
energies          = (1/2)*harmonic_energy(lengths,R0,from_axis=1) # 1/2 for double counting edges -> dedges
energies_gradient = array([edge_energy_gradient(g,CT_neighbours,R0) for g in CT_geometries])
energies_da       = np.sum(energies_gradient*CT_geometries_da,axis=-1)

exact_energies    = (R0-sqrt(2-2*cos(alphas)))**2 + (R0-1)**2 + (R0-1)**2
exact_energies_da = 2*(sqrt(2-2*cos(alphas))-R0)*sin(alphas)/sqrt(2-2*cos(alphas))
num_energies_da   = np.gradient(energies,dalpha)


################# ANGLES AND DERIVATIVES ##############
ALPHA0 = pi/3

cos_a        = array([a_corner_cos(Dab) for Dab in directions])
cos_bcdplus  = array([bcd_corner_cos(g,CT_neighbours,CT_next_on_face) for g in CT_geometries])
cos_bcdminus = array([bcd_corner_cos(g,CT_neighbours,CT_prev_on_face) for g in CT_geometries])

ang = np.mod(alphas,2*pi)
exact_angles = array([ang,pi/2-ang/2, pi/2-ang/2]).T
exact_cos_a  = cos(exact_angles)
exact_sin_a  = sin(exact_angles)

exact_cos_bcdplus  = np.roll(exact_cos_a,shift=1,axis=1) # Because it's just a triangle
exact_cos_bcdminus = np.roll(exact_cos_a,shift=-1,axis=1)# Because it's just a triangle

cos_a_grad   = array([a_corner_cos_gradient(lengths[i],directions[i]) for i in range(len(alphas))]) 
cos_bcd_grad = array([bcd_corner_cos_gradient(g,CT_neighbours,CT_next_on_face) for g in CT_geometries])

cos_energies = harmonic_energy(cos_a,ALPHA0,from_axis=1) # Because it's just one triangle, else needs bcdminus too 
cos_energies_grad = harmonic_energy_gradient(cos_a,ALPHA0,
                                             cos_a_grad+cos_bcd_grad,from_axis=1)


cos_a_da   = np.sum(cos_a_grad*CT_geometries_da[...,NA,:],axis=(-1))
cos_bcd_da = np.sum(cos_bcd_grad*CT_geometries_da[...,NA,:],axis=(-1))

cos_energies_da =  np.sum(cos_energies_grad*CT_geometries_da[...,NA,:],axis=(-1))
