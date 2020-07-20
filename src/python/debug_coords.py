from cube import *
from vedo import *
from C60ih import points_opt

C = read_cube("C60ih.cube")
Z,Q,X = log_atoms("C60ih.log")

cube_spheres = [Sphere(pos=x,r=2,c='r') for x in C['X']]
log_spheres  = [Sphere(pos=x,r=2,c='b') for x in X]


show([cube_spheres])
