import cube
import C60ih as c60
import numpy as np
import matplotlib.pyplot as plt
from vedo import *

C = cube.read_cube("C60ih.cube")
Nx, Ny, Nz = C['Ngrid'];

scalar_field = np.minimum(C['data'],1)
scalar_field[Nx//2:,:,:] = 0;


vol = Volume(scalar_field)
vol.color(["green", "pink", "blue"])
vol.alpha([0, 0.1, 0.1, 0.1, 0.15])

#vol.scale(1/3)
#vol.pos(C['X0'])

lego = vol.legosurface(vmin=0.01, vmax=0.25)
lego.addScalarBar3D()

print('numpy array from Volume:', 
       vol.getPointArray().shape, 
       vol.getDataArray().shape)

Xi = cube.world_to_grid_coords(C,C['X']);
spheres = [Sphere(pos=x,r=1.5) for x in Xi]

pentagons = Mesh([Xi,c60.pentagons],c='b',alpha=1)
hexagons  = Mesh([Xi,c60.hexagons],c='o',alpha=1)
#pentagons.phong()
#hexagons.phong()

#show([vol,pentagons,hexagons])
show([[vol,pentagons,hexagons]+spheres, [lego,pentagons,hexagons]+spheres], N=2, azimuth=10)
