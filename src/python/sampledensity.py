import cube
import numpy as np
import matplotlib.pyplot as plt
from vedo  import *
from C60ih import *
from geometry import *

C = cube.read_cube("C60ih.cube")
Nx, Ny, Nz = C['Ngrid'];
Xs = C['X']; 
Xi = cube.world_to_grid_coords(C,Xs);

density = np.minimum(C['data'],0.5);
print("cube:",density.shape)
print("Xi",Xi.min(),Xi.max())
    
Pcms, Pframes, Plams = face_planes(Xi,pentagons)
Pvertex_coords = plane_coordinates(Xi,pentagons)[...,:2]
PmidXi         = np.sum(Xi[pentagons], axis=-2)/5
Pxy_grid       = faces_grid(Pvertex_coords, pentagons)

Hcms, Hframes, Hlams = face_planes(Xi,hexagons)
Hvertex_coords = plane_coordinates(Xi,hexagons)[...,:2]
HmidXi         = np.sum(Xi[hexagons], axis=-2)/6
Hxy_grid       = faces_grid(Hvertex_coords, hexagons)

PXi = interpolate_on_grid(Pxy_grid,pentagons,Pvertex_coords,Xi,PmidXi)    #      12 x 100 x 100 x 3
HXi = interpolate_on_grid(Hxy_grid,hexagons, Hvertex_coords,Xi,HmidXi)   # (Nf-12) x 100 x 100 x 3

print("PXi:",PXi.min(),PXi.max())
print("HXi:",HXi.min(),HXi.max())

Pdensities = sample_trilinear(density, PXi); #     12  x 100 x 100 
Hdensities = sample_trilinear(density, HXi); # (Nf-12) x 100 x 100

Pmesh = [Mesh([Xi[p],[range(5)]],c='b',alpha=1) for p in pentagons]
Hmesh = [Mesh([Xi[h],[range(6)]],c='o',alpha=1) for h in hexagons]

def uvscale(coords, grid):
    xy0 = grid[ 0, 0]
    xy1 = grid[-1,-1]

    result = coords-xy0
    result[...,0] /= (xy1[0]-xy0[0])
    result[...,1] /= (xy1[1]-xy0[1])

    return result

def cmap(image, vmin=0, vmax=None):
    if vmax is None:
        vmax = image.max()

    image_normalized = 255*(image-vmin)/(vmax-vmin)
    r = image_normalized
    g = image_normalized
    b = image_normalized

    return np.stack([r,g,b],axis=-1).astype(np.uint8)

from PIL import Image
for i in range(12):
    texture = cmap(Pdensities[i])
    tname = "pentagon"+str(i)+".bmp"
    Image.fromarray(texture).save(tname)
    Pmesh[i].texture(tname,tcoords=uvscale(Pvertex_coords[i],Pxy_grid))

for i in range(len(Hmesh)):
    texture = cmap(Hdensities[i])
    tname = "hexagon"+str(i)+".bmp"
    Image.fromarray(texture).save(tname)
    Hmesh[i].texture(tname,tcoords=uvscale(Hvertex_coords[i],Hxy_grid))


PXiflat = [np.array([p for p in pts if (np.linalg.norm(p)>1e-6)]) for pts in PXi.reshape(12,-1,3)]

Ppts = [shapes.Points(pts) for pts in PXiflat]

scalar_field = np.minimum(density,1)
vol = Volume(scalar_field)
vol.color(["green", "pink", "blue"])
vol.alpha([0, 0.1, 0.1, 0.1, 0.15])


#show(Pmesh+[vol])
