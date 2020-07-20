import numpy as np
import numpy.linalg as la
NA = np.newaxis

def lstsq_plane(X):
    cm  = np.sum(X,axis=0)/len(X)
    Xcm = X-cm
    A = Xcm.T @ Xcm
    lams, U = la.eig(A)

    i = np.argmin(np.abs(lams))
    normal = U[:,i] / la.norm(U[:,i])
    a      = Xcm[0] - np.dot(Xcm[0],normal)*normal
    a     /= la.norm(a)
    b      = np.cross(normal,a)
    b     /= la.norm(b)         # For good measure, but not necessary as dot(n,a) = 0
    
    return cm, np.array([a,b,U[:,i]]), lams[i]

def face_planes(X, faces):
    N = len(faces)
    cms    = np.empty((N,3),  dtype=float)
    frames = np.empty((N,3,3),dtype=float)
    lams   = np.empty(N,      dtype=float)
    
    for i in range(N):
        cm, f, l  = lstsq_plane(X[faces[i]]);
        cms[i]    = cm
        frames[i] = f
        lams[i]   = l

    return cms, frames, lams

def plane_coordinates(X, faces, frames=None):
    N = len(faces)
    
    if frames is None:
        cms, frames, _ = face_planes(X,faces)
    else:
        cms = np.sum(X[faces],axis=0)/ N
        
    Xcm = X[faces] - cms[:,NA,:]

    coords = np.empty(faces.shape + (3,), dtype=float)

    for i in range(N):
        coords[i] = Xcm[i] @ frames[i].T

    return coords


def faces_grid(X, faces, nx=100, ny=100):
    # Face-plane coordinates of face vertices and face midpoints
    coords = plane_coordinates(X,faces)[...,:2]
    midpts = np.sum(coords,axis=-2)/coords.shape[-2]

    # Calculate bounding box large enough to fit any face
    box_x0, box_x1 = coords[...,0].min(), coords[...,0].max()
    box_y0, box_y1 = coords[...,1].min(), coords[...,1].max()

    # Generate rectangular grids from bounding boxes
    xs = np.linspace(box_x0, box_x1,nx)
    ys = np.linspace(box_y0, box_y1,ny)
    xy_grid = np.array(np.meshgrid(xs,ys)).transpose(1,2,0) # Shape (2,nx,ny) -> (nx,ny,2)

    return coords, xy_grid

# prop: values on vertices followed by values on midpoints
def interpolate_on_grid(xy_grid, faces, prop):
    image = np.empty(xy_grid.shape);

    Nf  = len(faces)
    deg = faces.shape[-1]
    
    for i in range(len(faces)):
        f = tris[i]             # e.g. 5x3x2
        for j in range(deg):
            xy1, xy2, xy3 = coords[i,j], coords[i,(j+1)%deg], mid[i]            
            p1,  p2,  p3  = prop[f[j]],  prop[f[(j+1)%deg]],  prop[Nf+i]

            a,b = barycentric(xy_grid, xy1, xy2, xy3)

            p             = a*p1 + b*p2 + (1-a-b)*p3;
            triangle_mask = (a>=0) & (b>=0) & ((1-a-b)>=0);
            image[triangle_mask] = p[triangle_mask]
            
    return image
    

# Input:
#  XY: (M,2) array of M 2D points to transform into barycentrics
#  xy1, xy2, xy3: 2D coordinates for the triangle vertices in a single triangle
# Output: (2,M) array of barycentric coordinates a,b (with c=1-a-b)
def barycentric(XY,xy1,xy2,xy3):
    T = np.array([xy1-xy3,xy2-xy3]).T
    r = XY-xy3
    
    return la.solve(T,r.reshape(-1,2).T)

# Input:
#  XY: (M,2) array of M 2D points to transform into barycentrics
#  XY1, XY2, XY3: (N,2) arrays of 2D points of the triangle vertices for N triangle
# Output: (M,N,2) array of barycentric coordinates a,b (with c=1-a-b)
def barycentrics(XY,XY1,XY2,XY3):
    T = np.stack([XY1-XY3,XY2-XY3], axis=2)
    r = XY-XY3
    return la.solve(T,r)


def interpolate_on_triangle(XY,tri,plane_points, prop, null_property=0):
    i1,  i2,  i3  = tri
    p1,  p2,  p3  = prop[i1],   prop[i2],   prop[i3]    
    xy1, xy2, xy3 = plane_points[i1], plane_points[i2], plane_points[i3]

    A, B = barycentric(XY, xy1, xy2, xy3);
    
    P = A*p1 + B*p2 + (1-A-B)*p3;
    P[(A<0) | (B<0) | (1-A-B<0)] = null_property;

    return P;

