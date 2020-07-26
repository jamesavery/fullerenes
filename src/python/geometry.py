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


def faces_grid(coords, faces, nx=100, ny=100):
    # Calculate bounding box large enough to fit any face
    box_x0, box_x1 = coords[...,0].min(), coords[...,0].max()
    box_y0, box_y1 = coords[...,1].min(), coords[...,1].max()

    # Generate rectangular grids from bounding boxes
    xs = np.linspace(box_x0, box_x1,nx)
    ys = np.linspace(box_y0, box_y1,ny)
    xy_grid = np.array(np.meshgrid(xs,ys)).transpose(1,2,0) # Shape (2,nx,ny) -> (nx,ny,2)

    return xy_grid

# prop: values on vertices followed by values on midpoints
def interpolate_on_grid(xy_grid, faces, face_coords, vertex_prop, facecenter_prop):
    images = np.zeros((len(faces),)+xy_grid.shape[:2]+vertex_prop.shape[-1:]);
#    print("images:",images.shape)

    Nf  = len(faces)
    deg = faces.shape[-1]
    
    for i in range(len(faces)):
        f = faces[i]
        mid = np.sum(face_coords[i])/len(f)
        for j in range(deg):
            xy1, xy2, xy3 = face_coords[i,j], face_coords[i,(j+1)%deg], mid
            p1,  p2,  p3  = vertex_prop[f[j]],  vertex_prop[f[(j+1)%deg]],  facecenter_prop[i]

#            print("ps:",p1.shape,p2.shape,p3.shape)
            a,b = barycentric(xy_grid, xy1, xy2, xy3)

            # TODO: Get rid of Batman
            p             = a[:,:,NA]*p1[NA,NA,:] + b[:,:,NA]*p2[NA,NA,:] + (1-a-b)[:,:,NA]*p3[NA,NA,:];
            triangle_mask = (a>=0) & (b>=0) & ((1-a-b)>=0);
            images[i][triangle_mask] = p[triangle_mask]
            
    return images
    

# Input:
#  XY: (M,2) array of M 2D points to transform into barycentrics
#  xy1, xy2, xy3: 2D coordinates for the triangle vertices in a single triangle
# Output: (2,M) array of barycentric coordinates a,b (with c=1-a-b)
def barycentric(XY,xy1,xy2,xy3):
    T = np.array([xy1-xy3,xy2-xy3]).T
    r = XY-xy3
    shape = (2,)+XY.shape[:2]
    # print("XY:",XY.shape)
    # print("xy1:",xy1)
    # print("xy2:",xy2)
    # print("xy3:",xy3)            
    # print("T:",T.shape)
    # print("r:",r.shape)    
    
    return la.solve(T,r.reshape(-1,2).T).reshape(shape)

def sample_trilinear(image,sample_points): # Sample points must already be in ijk-space
    ix,jy,kz = sample_points[...,0], sample_points[...,1], sample_points[...,2]
    ni,nj,nk = image.shape
    
    xminus,iminus= np.modf(ix-0.5);
    xplus, iplus = np.modf(ix+0.5);
    
    yminus,jminus= np.modf(jy-0.5);
    yplus, jplus = np.modf(jy+0.5);
    
    zminus,kminus= np.modf(kz-0.5);
    zplus, kplus = np.modf(kz+0.5); 

    iminus = np.maximum(iminus,0)
    jminus = np.maximum(jminus,0)
    kminus = np.maximum(kminus,0)    
    
    nx,ny,nz = image.shape[-3:];
    Dxx, Uxx = (iminus*nx*ny).astype(np.uint64), (iplus*nx*ny).astype(np.uint64);
    xDx, xUx = (jminus*ny   ).astype(np.uint64), (jplus*ny   ).astype(np.uint64);
    xxD, xxU = (kminus      ).astype(np.uint64), (kplus      ).astype(np.uint64);

    Dcc, Ucc = (1-xminus), xplus;
    cDc, cUc = (1-yminus), yplus;
    ccD, ccU = (1-zminus), zplus;
    
    I = image.reshape((nx*ny*nz));
    
    I_samples = Dcc*cDc*ccD*I[Dxx+xDx+xxD] \
              + Dcc*cDc*ccU*I[Dxx+xDx+xxU] \
              + Dcc*cUc*ccD*I[Dxx+xUx+xxD] \
              + Dcc*cUc*ccU*I[Dxx+xUx+xxU] \
              + Ucc*cDc*ccD*I[Uxx+xDx+xxD] \
              + Ucc*cDc*ccU*I[Uxx+xDx+xxU] \
              + Ucc*cUc*ccD*I[Uxx+xUx+xxD] \
              + Ucc*cUc*ccU*I[Uxx+xUx+xxU];

    # Check that everything is sane
    assert(((ix>=0) & (jy>=0) & (kz >= 0)).all())
    assert(((ix<ni) & (jy<nj) & (kz < nk)).all())            
#    assert(((xminus>=0) & (xplus>=0) & (yminus >= 0) & (yplus >= 0) & (zminus>=0) & (zplus >= 0)).all())
    assert(((xminus<=1) & (xplus<=1) & (yminus <= 1) & (yplus <= 1) & (zminus<=1) & (zplus <= 1)).all())        
    assert(((iminus>=0) & (iplus>=0) & (jminus >= 0) & (jplus >= 0) & (kminus>=0) & (kplus >= 0)).all())   
    
    volume = Dcc*cDc*ccD   \
             + Dcc*cDc*ccU \
             + Dcc*cUc*ccD \
             + Dcc*cUc*ccU \
             + Ucc*cDc*ccD \
             + Ucc*cDc*ccU \
             + Ucc*cUc*ccD \
             + Ucc*cUc*ccU;
    
    if(not(np.abs(volume-1)<1e-10).all()):
        print("Large sampling error: ",np.max(np.abs(volume-1)));
    
    return I_samples

# image_coords: (x,y,z) = (cx,cy,cz) + (dx,dy,dz)*(i,j,k)
# image_coords -> sample_coords: (i,j,k) = ((x,y,z)-(cx,cy,cz)) / (dx,dy,dz)
def sample_nearest(image,sample_points):
    ix,jy,kz  = sample_points[...,0], sample_points[...,1], sample_points[...,2]
    i, j, k   = np.round(ix).astype(int), np.round(jy).astype(int), np.round(kz).astype(int);
    print("Ix:",np.array([ix,jy,kz]).min(),np.array([ix,jy,kz]).max())
    print("I:", np.array([ix,jy,kz]).min(),np.array([i,j,k]).max())
    I_samples = image[k,j,i];   # kji vs ijk?!?!?!?!?!

    return I_samples

# Sample a general d-dimensional volume
def sample_volume(image,sample_points):
    # sample_points: ... x d
    # Xminus: ... x d float
    # Iminus: ... x d int
    Xminus, Iminus = np.modf(sample_points-0.5) 
    Xplus,  Iplus  = np.modf(sample_points+0.5)    

    Xfraction = np.stack([1-Xminus,Xplus],axis=0)

    output_shape = sample_points.shape[:-1]
    I_resampled = np.zeros(output_shape)

    for sign in range(8):       # 2^d
        term  = np.ones(output_shape)
        index = np.zeros(output_shape,dtype=np.uint64)
        
        for i in range(3):      # d
            sign_i = ((sign>>i)&1)
            term *= Xfraction[sign_i][...,i]


    
    
    
    LD = (iminus*nx+jminus).astype(np.uint64); # x-,y-
    LU = (iplus*nx +jminus).astype(np.uint64); # x-,y+
    RD = (iminus*nx+jplus).astype(np.uint64);  # x+,y-
    RU = (iplus*nx +jplus).astype(np.uint64);  # x+,y+
    
    I = image.reshape((-1,nx*ny));
    
    I_resampled = (1-xminus)*(1-yminus)*I[:,LD] \
                  +(1-xminus)*yplus     *I[:,LU] \
                  +xplus     *yplus     *I[:,RD] \
                  +xplus*(1-yminus)     *I[:,RU];
    
    return I_resampled.reshape((-1,xs.shape[0],xs.shape[1])) 
