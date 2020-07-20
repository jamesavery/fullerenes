import numpy as np
NA = np.newaxis

def read_cube(path):
    with open(path,'r') as f:
        title        = f.readline();
        density_type = f.readline();
        N, x0, y0, z0, Nval = f.readline().split();
        Nx, dx1, dx2, dx3   = f.readline().split();
        Ny, dy1, dy2, dy3   = f.readline().split();
        Nz, dz1, dz2, dz3   = f.readline().split();        

        N, Nval, Nx, Ny, Nz  = np.array([N,Nval,Nx,Ny,Nz],dtype=int);
        X0 = np.array([x0,y0,z0],dtype=float);
        dX = np.array([[dx1,dx2,dx3],[dy1,dy2,dy3],[dz1,dz2,dz3]],dtype=float);

        assert(Nval == 1);      # TODO: Deal with multivalued cubes

        Z = np.empty(N,dtype=int);
        Q = np.empty(N,dtype=float);
        X = np.empty((N,3),dtype=float);
        
        for i in range(N):
            tokens = f.readline().split();

            assert(len(tokens) == 5);             # Check that line is of the correct format
            Z[i] = int(tokens[0]);                # Atomic number
            Q[i] = float(tokens[1]);              # Charge
            X[i] = np.array(tokens[2:],dtype=float); # Atomic position

        data_string = f.read();
        data = np.array(data_string.split(), dtype=float).reshape(Nx,Ny,Nz);	# Data is given in z,y,x-order
        data = data.transpose(2,1,0);                                           # Reverse to get x,y,z-order

        return {
            'N' :    N,
            'Nval':  Nval,
            'Ngrid': np.array([Nx,Ny,Nz],dtype=int),
            'X0':    X0,
            'dX':    dX,             
            'Z':     Z,
            'Q':     Q,
            'X':     X,
	    'data': data
        };
            


def cube_to_hdf5(C,output_path):
    import h5py as h5

    f = h5.File(output_path,'w');

    f.create_dataset("Ngrid",data=C['Ngrid']);
    f.create_dataset("X0",data=C['X0']);
    f.create_dataset("dX",data=C['dX']);
    f.create_dataset("Z",data=C['Z']);
    f.create_dataset("Q",data=C['Q']);
    f.create_dataset("X",data=C['X']);
    f.create_dataset("data",data=C['data']);            
    
    f.close()

def log_atoms(filename):
    with open(filename,'r') as logfile:
        lines = logfile.readlines();
        
        i,N = 0,0;
        while not ("Standard orientation:" in lines[i]):
            i+=1;
        while not ("----" in lines[i+5+N]):
            N+=1;
        atom_lines = lines[i+5:i+5+N];

        Z = np.empty(N,     dtype=int);
        Q = np.empty(N,     dtype=float);
        X = np.empty((N,3), dtype=float);

        for i in range(N):
            id, sZ, sQ, x, y, z = atom_lines[i].split();
            assert(i == int(id)-1);
            
            Z[i] = int(sZ)
            Q[i] = float(sQ)+Z[i]
            X[i] = np.array([x,y,z],dtype=float)
        
        return Z,Q,X


        
        

def world_to_grid_coords(C,X):
    AAtoB = 1.88973
    X0, dX = C['X0'], np.diag(C['dX']);

    return (X-X0)/dX[NA,:]




