#!/usr/bin/python3
from gaussian_parse import *
from vedo import *
import jax.numpy as np
NA = np.newaxis
Angstrom = 1.88973

input_name, log_name = sys.argv[1:]
with open(input_name,'r') as input_file:
    meta, job, charge, atom_names, atom_coords, frozen = parse_input(input_file.readlines())
        
    with open(log_name,'r') as log_file:        
        atom_numbers, atom_names, Xss = read_log_geometries(log_file.read())

Xs = Xss[-1]*Angstrom
#shape = np.array([256,256,256])
shape = np.array([100,100,100])
(Nx,Ny,Nz) = shape

Xmin = Xs.min()
Xmax = Xs.max()

xs = np.linspace(Xmin-2,Xmax+2,Nx)
ys = np.linspace(Xmin-2,Xmax+2,Ny)
zs = np.linspace(Xmin-2,Xmax+2,Nz)


print("Xmins",Xmin)
print("Xmaxs",Xmax)

def pow(x,n): return x**n

def orb1s(Z,r):
    return 2*pow(Z,3/2)*np.exp(-Z*r)

def orb2s(Z,r):
    return 2*pow(Z/2,3/2)*(1-(Z/2)*r)*np.exp(-Z/2*r)

def orb2p(Z,r):
    return 1/np.sqrt(3)*pow(Z/2,3/2)*(Z*r)*np.exp(-Z/2*r)    
    
def orb3s(Z,r):
    return 2*pow(Z/3,3/2)*(1-(2/3)*Z*r +(2/27)*(Z*r)*(Z*r) )*np.exp(-Z/3*r)

def orb3p(Z,r):
    return 4/3*np.sqrt(2)*pow(Z/3,3/2)*(Z*r)*(1-(Z/6)*r)*np.exp(-Z/3*r)


def approximate_density(molecule, show_core=False, only_carbons=False):
    (atom_numbers,Xs) = molecule

    density = np.zeros((Nx,Ny,Nz),dtype=np.float32)
    for i in range(len(Xs)):
        Z  = atom_numbers[i]

        (x,y,z) = Xs[i]
        rx = (xs-x)**2
        ry = (ys-y)**2
        rz = (zs-z)**2

        r = np.sqrt(rx[:,NA,NA] + ry[NA,:,NA] + rz[NA,NA,:] )
       
        if(Z==6):               # Carbon
            if(show_core): density += 2*np.abs(orb1s(6,r))
            density += 2*np.abs(orb2s(4,r))
            density += 2*np.abs(orb2p(4,r))        

        if(only_carbons): continue
            
        if(Z==1):               # Hydrogen
            density += np.abs(orb1s(1,r))            
            
        if(Z==9):               # Fluorine
            if(show_core):
                density += 2*np.abs(orb1s(9,r))
            density += 2*np.abs(orb2s(7,r))
            density += 5*np.abs(orb2p(7,r))

        if(Z==17):              # Chlorine
            if(show_core):
                density += 2*np.abs(orb1s(17,r))
                density += 2*np.abs(orb2s(15,r))
                density += 6*np.abs(orb2p(15,r))
            density += 2*np.abs(orb3s(10,r))
            density += 5*np.abs(orb3p(10,r))                        

    return density

density = approximate_density([atom_numbers,Xs],show_core=False, only_carbons=False)
vol = Volume(density)

#show([vol.isosurface(threshold=4)])
vol.color(["green", "pink", "blue"])
vol.alpha([0, 0.15, 0.15, 0.15, 0.2, 0.25, 0.5,0.7])
show([vol])
