#!/usr/bin/env python3
import sys, re, numpy as np

elements = np.array([           # Translate from atomic number to element name
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P",
    "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru",
    "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce",
    "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf",
    "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",
    "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Uut",
    "Fl", "Uup", "Lv", "Uus", "Uuo"
],dtype=np.str)

# READ GEOMETRY AND OPTIMIZATION STEPS FROM GAUSSIAN16 LOG FILE
pattern_start = "orientation:\s*-*\s*Center\s*Atomic\s*Atomic\s*Coordinates\s*\(Angstroms\)\s*Number\s*Number\s*Type\s*X\s*Y\s*Z\s*-*"
pattern_end = "\s+---+\s+"
input_geometry_start = "Charge\s+\=\s+\d+\s+Multiplicity\s+\=\s+\d+\n"
input_geometry_end   = "\n\ \n"

def read_log_inputgeometry(txt):
    try:
        input_geometry_txt = re.split(input_geometry_end,re.split(input_geometry_start,txt)[1])[0]
    except:
        print("Log file does not contain input geometry",file=sys.stderr)
        raise

    lines = input_geometry_txt.split("\n")
    N = len(lines)    
    
    atom_names  = []
    atom_coords = ""
    frozen      = np.zeros(N,dtype=bool)

    for i in range(N):
        l = lines[i]
        print(i,l)        
        split_line = l.split()
        
        if(len(split_line)==4):
            a,x,y,z = split_line
        elif(len(split_line)==5):
            a,f,x,y,z = split_line
            frozen[i] = True
        else:
            print(f"Wrong length {len(split_line)} of line {split_line}")
            
        atom_names  += [a]
        atom_coords += f"{x} {y} {z} "
        
    return atom_names, np.fromstring(atom_coords,dtype=np.float,sep=' ').reshape(N,3), frozen

def read_log_geometry(txt):
    pattern = re.compile("Normal termination of Gaussian")
    if bool(re.search(pattern,txt)) == False: # Abnormal termination
        atom_numbers, atom_names, all_geometries = read_log_geometries(txt)
        if len(all_geometries)>0:
            return atom_names, all_geometries[-1]
        else: 
            return None, None, None
    
    tmp = re.split("0,1\\\\",txt)[-1]
    tmp = re.split("\\\\\\\\Version",tmp)[0]
    tmp1 = re.sub('[0-9]', '', tmp)
    tmp1 = re.sub('[-]', '', tmp1)
    tmp1 = re.sub('[.]', '', tmp1)
    tmp1 = re.sub('[,]', '', tmp1)
    tmp1 = re.sub('[\s]', '', tmp1)
    tmp1 = re.sub('[\\\\]', ',', tmp1)
    atoms = tmp1.split(",")
    
    tmp = re.sub('[A-Z]', '', tmp)
    tmp = re.sub("\s","",tmp)
    tmp = re.sub('[a-z]', '', tmp)
    tmp = tmp[1:]
    tmp = re.sub("\\\\",'',tmp)
    geometry = np.fromstring(tmp,sep=',').reshape(-1,3)
    return atoms, geometry

def read_log_geometries(txt, pattern_start=pattern_start, pattern_end=pattern_end):
    geometries = re.split(pattern_start,txt)[1:]
    geometries_array = []
    for geometry in geometries:
        lines = re.split(pattern_end,geometry)[0].split('\n')

        try:
            CAT = np.array([l.split()[:3] for l in lines if len(l.split())==6],dtype=int)
            XYZ = np.array([l.split()[3:] for l in lines if len(l.split())==6],dtype=float)
        except:
            print("Error parsing lines: ",lines)
        
        geometries_array += [XYZ]

    if(len(geometries)>0):
        atom_numbers = CAT[:,1]
        atom_names   = elements[atom_numbers]
        return atom_numbers, atom_names, np.array(geometries_array)
    else:
        print("No geometries found")
        sys.exit(-1)




