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

def read_log_geometry(txt):
    pattern = re.compile("Normal termination of Gaussian")
    if bool(re.search(pattern,txt)) == False: # Abnormal termination
        atoms, all_geometries = read_log_geometries(txt)
        if len(all_geometries)>0:
            return atoms, all_geometries[-1]
        else: 
            return None, None
    
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
        
        CAT = np.array([l.split()[:3] for l in lines],dtype=int)
        XYZ = np.array([l.split()[3:] for l in lines],dtype=float)
        
        geometries_array += [XYZ]

    if(len(geometries)>0):
        atom_numbers = CAT[:,1]
        atom_names   = elements[atomic_numbers]
        return atom_names, np.array(geometries_array)
    else:
        print("No geometries found")
        sys.exit(-1)


# Read and parse subset of Gaussian16 input files
def parse_input(input_lines):
    meta    = {}
    offsets = {}
    atom_names  = []
    atom_coords = []    
    frozen  = []
    job,comment = "",""

    i  = 0
    ai = 0
    input_section = 0
    for i in range(len(input_lines)):
        l = input_lines[i]
        if(input_section == 0): # meta section
            if(l[0] == '%'):
                key,val = l[1:].split('=',1)
                meta[key.strip()] = val.strip()
            elif(l[0] == '#'):
                job = l[1:]
            else:
                input_section = 1
                offsets["comment"] = i
                continue

        if(input_section == 1): # Comment segemtn
            comment += l
            if(l.strip() == ""):
                input_section = 2
                offsets["charge"] = i
                continue

        if(input_section == 2): # Charge
            charge = l
            input_section = 3
            offsets["molecule"] = i
            continue

        if(input_section == 3): # Molecular geometry
            ai += 1             # Gaussian numbers atoms 1,2,...
            split_line = l.split()

            if(len(split_line)==0):
                offsets["post_molecule"] = i
                input_section = 4
                continue
            elif(len(split_line)==4):
                a,x,y,z = split_line
            elif(len(split_line)==5):
                a,_,x,y,z = split_line
                frozen += [ai]
            else:
                print(f"Wrong length {len(split_line)} of line {split_line}")

            atom_names  += [a]
            atom_coords += [f"{x} {y} {z}"]

    return meta, job, charge, atom_names, atom_coords, frozen



