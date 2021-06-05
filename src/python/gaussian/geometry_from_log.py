#!/usr/bin/env python3
import sys, re, numpy as np

pattern_start = "orientation:\s*-*\s*Center\s*Atomic\s*Atomic\s*Coordinates\s*\(Angstroms\)\s*Number\s*Number\s*Type\s*X\s*Y\s*Z\s*-*"
pattern_end = "\s+---+\s+"

elements = np.array([
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

def read_log_geometry(txt):
    pattern = re.compile("Normal termination of Gaussian")
    if bool(re.search(pattern,txt)) == False: # Abnormal termination
        atoms, all_geometries = read_geometries(txt)
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

    atom_numbers = CAT[:,1]
    atom_names   = elements[atomic_numbers]
    return atom_names, np.array(geometries_array)


def write_geometry(X,atoms,frozen=[]):
    frozen_flat = [""] * len(X)
    for i in frozen:
        frozen_flat[i-1] = " -1 "
        
    output = ""
    for ((x,y,z), a,f) in zip(X,atoms,frozen_flat):
        output += " {:s} {:s} {:16.16f} {:16.16f} {:16.16f}\n".format(a,f,x,y,z)

    return output


default_options = {
        'task_options' : {
            'opt': None, # e.g. ReadFreeze
            'scf': 'yqc',
            'freq': None},
        'nprocs':32,
        'mem': 60,
        'model': "b3lyp/6-311++g**",
        'route_extra' : "p empiricaldispersion=gd3",
       };


def write_gaussian_input(calculation_name, X,atoms, tasks=['opt'],frozen=[],options=default_options):
    link0_section = f"%nprocshared={options['nprocs']}\r\n%mem={options['mem']}GB"
    
    route_section = f"# {options['model']} {options['route_extra']}"
    for t in ['scf','opt','freq']:
        if t in tasks:
            task_options=options['task_options'][t]

            if task_options:
                route_section += f" {t}={task_options} " 
            else:
                route_section += f" {t}"

    geometry_txt = write_geometry(X,atoms)

    return f"{link0_section}\r\n{route_section}\r\n\r\n{calculation_name}\r\n\r\n0 1\r\n{write_geometry(X,atoms,frozen)}"
    

if __name__ == "__main__":
    input_text = sys.stdin.read()
    atoms, coordinates = read_log_geometry(input_text)
    print(write_gaussian_input("Precursor folding", coordinates,atoms, ['opt']), file=sys.stdout)


