import sys

init_com = sys.argv[1]


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

    return meta, job, atom_names, atom_coords, frozen


with open(init_com,'r') as f:
    input_lines = f.readlines()

    print(parse_input(input_lines))
