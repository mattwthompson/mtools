def read_itp(charge_file):
    charges = dict()
    masses = dict()
    read_state = False
    with open(charge_file) as fi:
        for line in fi:
            if line.startswith(';') or line.rstrip() == '':
                continue
            if '[' in line:
                read_state = ('atoms' in line)
                continue
            if read_state == True:
                atom_data = line.split()
                if atom_data[4] in charges.keys():
                    continue
                else:
                    charges[atom_data[4]] = atom_data[6]
                    masses[atom_data[4]] = atom_data[7]
    return charges, masses
