import numpy as np

def read_itp(itp_file):
    """Takes itp file

    Returns
    -------
    charges : dict,
        Keys = atomtypes
        Elements = charges
    masses : dict,
        Keys = atomtypes
        Elements = mass
    """
    charges = dict()
    masses = dict()
    read_state = False
    with open(itp_file) as fi:
        for line in fi:
        #print read_state, line
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
