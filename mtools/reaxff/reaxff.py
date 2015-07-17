import numpy as np

def read_reaxc_bonds(bonds_file = 'bonds.reaxc', cutoff=0.3):
    """Takes bonds file, minimum 'bond order' to consider

    Returns
    -------
    bonds : np.ndarray, shape=(n, 3), dtype=(int, int, float)
        Columns: atom_1, atom_2, bond_order
    """
    bonds = np.ndarray(shape=(0,3))

    with open(bonds_file,'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                fields = line.split()
                nb = int(fields[2])
                for b in range(nb):
                    if float(fields[4+nb+b]) > cutoff:
                        atom_1 = int(fields[0])
                        if atom_1 in bonds[:,1]:
                            continue
                        else:
                            atom_2 = int(fields[3+b])
                            bond_order = float(fields[4+nb+b])
                            bonds = np.vstack([bonds, [atom_1, atom_2, bond_order]])
    return bonds


def calculate_coordination_number(bonds, cutoff=0.3):
    """Takes array of bond indices, minimum 'bond order' to consider
    
    Returns
    -------
    coordination_number : np.array, shape=(n, 1), dtype=int
    """
    coordination_number = list()
    for atom in np.unique(bonds[:,:2]):
        coordination_number.append(len([row for row in bonds[np.where(bonds[:,[0,1]] == atom),:][0] if row[2] > cutoff]))
    coordination_number = np.array(coordination_number)
    return coordination_number
