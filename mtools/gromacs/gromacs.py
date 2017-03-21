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


def write_custom_96(file_name):
    """Writes .xvg file for tabulated 9-6 potential

    See http://www.gromacs.org/@api/deki/files/94/=gromacs_nb.pdf

    Arguments
    -------
    file_name : string
    """
    dr = 0.001
    rcut = 2.0
    nbins = int((rcut + 1)/dr + 1)
    with open(file_name, 'w') as fo:
        for j in range(1, nbins):
            r = dr*j
            fo.write('{0} {1} {2} {3} {4} {5} {6}\n'.format(
                r, 1/r, 1/(r*r), -1/(r**6), -6/(r**7), 1/(r**9), 9/(r**10)))


def write_custom_lj_mn(file_name='table.xvg', m=12, n=6, dr=0.001, rcut=1.5):
    """Writes .xvg file for generic Lennard Jones-like n-m potential 

    See http://www.gromacs.org/@api/deki/files/94/=gromacs_nb.pdf

    Arguments
    -------
    file_name : string, path of table file
    m : int, order of repulsive term
    n : int, order of attractive term
    dr : float, spacing of values, nm
    rcut : float, cutoff distance, nm
    """
    dr = 0.001
    rcut = 2.0
    nbins = int((rcut + 1)/dr + 1)
    with open(file_name, 'w') as fo:
        for j in range(1, nbins):
            r = dr*j
            fo.write('{0} {1} {2} {3} {4} {5} {6}\n'.format(
                r, 1/r, 1/(r*r), -1/(r**n), -n/(r**(n+1)), 1/(r**m), m/(r**(m+1))))
