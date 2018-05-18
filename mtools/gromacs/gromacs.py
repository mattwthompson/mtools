import os
from copy import deepcopy

import numpy as np
import mdtraj as md
from mdtraj.core.element import virtual_site
from parmed import unit as u
from parmed.parameters import ParameterSet

def make_comtrj(trj):
    """Takes a trj and returns a trj with COM positions as atoms"""

    comtop = md.Topology()
    coords = np.ndarray(shape=(trj.n_frames, trj.n_residues, 3))

    for j, res in enumerate(trj.topology.residues):
        comtop.add_atom(res.name, virtual_site, comtop.add_residue(res.name, comtop.add_chain()))
        res_frame = trj.atom_slice([at.index for at in res.atoms])
        coords[:, j, :] = md.compute_center_of_mass(res_frame)

    comtrj = md.Trajectory(xyz=coords,
                           topology=comtop,
                           time=trj.time,
                           unitcell_angles=trj.unitcell_angles,
                           unitcell_lengths=trj.unitcell_lengths)

    return comtrj

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
        for j in range(nbins):
            r = dr*j
            fo.write('{0} {1} {2} {3} {4} {5} {6}\n'.format(
                r, 1/r, 1/(r*r), -1/(r**6), -6/(r**7), 1/(r**9), 9/(r**10)))


def write_custom_lj_mn(file_name='table.xvg', m=12, n=6, dr=0.001,
    r_cut=1.5, r_correct=0.1):
    """Writes .xvg file for generic Lennard Jones-like n-m potential 

    See http://www.gromacs.org/@api/deki/files/94/=gromacs_nb.pdf

    Arguments
    -------
    file_name : string, path of table file
    m : int, order of repulsive term
    n : int, order of attractive term
    dr : float, spacing of values, nm
    r_cut : float, cutoff distance, nm
    r_correct: float, distance to switch to exponential head, nm
    """
    
    nbins = int(1 + (r_cut)/dr)
    
    with open(file_name, 'w') as fo:
        for j in range(nbins):
            r = dr*float(j)
            # Fit exponential head correction
            Bm = np.log((r_correct/(r_correct - dr))**m)/dr
            Am = r_correct**-m * np.exp(Bm*r_correct)

            if n == 0:
                Bn = 0
                An = 0
            else:
                Bn = np.log((r_correct/(r_correct - dr))**n)/dr
                An = -1*r_correct**-n * np.exp(Bn*r_correct)

            if r == 0:
                fo.write('{0} {1} {2} {3} {4} {5} {6}\n'.format(
                    r,
                    1/(r+dr),
                    1/(r+dr)**2,
                    An*np.exp(-Bn*r),
                    An*Bn*np.exp(-Bn*r),
                    Am*np.exp(-Bm*r),
                    Am*Bm*np.exp(-Bm*r)))
            elif r < r_correct:
                fo.write('{0} {1} {2} {3} {4} {5} {6}\n'.format(
                    r,
                    1/r,
                    1/(r*r),
                    An*np.exp(-Bn*r),
                    An*Bn*np.exp(-Bn*r),
                    Am*np.exp(-Bm*r),
                    Am*Bm*np.exp(-Bm*r)))
            else:
                fo.write('{0} {1} {2} {3} {4} {5} {6}\n'.format(
                    r,
                    1/r,
                    1/(r*r),
                    -1/(r**n),
                    -n/(r**(n+1)),
                    1/(r**m),
                    m/(r**(m+1))))


def unwrap_trj(filename):
    """Wrapper around ```gmx trjconv``` that unwraps GROMACS trajectories.

    See http://manual.gromacs.org/programs/gmx-trjconv.html

    Arguments
    -------
    filename : string, path of table file
    """
    os.system('gmx trjconv -f {0} -o {1}_unwrapped.{2} -pbc nojump'.format(
        filename,
        filename.rsplit('.')[:-1][0],
        filename.rsplit('.')[-1]))


def parse_nonbond_params(top_file, gmxtop):
    """Parse the `nonbond_params` directive in a GROMACS topology and store the
    information in a parmed gmxtop.

    Code is largely based off of ParmEd's reading of this directive:
    ParmEd/ParmEd/blob/master/parmed/gromacs/gromacstop.py

    """

    params = deepcopy(gmxtop.parameterset)
    atom_types = [key for key in params.atom_types.keys()]
    current_section = None
    with open(top_file, 'r') as top:
        for line in top:
            line = line.strip()
            if not line:
                continue
            if line[0] == ';':
                continue
            if line[0] == '[':
                current_section = line[1:-1].strip()
                continue
            if current_section == 'nonbond_params':
                words = line.split()
                a1, a2 = words[:2]
                if a1 in atom_types and a2 in atom_types:
                    sig, eps = (float(x) for x in words[3:5])
                    sig *= 10 # Convert to Angstroms
                    eps *= u.kilojoule.conversion_factor_to(u.kilocalorie)
                    params.nbfix_types[(a1, a2)] = (eps, sig*2**(1/6))
                    params.nbfix_types[(a2, a1)] = (eps, sig*2**(1/6))
                    params.atom_types[a1].add_nbfix(a2, sig*2**(1/6), eps)
                    params.atom_types[a2].add_nbfix(a1, sig*2**(1/6), eps)
    gmxtop.parameterset = deepcopy(params)
    return deepcopy(gmxtop)
