import numpy as np
import mdtraj as md
import parmed as pmd
import warnings
import mdtraj.core.element as Element


def calc_charge_density(coord_file, trj_file, top_file, bin_width, area,
    dim, box_range, data_path):
    """
    Return a 1-D charge density profile across a channel.

    """

    system = pmd.load_file(top_file)

    first_frame = md.load_frame(trj_file, top=coord_file, index=0)

    a = [x.type for x in system.atoms]
    indices = np.unique(a, return_index=True)[1]
    parmed_types = [a[index] for index in sorted(indices)]

    a = [x.name for x in first_frame.topology.atoms]
    indices = np.unique(a, return_index=True)[1]
    mdtraj_types = [a[index] for index in sorted(indices)]

    atom_density = dict()
    charge_density = list()

    for pt, mt in zip(parmed_types, mdtraj_types):
        charge = np.unique([x.charge for x in system.atoms if x.type == pt])

        ids = first_frame.topology.select('name {}'.format(mt))
        traj = md.load(trj_file, top=coord_file,
            atom_indices=first_frame.topology.select('name {}'.
            format(mt)))

        atom_density[mt] = traj.xyz[:, :, 2]

        x = np.histogram(atom_density[mt].reshape((-1, 1)),
            bins=np.linspace(box_range[0], box_range[1],
            num=1+round((box_range[1]-box_range[0])/bin_width)))

        charge_density.append(x[0]*charge/(area*bin_width*len(traj)))

    np.savetxt('{}/charge-density.txt'.format(data_path),
        np.transpose(np.vstack([x[1][:-1]+np.mean(x[1][:2]),
        np.sum(np.array(charge_density), axis=0)])))

def calc_number_density(coord_file, trj_file, bin_width, area, dim, box_range, data_path):
    """
    Return a 1-D number density profiles for each residue across a channel.

    """

    first_frame = md.load_frame(trj_file, top=coord_file, index=0)
    
    resnames = np.unique([x.name for x in first_frame.topology.residues])

    open('{0}/resnames.txt'.format(data_path), 'w').close()

    for resname in resnames:
        traj = md.load(trj_file, top=coord_file,
            atom_indices=first_frame.topology.select('resname {}'
                .format(resname)))

        indices = [[at.index for at in compound.atoms]
            for compound in list(traj.topology.residues)]

        if 0 in [x.mass for x in
            [atom.element for atom in traj.topology.atoms]]:
            warnings.warn("mdtraj found zero mass, setting element to hydrogen", UserWarning)
            for atom in traj.topology.atoms:
                if atom.element in [Element.virtual, Element.virtual_site]:
                    atom.element = Element.hydrogen

        com = list()
        for i, ids in enumerate(indices):
            sub_traj = traj.atom_slice(ids)
            com.append(md.compute_center_of_mass(sub_traj))
        com = np.array(com)

        x = np.histogram(com[:, 1:, dim].reshape((-1, 1)),
            bins=np.linspace(box_range[0], box_range[1],
            num=1+round((box_range[1]-box_range[0])/bin_width)))

        np.savetxt('{0}/{1}-number-density.txt'.format(data_path, resname),
            np.vstack([x[1][:-1]+np.mean(x[1][:2])-box_range[0],
            x[0]/(area*bin_width*(len(traj)-1))]).transpose())

        with open('{0}/resnames.txt'.format(data_path), "a") as myfile:
            myfile.write(resname + '\n')


def calc_msd(coord_file, trj_file, dims=[1, 1, 1]):
    """Calculate the MSD and diffuvisity of a bulk simulation

    Parameters
    ----------
    coord_file : String
        Path to topology-containing coordinate file
    trj_file : String
        Path to unwrapped trajectory file

    Returns
    -------
    D : Float
        Bulk 3-D self-diffusvity
    msd : np.ndarray
    """

    traj = md.load(trj_file, top=coord_file)

    msd = [np.sum([(row - traj.xyz[0, :, dim])**2/int(traj.n_atoms) for dim in [1, 1, 1]]) for row in traj.xyz]
    y_fit = msd
    x_fit = [val - traj.time[0] for val in traj.time]

    fit = np.polyfit(x_fit[int(np.round(len(msd)/10, 0)):],
        y_fit[int(np.round(len(msd)/10, 0)):],
        1)

    D = fit[0]/6 * 1e-6

    return D, msd, x_fit, y_fit


def calc_density(traj, units='macro'):
    vol = np.product(traj.unitcell_lengths, axis=1)
    total_mass = np.sum([x.element.mass for x in traj.topology.atoms])
    if units == 'nano':
        rho = total_mass/vol
    elif units == 'macro':
        rho = total_mass/vol * 1.66054 # Convert from amu/nm^3 to kg/m^3
    else:
        raise ValueError('Unsupported units!')
    return rho
