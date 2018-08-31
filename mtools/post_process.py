import numpy as np
import mdtraj as md
import parmed as pmd
import warnings
import mdtraj.core.element as Element
import matplotlib.pyplot as plt
from scipy import stats

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


def calc_msd(traj, dims=[1, 1, 1], fit_with='scipy'):
    """Calculate the MSD and diffuvisity of a bulk simulation

    Parameters
    ----------
    traj : md.Trajectory
        mdtraj Trajectory

    dims : list of int, default = [1, 1, 1]
        Dimensions to consider, in order of `x, y, z`. A value of 1 will result
        in that dimension being considered.

    fit_with : str, default='scipy.stats'
        Package to use for linear fitting. Options are 'numpy' (`numpy.polyfit`) and
        'scipy' (`scipy.stats.linregress`)

    Returns
    -------
    D : Float
        Bulk 3-D self-diffusvity
    msd : np.ndarray
    """


    if dims == [1, 1, 1]:
        msd = [np.sum((traj.xyz[index, :, :] - traj.xyz[0, :, :]) ** 2) / traj.n_atoms for index in range(traj.n_frames)]
    else:
        msd = np.zeros(shape=traj.n_frames)
        for dim, check in enumerate(dims):
            if check == 1:
                msd += [np.sum((traj.xyz[index, :, dim] - traj.xyz[0, :, dim]) ** 2) / traj.n_atoms for index in range(traj.n_frames)]
            elif check != 0:
                raise ValueError('Indices of dim must be 0 or 1!')

    y = msd
    x = traj.time - traj.time[0]

    if fit_with == 'scipy':
        slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
        D = slope / (2*np.sum(dims)) * 1e-6

        x_fit = x[int(np.round(len(msd)/10, 0)):]
        y_fit = slope * x_fit + intercept

    elif fit_with == 'numpy':
        fit = np.poly1d(
                np.polyfit(x[int(np.round(len(msd)/10, 0)):],
                           y[int(np.round(len(msd)/10, 0)):],
                           1)
                )

        D = fit[0]/(2*np.sum(dims)) * 1e-6

        x_fit = x[int(np.round(len(msd)/10, 0)):]
        y_fit = fit(x_fit)

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


def slice_and_chunk(trj_file = 'traj_unwrapped.xtc', top_file='confound.pdb', chunk_len=100, skip=1, dims=[1, 1, 1],
                x_range=[0, 3], y_range=[0, 3], z_range=[0, 3],
                msd_file='msd.txt', img_file='msd.pdf'):
    for n, chunk in enumerate(md.iterload(trj_file, top=top_file, chunk=chunk_len, skip=skip)):
        print(n)

        if len(chunk) < chunk_len:
            continue

        indices = [[at.index for at in compound.atoms] for compound in list(chunk.topology.residues)]

        nlist = list()

        print('building nlist ...')

        #for i, ids in enumerate(indices):
        for ids in indices:
            sub_chunk = chunk[0].atom_slice(ids)
            sub_com = md.compute_center_of_mass(sub_chunk)[0]
            if (x_range[0] < sub_com[0] < x_range[1]):
                if (y_range[0] < sub_com[1] < y_range[1]):
                    if (z_range[0] < sub_com[2] < z_range[1]):
                        nlist.append(ids)

        nlist = [val for sublist in nlist for val in sublist]

        if len(nlist) == 0:
            continue

        slice = chunk.atom_slice(nlist)
        D, msd, x_fit, y_fit = calc_msd(slice, dims=dims)
        if n == 0:
            master_msd = msd*slice.n_atoms
            master_n_atoms = slice.n_atoms
            fix, ax = plt.subplots()
        else:
            master_msd += msd*slice.n_atoms
            master_n_atoms += slice.n_atoms
            ax.plot(slice.time-slice.time[0], msd, 'b--', alpha=0.2)
            ax.plot(x_fit, y_fit, 'k-', alpha=0.2)

    ax.plot(slice.time-slice.time[0], [x/master_n_atoms for x in master_msd], 'r-')

    x = [val - slice.time[0] for val in slice.time]
    y = [x/master_n_atoms for x in master_msd]

    fit = np.polyfit(x[int(np.round(len(master_msd)/5, 0)):],
            y[int(np.round(len(master_msd)/5, 0)):], 1)

    D = fit[0]/(2*np.sum(dims)) * 1e-6

    x_fit = x[int(np.round(len(msd)/5, 0)):]
    y_fit = [fit[0]*x + fit[1] for x in x_fit]

    ax.plot(x_fit, y_fit, 'r--')
    ax.set_title(D)
    plt.savefig(img_file)
    with open(data_file, "a") as myfile:
            myfile.write(D)


def compute_cn(r, g_r, rho_j):
    """
    Compute the coordination number of `j` around `i`.

    Note: The coordination number is not pairwise symmetric with respect to
    `i` and `j`, unlike the radial distribution funciton. For example, in a
    typical aqueous electrolyte system, an ion has more waters surrounding it
    than a water has ions surrounding it.

    Parameters
    ----------
    r : array-like
        Radii values corresponding to the centers of the RDF bins
    g_r : array-like
        Radial distribution function (RDF) values of `ij`
    rho_j : float
        Partial density of selection `j`. It is importat

    Returns
    -------
    c_n : array-like

    """

    factor = 4 * np.pi * rho_j
    c_n = [np.trapz(g_r[:i] * r[:i] ** 2, r[:i], r) for i in range(len(r))]
    c_n *= factor
    return c_n
