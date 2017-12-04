import numpy as np
import mdtraj as md
from scipy.optimize import curve_fit


def calc_pairing(trj, cutoff, names, chunk_size=100, normalize=False):
    """Calculate the number of molecular pairs over a trajectory."""
    c = np.zeros(len(trj))
    for i, frame in enumerate(trj):
        if i % chunk_size == 0:
            pairs = build_initial_state(frame, frame_index=0,
                                        names=names, cutoff=cutoff)
        # If no pairs were found, set number of pairs to 0
        if len(pair) == 0:
            c[i] = 0
            continue
        for pair in pairs:
            if pair[2]:
                pair[2] = get_paired_state(frame, pair[0], pair[1],
                                           frame_index=0, cutoff=0.8)
        c[i] = [frame.time[0], np.sum(np.array(pairs)[:, 2])]

    hbar = np.zeros(shape=(chunk_size, 2))

    for chunk in chunks(c, chunk_size):
        if len(chunk) < chunk_size:
            continue
        hbar[:, 0] = chunk[:, 0] - chunk[0, 0]
        hbar[:, 1] += chunk[:, 1]

    if normalize:
        hbar[:, 1] /= hbar[0, 1]

    return hbar


def calc_caging(trj, cutoff, names, chunk_size=100, normalize=False):
    """Calculate the number of molecular cages over a trajectory."""
    c = np.zeros(len(trj))
    for i, frame in enumerate(trj):
        if i % chunk_size == 0:
            cages = build_initial_state(frame, frame_index=0,
                                        names=names, cutoff=cutoff)
        for cage in cages:
            cage = check_cage(cage)
            if not cage[-1]:
                cages.remove(cage)
        c[i] = [frame.time[0], len(cages)]

    hbar = np.zeros(shape=(chunk_size, 2))

    for chunk in chunks(c, chunk_size):
        if len(chunk) < chunk_size:
            continue
        hbar[:, 0] = chunk[:, 0] - chunk[0, 0]
        hbar[:, 1] += chunk[:, 1]

    if normalize:
        hbar[:, 1] /= hbar[0, 1]

    return hbar


def get_paired_state(trj, id_i, id_j, frame_index=0, cutoff=1):
    """Check to see if a given pair is still paired."""
    dist = np.sum(np.sqrt((trj.xyz[frame_index, id_i] -
                           trj.xyz[frame_index, id_j]) ** 2))
    if dist < cutoff:
        paired = True
    else:
        paired = False
    return paired


def build_initial_state(trj, names, frame_index=0, cutoff=1):
    """Build initial pair list. See 10.1021/acs.jpclett.5b00003 for a
    definition. The re-forming of pairs is supported implicitly."""
    atom_ids = [a.index for a in trj.topology.atoms if a.name in names]
    pairs = [prod for prod in itertools.combinations(atom_ids, r=2)]
    pairs = [list([*pair, False]) for pair in pairs]

    for i, pair in enumerate(pairs):
        if pair[0] == pair[1]:
            continue
        pair[2] = get_paired_state(trj, pair[0], pair[1],
                                   frame_index=frame_index)
    pairs = [pair for pair in pairs if pair[2] == True]

    return pairs


def build_initial_cages(trj, names, frame_index=0, cutoff=1):
    """Build initial cage list. See 10.1021/acs.jpclett.5b00003 for a
    definition. The re-forming of cages is not permitted."""
    atom_ids = [a.index for a in trj.topology.atoms if a.name in names]

    cages = list()

    for id_i in atom_ids:
        current_cage = list()
        current_cage.append(id_i)
        for id_j in atom_ids:
            pair_check = get_paired_state(id_i, id_j, frame_index=frame_index,
                                          cutoff=cutoff)
            if pair_check:
                current_cage.append(id_j)
        current_cage.append(False)
        cages.append(current_cage)

    return cages


def check_cage(cage):
    """Check if a given cage still meets its defined criteria."""
    id_i = cage[0]
    for id_j in cage[1:-2]:
        check = get_paired_state(id_i, id_j, frame_index=0, cutoff=0.8)
        if check is False:
            cage[-1] = check
            return cage
    return cage


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


def stretched_exp(x, a, b):
    """Define a stretched exponential function."""
    f = np.exp(-1 * b * x ** a)
    return f


def fit(func, t, n_pairs):
    """Fit pairing data to a stretched exponential function"""
    popt, pcov = curve_fit(stretched_exp, t, n_pairs)
    return popt
