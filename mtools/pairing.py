import numpy as np
import mdtraj as md
import itertools
from scipy.optimize import curve_fit


def calc_pairing(trj, cutoff, names, chunk_size=100,
                 check_reform=False, normalize=False):
    """Calculate the number of molecular pairs over a trajectory."""
    c = np.zeros(shape=(len(trj), 2))
    for i, frame in enumerate(trj):
        if i % chunk_size == 0:
            pairs = build_initial_state(frame, frame_index=0,
                                        names=names, cutoff=cutoff)
        # If no pairs were found, set number of pairs to 0
        if len(pair) == 0:
            c[i] = [frame.time[0], 0]
            continue
        for pair in pairs:
            # Check unless both pair[2] and check_reform are False
            if pair[2] or check_reform:
                pair[2] = get_paired_state(frame, pair[0], pair[1],
                                           frame_index=0, cutoff=0.8)
        c[i] = [frame.time[0], np.sum(np.array(pairs)[:, 2])]

    hbar = np.zeros(shape=(chunk_size, 2))

    for chunk in chunks(c, chunk_size):
        if len(chunk) < chunk_size:
            # Neglect data in remainder of modulus
            continue
        hbar[:, 0] = chunk[:, 0] - chunk[0, 0]
        hbar[:, 1] += chunk[:, 1]

    if normalize:
        hbar[:, 1] /= hbar[0, 1]

    return hbar


def calc_caging(trj, cutoff, names, chunk_size=100, normalize=False):
    """Calculate the number of molecular cages over a trajectory."""
    c = np.zeros(shape=(len(trj), 2))
    for i, frame in enumerate(trj):
        if i % chunk_size == 0:
            cages = build_initial_cages(frame, frame_index=0,
                                        names=names, cutoff=cutoff)
        for cage in cages:
            if not check_cage(frame, cage, names):
                cages.remove(cage)
        c[i] = [frame.time[0], len(cages)]
        print(c[i])
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
    definition. The re-forming of pairs is supported with a flag."""
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
            pair_check = get_paired_state(trj, id_i, id_j, frame_index=frame_index,
                                          cutoff=cutoff)
            if pair_check:
                current_cage.append(id_j)
        if len(current_cage) > 1:
            current_cage.append(True)
            cages.append(current_cage)

    return cages


def check_cage(trj, cage, names):
    """Check if a given cage still meets its defined criteria."""
    atom_ids = [a.index for a in trj.topology.atoms if a.name in names]

    # Check to see if any ions left the cage
    for id_j in cage[1:-2]:
        # Verify ions still exist in shell
        check = get_paired_state(trj, cage[0], id_j, frame_index=0, cutoff=0.8)
        if not check:
            return False
    # See if any new ions entered the shell
    for id_k in atom_ids:
        if id_k in cage[:-2]:
            continue
        pair_check = get_paired_state(trj, cage[0], id_k, frame_index=0, cutoff=0.8) 
        if pair_check:
            return False
    return True


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
