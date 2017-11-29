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

    
def get_paired_state(trj, id_i, id_j, frame_index=0, cutoff=1):
    """Check to see if a given pair is still paired."""
    dist = np.sum(np.sqrt((trj.xyz[frame_index, id_i] - trj.xyz[frame_index, id_j])**2))
    if dist < cutoff:
        paired = True
    else:
        paired = False
    return paired

def build_initial_state(trj, names, frame_index=0, cutoff=1):
    """Build initial pair list. See 10.1021/acs.jpclett.5b00003 for a definition. The re-forming of pairs is supported implicitly."""
    atom_ids = [a.index for a in trj.topology.atoms if a.name in names]
    pairs = [prod for prod in itertools.combinations(atom_ids, r=2)]
    pairs = [list([*pair, False]) for pair in pairs]

    for i, pair in enumerate(pairs):
        if pair[0] == pair[1]:
            continue
        pair[2] = get_paired_state(trj, pair[0], pair[1], frame_index=frame_index)
    pairs = [pair for pair in pairs if pair[2] == True]

    return pairs

def chunks(l, n): 
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n): 
        yield l[i:i + n]

def stretched_exp(x, a, b): 
    """Define a stretched exponential function."""
    f = np.exp(-1 * b * x ** a)
    return f

def fit(func, t, n_pairs):
    """Fit pairing data to 
    popt, pcov = curve_fit(stretched_exp, t, n_pairs)
    return popt

