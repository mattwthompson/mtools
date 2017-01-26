import mdtraj as md
import numpy as np
import pdb
import matplotlib.pyplot as plt

top_file = 'test-data/npt.pdb'
trj_file = 'test-data/npt.xtc'

selection_box = [[0, 2], [0, 2], [0, 2]]

fig, ax = plt.subplots()

#iterate over trajectory in chunks
for n, chunk in enumerate(md.iterload(trj_file, top=top_file,
    chunk=100, stride=10, skip=100)):
    nlist = list()
    print('Loaded chunk {}'.format(n))
    indices = [[at.index for at in compound.atoms] for compound in list(chunk.topology.residues)]

#for first frame in each chunk, calculate COM of each molecule
    for i, ids in enumerate(indices):
        sub_chunk = chunk.atom_slice(ids)
        sub_com = md.compute_center_of_mass(sub_chunk[0])
        if (selection_box[0][0] < sub_com[0][0] < selection_box[0][1]):
            if (selection_box[1][0] < sub_com[0][1] < selection_box[1][1]):
                if (selection_box[2][0] < sub_com[0][2] < selection_box[2][1]):
                    nlist.append(ids)
    
    flattened = [val for sublist in nlist for val in sublist]
    nlist = flattened
#for molecules with COM inside selection_box, compute MSD along chunk
    msd = []
    msd_slice = chunk.atom_slice(nlist)
    for i in range(msd_slice.n_frames):
        temp = 0
        
        #loops over only the indices of atoms satisfying criteria
        for j in [atom.index for atom in msd_slice.topology.atoms]:
            for k in range(3):
                temp += np.square(msd_slice.xyz[i,j,k] - msd_slice.xyz[0,j,k])
        temp /= msd_slice.n_atoms
        msd.append(temp)
        
    ax.plot([t - msd_slice[0].time for t in msd_slice.time], msd, label='{}'.format(n))
    x = np.array([t - msd_slice[0].time for t in msd_slice.time])
    x = x.flatten()
    y = msd
    D = np.polyfit(x[10:], y[10:], deg=1)
    print(D[1]*10**-9)
    ax.plot(x, y, label='{}'.format(n))

ax.legend(loc=0)
fig.savefig('img/msd.pdf')
