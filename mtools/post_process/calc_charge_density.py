import numpy as np
import matplotlib.pyplot as plt
import mdtraj as md
import seaborn as sns

gro_file = '../../data/RAH/neutral/equil.gro'
trj_file = '../../data/RAH/neutral/equil.xtc'
bin_width = 0.01

emim, NULL = read_charges('../../forcefields/emim.itp')
tfsi, NULL = read_charges('../../forcefields/Tf2N.itp')
charges = dict(emim, **tfsi)

atom_types = charges.keys()

first_frame = md.load_frame(trj_file, index=0, top=gro_file)

atom_density = dict()
charge_density = list()

plt.figure()
for atom in atom_types:
    ids = first_frame.topology.select('name {}'.format(atom))
    traj = md.load(trj_file, top=gro_file, atom_indices = ids)
    print('loaded atom type', atom, np.shape(traj.xyz))
    atom_density[atom] = traj.xyz[:,:,1]

    z = np.histogram(atom_density[atom].reshape((-1, 1))-0.682, bins=np.linspace(0, 6.7, num=6.7/bin_width + 1))
    plt.plot(z[1][:-1]+np.mean(z[1][:2]), z[0]*float(charges[atom])/(4.7205**2*bin_width*len(traj)))
    charge_density.append(z[0]*float(charges[atom])/(4.7205**2*bin_width*len(traj)))
plt.savefig('tmp.pdf')

plt.ylabel('Ion Number Densities (molecule/nm^3)')
plt.xlabel('Z position')
plt.legend(atom_types)

plt.figure()
tmp = np.array(charge_density)
plt.plot(z[1][:-1]+np.mean(z[1][:2]), np.sum(tmp, axis=0))
np.savetxt('data/charge-density.dat', np.vstack([z[1][:-1]+np.mean(z[1][:2]), np.sum(tmp, axis=0)]))
plt.ylabel('Charge Density (e/m^3)')
plt.xlabel('Z position')
plt.savefig('img/charge-density.pdf', bbox_inches='tight')
