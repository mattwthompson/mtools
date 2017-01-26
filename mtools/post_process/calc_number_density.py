import numpy as np
import mdtraj as md
import pdb
import matplotlib.pyplot as plt

#gro_file = 'test-data/confout.gro'
#trj_file = 'test-data/traj.xtc'
#
#bin_width = 0.067/2
#area = 6.15*5.538
#
#first_frame = md.load_frame(trj_file, top=gro_file, index=0)
#
#ions = dict()
#ions['EMIM'] = first_frame.topology.select('resname EMIM')
#ions['Tf2N'] = first_frame.topology.select('resname Tf2N')
#atoms_per_ion = dict()
#atoms_per_ion['EMIM'] = 19
#atoms_per_ion['Tf2N'] = 15

def calc_number_density(gro_file, trj_file, ions,
    bin_width, area, dim, box_range, data_path):
    
    for ion in np.sort(list(ions.keys())):
        #first_frame = md.load_frame(trj_file, top=gro_file,
        #    index=0, atom_indices=ions[ion])
        traj = md.load(trj_file, top=gro_file)#, atom_indices=ions[ion])
        print(traj)
        indices = [[at.index for at in compound.atoms] for compound in list(
            traj.topology.residues)]
        com = list()
        for i, ids in enumerate(indices):
            sub_traj = traj.atom_slice(ids)
            com.append(md.compute_center_of_mass(sub_traj))
        com = np.array(com)
        print(np.shape(com))
        x = np.histogram(com[:,1:,dim].reshape((-1, 1)),
            bins=np.linspace(box_range[0], box_range[1],
            num=1+round((box_range[1] - box_range[0])/bin_width)))
        np.savetxt('{0}/{1}-number-density.txt'.format(data_path, ion),
            np.vstack([x[1][:-1]+np.mean(x[1][:2])-box_range[0],
            x[0]/(area*bin_width*(len(traj)-1))]))
