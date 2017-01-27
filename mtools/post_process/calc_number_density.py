import numpy as np
import mdtraj as md

def calc_number_density(coord_file, trj_file, ions,
    bin_width, area, dim, box_range, data_path):
    
    for ion in np.sort(list(ions.keys())):
        #first_frame = md.load_frame(trj_file, top=coord_file,
        #    index=0, atom_indices=ions[ion])
        traj = md.load(trj_file, top=coord_file)#, atom_indices=ions[ion])
        indices = [[at.index for at in compound.atoms] for compound in list(
            traj.topology.residues)]
        com = list()
        for i, ids in enumerate(indices):
            sub_traj = traj.atom_slice(ids)
            com.append(md.compute_center_of_mass(sub_traj))
        com = np.array(com)
        x = np.histogram(com[:,1:,dim].reshape((-1, 1)),
            bins=np.linspace(box_range[0], box_range[1],
            num=1+round((box_range[1] - box_range[0])/bin_width)))
        np.savetxt('{0}/{1}-number-density.txt'.format(data_path, ion),
            np.vstack([x[1][:-1]+np.mean(x[1][:2])-box_range[0],
            x[0]/(area*bin_width*(len(traj)-1))]))
