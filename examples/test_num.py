from mtools.post_process import calc_number_density

calc_number_density.calc_number_density(coord_file='water-box/confout.pdb',
    trj_file='water-box/traj.xtc', ions={'SOL': ['OW', 'HW1', 'HW2']},
    bin_width=0.01, area=3.0129*3.0129, dim=0, box_range=[0, 3.0129],
    data_path='data')
