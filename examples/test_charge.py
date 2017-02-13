from mtools.post_process import calc_charge_density

calc_charge_density.calc_charge_density(coord_file='aqueous-channel/em.gro',
    trj_file = 'aqueous-channel/equil.xtc', itp_path = 'aqueous-channel/',
    resnames=['Na', 'Cl'], bin_width=0.01, area=5.166*5.112,
    dim=2, box_range=[0.682, 5.682], data_path='data')
