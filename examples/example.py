import numpy as np
import matplotlib.pyplot as plt
from mtools.post_process import calc_number_density
from mtools.post_process import calc_charge_density

resnames = ['Na', 'Cl', 'HOH']
coord_file='aqueous-channel/prod.gro'
trj_file='aqueous-channel/prod.xtc'
bin_width=0.05
area=5.166*5.122
dim=2
box_range=[0.682, 5.682]
data_path='data/'

calc_number_density.calc_number_density(coord_file, trj_file,
    resnames, bin_width, area, dim, box_range, data_path)

fig = plt.figure(figsize=(6,4))
           
for res in resnames:
    data = np.loadtxt('{}/{}-number-density.txt'.format(data_path,res))
    if res == 'HOH':
        plt.plot([x-0.682 for x in data[:,0]],
            [x/20 for x in data[:,1]], '-', label='Water/20')
    else:
        plt.plot([x-0.682 for x in data[:,0]], data[:,1], '-', label=res)

plt.legend(loc=0)
plt.ylabel(r'Ion Number Densities $\left(\mathsf{\frac{molecule}{nm^3}}\right)$')
plt.xlabel('Position in channel, $nm$')
fig.savefig('img/number_density.pdf')

calc_charge_density.calc_charge_density(coord_file, trj_file,
    itp_path='aqueous-channel/', resnames=['Na', 'Cl'],
    bin_width=bin_width, area=area, dim=dim, box_range=box_range,
    data_path=data_path)

fig = plt.figure(figsize=(6,4))

data = np.loadtxt('{}/charge-density.txt'.format(data_path))
plt.plot([x-0.682 for x in data[:,0]], data[:,1], '-')

plt.legend(loc=0)
plt.ylabel(r'Charge Density $\left(\mathsf{\frac{e^-}{nm^3}}\right)$')
plt.xlabel('Position in channel, $nm$')
fig.savefig('img/charge-density.pdf')
