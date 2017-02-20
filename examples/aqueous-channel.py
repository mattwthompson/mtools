import numpy as np
import matplotlib.pyplot as plt
from mtools.post_process import calc_charge_density
from mtools.post_process import calc_number_density

resnames = ['Na', 'Cl', 'HOH']
coord_file='aqueous-channel/prod.gro'
trj_file='aqueous-channel/prod.xtc'
top_file='aqueous-channel/prod.tpr'
bin_width=0.01
area=5.166*5.122
dim=2
box_range=[0.682, 5.682]
data_path='data/'

calc_number_density(coord_file, trj_file, bin_width, area, dim, box_range, data_path)

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

calc_charge_density(coord_file, trj_file, top_file, bin_width, area, dim, box_range, data_path)

fig, ax = plt.subplots(figsize=(6,4))

data = np.loadtxt('{}/charge-density.txt'.format(data_path))
ax.plot([x-0.682 for x in data[:,0]], data[:,1], '-')

ax.set_ylabel(r'Charge Density $\left(\mathsf{\frac{e^-}{nm^3}}\right)$')
ax.set_xlabel('Position in channel, $nm$')
fig.savefig('img/charge-density.pdf')
