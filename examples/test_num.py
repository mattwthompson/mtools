from mtools.post_process import calc_number_density

calc_number_density.calc_number_density(coord_file='aqueous-channel/equil.gro',
    trj_file='aqueous-channel/prod.xtc', resnames=['Na', 'Cl'],
    bin_width=0.01, area=5.112*5.166, dim=2, box_range=[0.682, 5.682],
    data_path='data')

fig = plt.figure(figsize=(6,4))

for ion in resnames:
    data = np.loadtxt('data/{}-number-density.txt'.format(ion))
    plt.plot([x-0.682 for x in data[:,0]], data[:,1], '-', label=ion)

plt.legend(loc=0)
plt.ylabel(r'Ion Number Densities $\left(\mathsf{\frac{molecule}{nm^3}}\right)$')
plt.xlabel('Position in channel, $nm$')
fig.savefig('img/number_density.pdf')
