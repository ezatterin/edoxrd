"""
Set of functions to simulate CTRs and RSMs. Updated 7.11.
"""

import numpy as np
import matplotlib.pyplot as plt
from pkg_resources import resource_filename
from edoxrd.utils import asf, mat_dict, tt2q, find_osc
from edoxrd.read import read_data, read_rsm_data

def calc_str_fac(material,ttheta):

	"""
	Calculate the Structure Factor of a material belonging to the Perovskite
	Oxides family. Materials: STO,PTO,MTO,SRO,KTO,DSO.

	Parameters
	----------
	material: string
	ttheta: ndarray

	Returns
	-------
	F: ndarray

	Example
	-------
	>>> F = calc_str_fac('PTO',xdata)

	"""

	# Initialise
	mat_param, mat_pos = mat_dict()
	q = tt2q(ttheta)
	f = asf(q)

	ls = dict(
	    STO = (['Sr', 'Ti'] + ['O']*3),
	    PTO = (['Pb', 'Ti'] + ['O']*3),
	    MTO = (['Mn', 'Ti'] + ['O']*3),
	    SRO = (['Sr', 'Ti'] + ['O']*3),
	    KTO = (['K', 'Ta'] + ['O']*3),
	    DSO = (['Dy', 'Sc'] + ['O']*3)
	        )

	if material == 'SRO_b' or material == 'SRO_t':
	    atom_list = ls['SRO']
	else:
	    atom_list = ls[material]

	F = 0
	for i, atom in enumerate(atom_list):
	    F += f[atom] * np.exp(1j * ( q * mat_pos[material][i,2]*mat_param[material][2]))

	return F

def calc_thickness(fname, threshold=1e-4, distance=20, side='r',
					comm='*'):

	"""
	Calculate the thickness of a thin film on top of a substrate from its
	Laue oscillations.

	Parameters
	----------
	fname: string
		The filename in .ras format.
	thres : float between [0., 1.]
	    Normalized threshold. Only the peaks with amplitude higher than the
	    threshold will be detected. Default: 1e-4
	min_dist : int
	    Minimum distance between each detected peak. The peak with the highest
	    amplitude is preferred to satisfy this constraint. Default: 20.
	side: string
		Side of the film peak where to look for oscillations. Either 'l'
		(left), or 'r' (right). Default: 'r'

	Returns
	-------
	t: float
		Thickness of the film.

	Example
	-------
	>>> t = calc_thickness(dset, d, threshold=1e-5, distance=15, side='l')

	"""

	xpeaks, ypeaks = find_osc(fname, threshold=threshold,
							  m_distance=distance, peak_side=side, comm=comm)
	x = np.asarray([x for x in range(len(xpeaks))]) # oscillation peak order
	y = 4 * np.pi * np.sin(np.deg2rad(xpeaks/2)) / 0.15406 # q_m's
	m, b = np.polyfit(x, y, 1) # linear fit
	t = abs(2 * np.pi / m)
	fig, ax = plt.subplots(1,2, figsize=(12,4))
	fig.suptitle('The thickness of sample {0} is {1:.3f} nm with accuracy {2:.3f}'\
	      .format(fname[:6], t, 1./x.max() ))
	ax[0].scatter(xpeaks, ypeaks,c='red')
	ax[0].plot(*read_data(fname, comments=comm)); ax[0].set_yscale('log')
	ax[1].scatter(x, y)
	ax[1].plot(x, x*m+b, c='red', label='$q_m = {0:.3f}m + {1:.3f}$'.format(m, b))
	ax[1].legend(); plt.show()
	return t

def l_prof_max(filelist, d, L=None, q=False):

	""" TODO! """

	# select correct dataset
	pto = filelist[1]

	# read data
	h, l, i = read_rsm_data(pto, d, scale='lin')

	if L==None:
	    L = l[np.where(i==i.max())]

	# find indexes of values close to selected L
	a = abs(l-L)
	mask_a = ((a > 0) & (a < 1e-4))
	lst_a = []
	for number in a[mask_a]:
	    lst_a.append([(int(np.where(a==number)[0])),int(np.where(a==number)[1])])
	idx_l = np.array(lst_a)

	I = i[idx_l[:,0], idx_l[:,1]][::-1]
	H = h[idx_l[:,0], idx_l[:,1]][::-1]

	if q == True:
	    H = (2*np.pi/0.399)*(H +1)

	return H, I
