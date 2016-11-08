"""
Set of functions to perform a variety of tasks related to Diffraction data
analysis. Updated 11.7.
"""

import numpy as np
import peakutils as pk
from pkg_resources import resource_filename
from edoxrd.read import read_data


def mat_dict():

	"""
	Generates a tuple containing dictionaries with Unit Cell Parameters and
	Atomic coordinates for a range of materials.

	Returns
	-------
	out: tuple
	A tuple of two dictionaries, Material_name:Parameters and
	Material_name:Atomic_coordinates

	Example
	-------
	>>> param, pos = mat_dict()
	>>> print(param['PTO'])

	"""

	# Load data files
	mat_param = resource_filename('edoxrd','/data/material_parameters.txt')
	atom_positions = resource_filename('edoxrd', '/data/atom_positions.txt')

	# Generate arrays of parameters and material names
	param = np.genfromtxt(mat_param, usecols=(0,1,2))
	names = list(np.genfromtxt(mat_param, usecols=3, dtype='|S5'))

	# Make a dictionary Mat_name:Parameters(a,b,c)
	mat_param = {}
	for n, a in zip(names, param):
	    mat_param[n] = a

	# Generate array of atomic positions for each material
	pos = np.genfromtxt(atom_positions).reshape(((len(mat_param)+1),5,3))

	# Make a dicitonary Mat_name:Positions(x,y,z)
	mat_pos = {}
	for name, array in zip(mat_param.iterkeys(), pos):
	    mat_pos[name] = array

	return mat_param, mat_pos

def tt2q(ttheta):

	"""
	Convert TwoTheta values to Q (Scattering vector).

	Parameters
	----------
	ttheta: 1D array of TwoTheta angles

	Returns
	-------
	out: ndarray
	1D Array of Q values

	Example
	-------
	>>> Qx = tt2q(two_theta_x)
	"""

	q = 4 * np.pi * np.sin(np.deg2rad(ttheta/2))/1.5406

	return q

def asf(q):
	"""
	Function to calculate the Atomic Scattering Factor of several atoms.

	Parameters
	----------
	q: ndarray
		Array of scattering vectors.

	Returns
	-------
	out: dictionary
		Each entry of the output dictionary contains the ASF (ndarray) for the specified
		material key.

	Example
	-------
	>>> F_Ti = asf(q)['Ti']
	"""

	# Load data
	filepath = resource_filename('edoxrd','/data/asf_atom_constants.txt')

	# Generate list of scattering constants for each atom
	atom_const = np.genfromtxt(filepath)
	atom_list = ['Mn', 'O', 'Pb', 'Ru', 'Sr', 'Ti', 'K', 'Ta', 'Mg', 'Dy', 'Sc']
	atom_asf = {}
	for name, array in zip(atom_list, atom_const):
	    atom_asf[name] = array

	# Calc the ASF for each atom and append in dictionary
	fmat = {}
	for name, m in atom_asf.iteritems(): #each m is the atom array
	    f = m[0] * np.exp(-m[1]*q/4/np.pi) +\
	        m[2] * np.exp(-m[3]*q/4/np.pi) +\
	        m[4] * np.exp(-m[5]*q/4/np.pi) +\
	        m[6] * np.exp(-m[7]*q/4/np.pi) +\
	        m[8]
	    fmat[name] = f

	return fmat

def find_osc(fname, d, threshold=0.000006, m_distance=10, peak_side='r',
			comm='*'):
	"""TODO!"""

	xdata, ydata = read_data(fname,d,comments=comm)
	idxs = pk.indexes(ydata, thres=threshold, min_dist=m_distance)
	film_peak = ydata[idxs].argsort()[::-1][1]
	peaks = ydata[idxs].argsort()[::-1][2:]
	xpeaks, ypeaks = xdata[idxs][peaks],ydata[idxs][peaks]
	if peak_side == 'r':
	    xpeaks, ypeaks = xpeaks[xpeaks > xdata[idxs][film_peak]], ypeaks[xpeaks > xdata[idxs][film_peak]]
	else:
	    xpeaks, ypeaks = xpeaks[xpeaks < xdata[idxs][film_peak]], ypeaks[xpeaks < xdata[idxs][film_peak]]
	return xpeaks[0:12], ypeaks[0:12]