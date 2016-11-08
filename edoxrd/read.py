"""
Set of functions to read CTR and RSM data as generated by a Rigaku SmartLab
Diffractometer. Updated 7.11.
"""

import numpy as np
import pandas as pd

def read_data(fname, datadir,comments='*'):
	"""
	Reads .ras Rigaku files.

	Parameters
	----------
	fname: string
		Filename in .ras extension.
	datadir: string
		Directory containing fname.

	Returns
	-------
	xdata: ndarray
		Array of Two Theta values.
	ydata: ndarray
		Array of Intensity values.

	Example
	-------
	>>> xdata, ydata = read_data('e16014_01_t2t_001.ras','~/Desktop/')

	"""

	# Define path to file
	path = datadir + fname

	# Read it into an array
	data = np.genfromtxt(path, delimiter=" ", comments=comments)
	xdata = data[:,0] # 2theta in degrees
	ydata = data[:,1] # Intensity
	offset = data[:,2] # Offset

	# Correct for detector attenuation
	if np.any(offset > 1):
		mask = offset > 1
		indeces = np.where(mask)
		ydata[indeces] = ydata[indeces] * offset[indeces]

	return xdata, ydata

def read_rsm_data(fname, datadir, scale='log', coordinates='hkl'):

	"""
	Read data from a Reciprocal Space Map, in .asc format.

	Parameters
	----------
	fname: string
		Filename in .asc extension.
	datadir: string
		Directory containing fname.
	scale: string, optional
		Scale of Intensity values, linear ('lin') or logarithmic ('log').
		Default is 'log'.
	coordinates: string
		Coordinate frame for Intensity mapping. Can be in 2theta, omega space
		('ttomega') or reciprocal space ('hkl'). Default is 'hkl'.

	Returns
	-------
	h: ndarray
		Array of H values.
	l: ndarray
		Array of L values.
	I: ndarray
		Array of Intensity values scaled as specified.

	Example
	-------
	>>> h,k,I = read_rsm_data('e16019_01_-103_KTO_RSM_2-Theta.asc',d)

	"""

	# full path to the data file
	peak = datadir + fname
	# read data as lines
	d = np.genfromtxt(peak, dtype='str', delimiter='\n')

	# Wavelength
	wave = [line for line in d if '*WAVE_LENGTH1' in line]
	lambdaone = float(wave[0].split()[2])

	# Scan speed
	speed_name = [line for line in d if '*SPEED' in line]
	speed = int(speed_name[1].split()[2])

	# build array of measured Twotheta values
	start = [line for line in d if '*START' in line]
	twotheta_start = float(start[0].split()[2])
	stop = [line for line in d if '*STOP' in line]
	twotheta_stop = float(stop[0].split()[2])
	step = [line for line in d if '*STEP' in line]
	twotheta_step = float(step[0].split()[2])
	count = [line for line in d if '*COUNT' in line]
	no_twotheta = int(count[1].split()[2])  # there is a 'counter' line
	# 1D array of meas twotheta
	twotheta = np.arange(twotheta_start, twotheta_stop, twotheta_step)

	# build array of measured Omega values
	sec_count = [line for line in d if '*SEC_COUNT' in line]
	no_omega = int(sec_count[0].split()[2])
	omega_lines = [line for line in d if '*OFFSET' in line]
	# 1D array of meas omega
	omega = np.zeros(len(omega_lines, ))
	for index, line in enumerate(omega_lines):
		omega[index] = line.split()[2]

	# build array of measured Intensity values
	d2 = pd.read_csv(peak, header=None, comment='*')
	d2 = d2.values.flatten()  # 1D nparray of I
	d2 = d2[np.logical_not(np.isnan(d2))]  # delete NaN values
	d2 = d2.reshape(no_twotheta, no_omega, order='F')  # matlab order
	d2 = d2 / speed # Normalise to counts per second
	if scale=='log':
		I = np.log(d2)
	elif scale=='lin':
		I = d2
	I[I<0] = 0.0 # avoid log scale from making negative intensities

	if coordinates=='ttomega':
		return twotheta, omega, I

	elif coordinates=='hkl':

		# make mesh
		xx, yy = np.meshgrid(twotheta, omega, indexing='ij')
		# to radiants
		ttrad = np.deg2rad(xx)
		omrad = np.deg2rad(yy)

		# to h, l
		h = 3.988 / lambdaone * (np.cos(omrad) - np.cos(ttrad - omrad))
		l = 3.988 / lambdaone * (np.sin(omrad) + np.sin(ttrad - omrad))

		return h, l, I
