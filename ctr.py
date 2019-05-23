import numpy as np
from pkg_resources import resource_filename
from edoxrd.utils import asf, mat_dict, tt2q
from edoxrd.read import read_data
from edoxrd.calc import calc_str_fac

def calc_ctr(fname, sub, film, Nfilm, c_film, Nelectrode_b=0, Nelectrode_t=0,
			c_electrode_b=0, c_electrode_t=0, scale=1e7, Nsub=1e4, c_sub=None,
			comm='*'):

	"""
	Calculate the thickness of a thin film on top of a substrate from its
	Laue oscillations.

	Parameters
	----------
	fname: string
		The filename in .asc format.
	datadir: string
		The directory containing fname.
	sub: string
		The name of the substrate material.
	film: string
		The name of the film material
	Nfilm: int
		The thickness of the film in unit cells.
	c_film: float
		The out-of-plane lattice parameter of the film in A.
	Nelectrode_b: int
		The thickness of the Bottom electrode (SRO_b) in unit cells.
		Default: 0 (No bottom Electrode).
	Nelectrode_t: int
		The thickness of the Top electrode (SRO_t) in unit cells.
		Default: 0 (No top Electrode).
	c_electrode_b: float
		The out-of-plane lattice parameter of the Bottom electrode in A.
		Default: 0 (No bottom Electrode).
	c_electrode_t: float
		The out-of-plane lattice parameter of the Top electrode in A.
		Default: 0 (No top Electrode).
	scale: int
		Overall multiplier for the calculated intensity.
		Default: 1e7.
	Nsub: int
		Thickness of the Substrate layer.
		Default: 1e4.

	Returns
	-------
	xdata: ndarray
		The original measured twotheta range.
	I: ndarray
		The calculated Intensity.

	Example
	-------
	>>> x, y = xrd.calc_ctr('e16014_09_t2t_002.ras', d, 'KTO', 'PTO', 60, 3.92)

	"""

	# Constants
	wave = 1.5406
	mu = 4e4
	elec_b = 'SRO_b'
	elec_t = 'SRO_t'
	param, pos = mat_dict()
	if c_sub is None:
		c_sub = param[sub][2]
	elif c_sub is not None:
		print('using c_sub as {0}, the "original" value is {1}'
				.format(c_sub, param[sub][2]))

	# Read the data
	xdata, ydata = read_data(fname)
	q = tt2q(xdata)

	# No electrodes
	if Nelectrode_b <=0 and Nelectrode_t <= 0:

		# Load data
		mats = [sub, film]
		t = (Nsub*c_sub + Nfilm*c_film)

		# Calc structure factor and build shape function dictionary
		F = {}
		for material in mats:
			F[material] = calc_str_fac(material, xdata) # does not want q, fix?
			S = {sub:0, film:0}
		# Calc shape function for each material
		for l in range(0, int(Nsub)):
			S[sub] += np.exp((t-l*c_sub) * (1j*q-(4*np.pi)/(mu*wave*q)))
		for l in range(0, int(Nfilm)):
			S[film] += np.exp((t - (l*c_film + Nsub*c_sub)) *
			(1j*q - (4*np.pi) / (mu*wave*q)))

	# Top electrode
	elif Nelectrode_b <= 0 and Nelectrode_t > 0:

		# Load data
		mats = [sub, elec_t, film]
		t = (Nsub*c_sub + Nelectrode_t*c_electrode_t + Nfilm*c_film)

		# Calc structure factor and build shape function dictionary
		F = {}
		for material in mats:
			F[material] = calc_str_fac(material, xdata)
			S = {sub:0, film:0, elec_t:0}

		# ?
		for l in range(0,int(Nsub)):
			S[sub] += np.exp((t-l*c_sub) * (1j*q-(4*np.pi)/(mu*wave*q)))
		for l in range(0, int(Nfilm)):
			S[film] += np.exp((t - (l*c_film + Nsub*c_sub)) *
			(1j*q - (4*np.pi) / (mu*wave*q)))
		for l in range(0, int(Nelectrode_t)):
			S[elec_t] += np.exp((t - (l*c_electrode_t + Nfilm*c_film
			+ Nsub*c_sub)) * (1j*q - (4*np.pi) / (mu*wave*q)))

	# Bottom electrode
	elif Nelectrode_b > 0 and Nelectrode_t <= 0:

		mats = [sub, elec_b, film]
		t = (Nsub*c_sub + Nelectrode_b*c_electrode_b + Nfilm*c_film)

		F = {}
		for material in mats:
			F[material] = calc_str_fac(material, xdata)
			S = {sub:0, film:0, elec_b:0}

		# ?
		for l in range(0,int(Nsub)):
			S[sub] += np.exp((t-l*c_sub) * (1j*q-(4*np.pi)/(mu*wave*q)))
		for l in range(0, int(Nelectrode_b)):
			S[elec_b] += np.exp((t - (l*c_electrode_b + Nsub*c_sub)) *
			(1j*q - (4*np.pi) / (mu*wave*q)))
		for l in range(0, int(Nfilm)):
			S[film] += np.exp((t - (l*c_film + Nelectrode_b*c_electrode_b
			+ Nsub*c_sub)) * (1j*q - (4*np.pi) / (mu*wave*q)))

	# Both electrodes
	elif Nelectrode_b > 0 and Nelectrode_t > 0:

		mats = [sub,film,elec_b,elec_t]
		t = (Nsub*c_sub + Nelectrode_b*c_electrode_b + Nfilm*c_film +
			 Nelectrode_t*c_electrode_t)

		F = {}
		for material in mats:
			F[material] = calc_str_fac(material, xdata)
			S = {sub:0, film:0, elec_b:0, elec_t:0}

		# ?
		for l in range(0,int(Nsub)):
			S[sub] += np.exp((t-l*c_sub) * (1j*q-(4*np.pi)/(mu*wave*q)))
		for l in range(0, int(Nelectrode_b)):
			S[elec_b] += np.exp((t - (l*c_electrode_b + Nsub*c_sub)) *
			(1j*q - (4*np.pi) / (mu*wave*q)))
		for l in range(0, int(Nfilm)):
			S[film] += np.exp((t - (l*c_film + Nelectrode_b*c_electrode_b
			+ Nsub*c_sub)) * (1j*q - (4*np.pi) / (mu*wave*q)))
		for l in range(0, int(Nelectrode_t)):
			S[elec_t] += np.exp((t - (l*c_electrode_t + Nfilm*c_film +
			Nelectrode_b*c_electrode_b + Nsub*c_sub)) *
			(1j*q - (4*np.pi) / (mu*wave*q)))


	g = 0
	for material in mats:
		g += F[material]*S[material]

	# Normalisation and Intensity
	g = g / g.max()
	I = scale * g * np.conj(g)

	return xdata, I
