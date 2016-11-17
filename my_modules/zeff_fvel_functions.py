#! /usr/bin/env python

# Functions to help with calculating the effective redshift and logarithmic derivative
# of the growth function (f_vel) for a sample of galaxies with redshifts, or just f_vel
# for a sample of galaxies/halos/DM particles that have a known effective redshift. I'll 
# also include one for getting the average bias

def distance_interpolator(z_max=2.0, nz=101, **cosmology) :
	import numpy as np
	from scipy.interpolate import UnivariateSpline
	from astropy.cosmology import w0waCDM
	## This is setting up our cosmology
	cosmo = w0waCDM(100.*cosmology['h0'], cosmology['omega_m'], \
	1.-cosmology['omega_k']-cosmology['omega_m'], w0=cosmology['w'], wa=cosmology['wa'])
	
	## Set up a grid of redshifts at which to evaluate the distance to set up the spline
	table_z = np.linspace(0.0, z_max, num=nz)
	table_r = cosmo.comoving_distance(table_z).value
	
	## Create and return the spline object. Note that a cubic spline works better than
	## a linear spline here
	distance = UnivariateSpline(table_z, table_r, s=0)
	return distance

def Hubble_interpolator(z_max=2.0, nz=101, **cosmology) :
	import numpy as np
	from scipy.interpolate import UnivariateSpline
	from astropy.cosmology import w0waCDM
	## This is setting up our cosmology
	cosmo = w0waCDM(100.*cosmology['h0'], cosmology['omega_m'], \
	1.-cosmology['omega_m']-cosmology['omega_k'], w0=cosmology['w'], wa=cosmology['wa'])
	
	## Set up a grid of redshifts at which to evaluate the Hubble parameter to set up 
	## the spline
	table_z = np.linspace(0.0, z_max, num=nz)
	table_H = cosmo.H(table_z).value
	
	## Create and return the spline object. Note that a cubic spline works better than
	## a linear spline here
	Hofz = UnivariateSpline(table_z, table_H, s=0)
	return Hofz

## This one needs real data: the redshifts and error(s) of your sample in order to 
## generate the redshift distribution
def phi_interpolator(data_z, data_zerr, z_max=2.0, nz=101) :
	import numpy as np
	from scipy.interpolate import UnivariateSpline
	from astropy.modeling.models import Gaussian1D
	
	assert hasattr(data_z, "__len__") and not isinstance(data_z, str), \
	'Data redshifts must be array-like'
	data_z = np.array(data_z)
	assert isinstance(data_zerr, float) or (hasattr(data_zerr, "__len__") and not \
	isinstance(data_zerr, str)), 'Data redshift error(s) must either be scalar float or \
	array-like'
	if hasattr(data_zerr, "__len__") :
		data_zerr = np.array(data_zerr)
		assert data_zerr.size == data_z.size, 'Data redshift errors must have same \
		length as data redshifts if not scalar float'
	
	## Set up the Gaussians for the redshifts
	if isinstance(data_zerr, np.ndarray) :
		amp = 1./np.sqrt(2.*np.pi*np.power(data_zerr.copy(), 2))
		g = Gaussian1D(amplitude=amp.copy(), mean=data_z.copy(), stddev=data_zerr.copy())
		del amp
	else :
		amp = 1./np.sqrt(2.*np.pi*np.power(data_zerr, 2))
		g = Gaussian1D(amplitude=amp, mean=data_z.copy(), stddev=data_zerr)
		del amp
	
	## Set up a grid of redshifts at which to evaluate the redshift distribution to set
	## up the spline. We will also need the bin size to normalize properly.
	table_z, delta_z = np.linspace(0.0, z_max, num=nz, retstep=True)
	table_g = np.array([np.sum(g(zi)) for zi in table_z])
	table_g /= (table_g.sum()*delta_z)
	
	## Create and return the spline object. Note that a linear spline is safer here.
	phi = UnivariateSpline(table_z, table_g, k=1, s=0)
	return phi

def zeff_numerator_integrand(z, H, xi, phi, interp_flag=True) :
	import numpy as np
	## This one is straight forward. However, if we don't use splines, we need to get 
	## just the values
	if interp_flag :
		return z*np.power(phi(z), 2)*(H(z)/np.power(xi(z), 2))
	else :
		return z*np.power(phi(z), 2)*(H(z).value/np.power(xi(z).value, 2))

def zeff_denominator_integrand(z, H, xi, phi, interp_flag=True) :
	import numpy as np
	## This one is straight forward. However, if we don't use splines, we need to get 
	## just the values
	if interp_flag :
		return np.power(phi(z), 2)*(H(z)/np.power(xi(z), 2))
	else :
		return np.power(phi(z), 2)*(H(z).value/np.power(xi(z).value, 2))

## Adding the variable interp_flag to determine whether we are okay with simply using 
## the splines for the calculation or if we want to use the exact functions. The 
## exception is the redshift distribution, can only be done with a spline anyways
def calc_zeff(data_z, data_zerr, table_z_max=2.0, nz=101, z_min=0.6, z_max=0.8, \
interp_flag=True, **cosmology) :
	import numpy as np
	from scipy.integrate import quad
	
	assert table_z_max > z_max, 'To avoid extrapolation, table_z_max must be larger \
	than z_max: table_z_max={}, z_max={}'.format(table_z_max, z_max)
	
	## Get the redshift distribution
	phi = phi_interpolator(data_z, data_zerr, table_z_max, nz)
	
	if interp_flag :
		## In this case, get the distance and Hubble splines
		xi = distance_interpolator(table_z_max, nz, **cosmology)
		H = Hubble_interpolator(table_z_max, nz, **cosmology)
		## Calculate the numerator and denominator
		zeff_num = quad(zeff_numerator_integrand, z_min, z_max, args=(H, xi, phi, \
		interp_flag))[0]
		zeff_denom = quad(zeff_denominator_integrand, z_min, z_max, args=(H, xi, phi, \
		interp_flag))[0]
	else :
		## Here we need to set up the cosmology instead
		from astropy.cosmology import w0waCDM
		cosmo = w0waCDM(100.*cosmology['h0'], cosmology['omega_m'], \
		1.-cosmology['omega_k']-cosmology['omega_m'], w0=cosmology['w'], \
		wa=cosmology['wa'])
		## Calculate the numerator and denominator
		zeff_num = quad(zeff_numerator_integrand, z_min, z_max, args=(cosmo.H, \
		cosmo.comoving_distance, phi, interp_flag))[0]
		zeff_denom = quad(zeff_denominator_integrand, z_min, z_max, args=(cosmo.H, \
		cosmo.comoving_distance, phi, interp_flag))[0]
	
	## Now we just return zeff = zeff_num/zeff_denom
	return zeff_num/zeff_denom

## Here we will just assume there is a look-up table for fvel evaluated at different
## redshifts that we can read in and create a spline from
def calc_fvel(fvel_lookup_path, zeff) :
	import numpy as np
	from scipy.interpolate import UnivariateSpline
	
	## Read in the grid of z and fvel
	table_z, table_f = np.loadtxt(fvel_lookup_path, unpack=True)
	
	## Create the spline object. I don't know if this is better linear or cubic yet,
	## so I'll assume linear until I double check
	fvel = UnivariateSpline(table_z, table_f, k=1, s=0)
	
	## Return the value at the given effective redshift
	return fvel(zeff)

## Some helper functions for calculating the bias. First the numerator integrand
def bias_numerator_integrand(logm, a, delta, cd) :
	import cosmocalc
	import numpy as np
	## Set the cosmocalc cosmology
	cosmocalc.set_cosmology(cd)
	## Need the mass in the integral, not the log10(mass)
	m = np.power(10., logm)
	## Return the integrand
	return cosmocalc.tinker2010_mass_function(m, a, delta)*\
	cosmocalc.tinker2010_bias(m, a, delta)*m

## Now the denominator integrand
def bias_denominator_integrand(logm, a, delta, cd) :
	import cosmocalc
	import numpy as np
	## Set the cosmocalc cosmology
	cosmocalc.set_cosmology(cd)
	## Need the mass in the integral, not the log10(mass)
	m = np.power(10., logm)
	## Return the integrand
	return cosmocalc.tinker2010_mass_function(m, a, delta)*m

## Here we'll calculate the average bias at a redshift above a given mass threshold. Note
## that min_mass is given as the base 10 log of the minimum mass in Msun/h
def calc_bias(zeff, min_mass=13.5, **cosmology) :
	import numpy as np
	from scipy.integrate import quad
	
	## Here's the cosmocalc cosmology dictionary
	cd = {"om" : cosmology['omega_m'], "ob" : cosmology['omega_b'], \
	"ol" : 1.-cosmology['omega_k']-cosmology['omega_m'], "ok" : cosmology['omega_k'], \
	"h" : cosmology['h0'], "s8" : cosmology['sigma8_input'], "ns" : cosmology['n_s'], \
	"w0" : cosmology['w'], "wa" : cosmology['wa']}
	
	## cosmocalc uses a(z) rather than redshift, so calculate it for our zeff
	a = 1./(1. + zeff)
	
	## Calculate the numerator and denominator with delta=200
	bias_num = quad(bias_numerator_integrand, min_mass, 16., args=(a, 200, cd))[0]
	bias_denom = quad(bias_denominator_integrand, min_mass, 16., args=(a, 200, cd))[0]
	
	## Now return the bias = bias_num/bias_denom
	return bias_num/bias_denom