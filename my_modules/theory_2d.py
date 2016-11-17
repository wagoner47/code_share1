#! /usr/bin/env python

# A compilation of functions used to convert from a 3D power spectrum P(k, z) to a 2D correlation
# function w(theta) (or w(Rperp) in future versions). These function declarations can be combined
# for direct conversion or individually for conversion from precomputed Cl's

## Sub-function to be used in converting P(k, z) to Cl #1: integrand
## Inputs :
## z = redshift at which to evaluate
## distance = distance function/interpolator that takes as an argument z and returns the
## comoving distance
## select = the normalized redshift selection function/interpolator that takes z as an
## argument
## power = the function/interpolator for the 3D power P(k, z)
## l = the angular mode at which to evaluate
## cosmol = the cosmology object that contains H(z)
def power_integrand(z, distance, select, power, l, cosmol) :
	import numpy as np
	## define speed of light c
	c = 3.e5
	return ((cosmol.H(z).value/c)*np.power(select(z)/distance(z), 2)*power(l/distance(z), z))

## Sub-function for converting P(k, z) to Cl #2: integrator
## Inputs :
## l = angular mode at which to evaluate
## select = the normalized redshift selection function/interpolator that takes z as an
## argument
## power = the function/interpolator for the 3D power P(k, z)
## distance = distance function/interpolator that takes as an argument z and returns the
## comoving distance
## cosmol = the cosmology object that contains H(z)
## zmin (optional) = minimum redshift for integration (default: 0.0)
## zmax (optional) = maximum redshift for integration (default: 1.0)
def power_Cl(l, select, power, distance, cosmol, zmin=0.0, zmax=1.0) :
	from scipy.integrate import quad
	return (quad(power_integrand, zmin, zmax, args=(distance, select, power, l, cosmol))[0])

## P(k, z) to Cl function
## This function essentially calls the integration function in a loop over values of l,
## and then saves the output to a file and returns the array of Cl
## Inputs :
## table_l = numpy array of angular modes at which to evaluate
## select = the normalized redshift selection function/interpolator that takes z as an
## argument
## power = the function/interpolator for the 3D power P(k, z)
## distance = distance function/interpolator that takes as an argument z and returns the
## comoving distance
## cosmol = the cosmology object that contains H(z)
## path = path at which to save results
## zmin (optional) = minimum redshift for integration (default: 0.0)
## zmax (optional) = maximum redshift for integration (default: 1.0)
def Pk_to_Cl(table_l, select, power, distance, cosmol, path, zmin=0.0, zmax=1.0) :
	import numpy as np
	
	table_C = np.empty_like(table_l)
	for i, l in zip(range(table_l.size), table_l) :
		table_C[i] = power_Cl(l, select, power, distance, cosmol, zmin, zmax)
	
	np.savetxt(path, np.array([table_l, table_C]).T, fmt='%-25.18e', header='{:<25s} {:<25s}'.format('# l', 'Cl'), comments='')
	return table_C

## Cl to w(theta) function
## Inputs :
## nz = number of redshift bins (integer)
## table_l = table of l values at which Cl is defined
## table_Cl = table of Cl's with shape (nz, nl)
## theta = table of thetas in degrees
## path = file in which to save the results
def Cl_to_wtheta(nz, table_l, table_Cl, theta, path, nu=0, N=500, h=0.005) :
	assert isinstance(nz, int), 'Invalid nz: {}. Must be an integer'.format(nz)
	import os
	from hankel import HankelTransform as HT
	import numpy as np
	from scipy.interpolate import UnivariateSpline
	
	header = '{:<25s}'.format('# r (deg)')
	for i in range(nz) :
		header += '{:<25s}'.format('w_{}'.format(i))
	
	ht = HT(nu=nu, N=N, h=h)
	
	w = np.empty((nz, theta.size))
	for i in range(nz) :
		# Fit Cl's with Spline
		Cl = UnivariateSpline(table_l, table_Cl[i], s=0, k=1)

		# Do the Hankel Transform to get the correlation function
		for j, thetai in zip(range(theta.size), np.deg2rad(theta)) :
			f = lambda x: (x*Cl(x/thetai))/(2.*np.pi*np.power(thetai,2))
			w[i,j] = ht.transform(f)[0]
	
	np.savetxt(path, np.vstack((theta, w)).T, fmt='%-25.18e', header=header, comments='')
	return w