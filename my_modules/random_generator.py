
# coding: utf-8

# # Random Generator
# 
# ## Updated: June 10, 2015
# 
# ### [Erika L Wagoner](mailto:wagoner47@email.arizona.edu 'Email')
# 
# This file contains functions useful for creating a random catalog from a HEALPix map for a given data catalog. It can be called in scripts like a module to run the various processes needed to create the random catalog.
# 
# ##### Update on June 10, 2015
# To attempt to make the density plots, Eduardo and I have devised a new scheme. We will create our own pixels (using correct x and y values to preserve area). On each healpix pixel, we will calculate the weighted number of points (N$_{\mathrm{pix}}$/FRACGOOD), and on each of our larger pixels, we will add up the weighted counts of all healpix pixels that lie within it and divide by the sum of the pixel areas to get the density on our larger pixels. Using this method, we may be able to use less memory and create an array version of the pixel density that can be used in color map type matplotlib functions.

# #### HEALPix mask: FRACGOOD
# 
# Given the location of a HEALPix file, extension, and Boolean for header (used in astropy.io.fits.getdata), create a mask that has FRACGOOD as the value of the mask at all pixels with observations and 0 everywhere else. This mask would be useful for the type of map that uses 1-FRACGOOD points at each pixel rather than each pixel simply being good or bad. The other values used in generating the random catalog will also be read out here while the HEALPix file is being read.
# 
# The function will return (in this order):
# 
# * nside
# * nest
# * mask = The array holding the FRACGOOD for all pixels
# * elapsed time
# 
# Default value for extension is ext=1, and for header is header=True, so these are optional.

# In[2]:

def hp_mask_fracgood(hp_file_path, ext=1, header=True) :
	# Import modules needed
	import numpy as np
	from astropy.io import fits
	import healpy as hp

	# Read data from HEALPix file
	stg, hdr = fits.getdata(hp_file_path, ext=ext, header=header)

	# Extract NSIDE and NEST
	nside = hdr['NSIDE']
	nest = hdr['NEST']

	# Initialize mask: array of floating point 0's, shape should be healpy.nside2npix(nside)
	mask = np.zeros(hp.nside2npix(nside), dtype=np.float32)

	# Set FRACGOOD for observed pixels.
	mask[stg['HPIX']] = stg['FRACGOOD']

	return (nside, nest, mask)


# #### Generate points
# 
# Given location of file with galaxy data for which randoms are being created and multiplier N, generate N$\cdot$(size of data file) random points uniformly distributed on the surface of a unit sphere.
# 
# Returned values are, in order:
# 
# * ra = array of random $\alpha$'s
# * dec = array of random $\delta$'s
# * size of data array
# * elapsed time
# 
# Note: See [Wolfram Mathworld Sphere Point Picking page](http://mathworld.wolfram.com/SpherePointPicking.html "Wolfram Mathworld: Sphere Point Picking") for random point generation algorithm.

# In[3]:

def generate_points(N, seed=0) :
    # Import modules
    from astropy.io import fits
    import numpy as np
    import sys
    
    # Generate num_points random points in (0,1) for u and v
    np.random.seed(0)
    np.random.seed(seed)
    u = np.random.rand(N)
    while np.count_nonzero(u) != u.size :
        u[(u == 0)] = np.random.rand(u[(u == 0)].size)
    v = np.random.rand(N)
    while np.count_nonzero(v) != v.size :
        v[(v == 0)] = np.random.rand(v[(v == 0)].size)
    
    # Determine RA and DEC from u and v
    ra = np.degrees(2.0*np.pi*u)
    dec = 90.0 - np.degrees(np.arccos(2.0*v - 1.0))
    
    # Return values
    return (ra, dec)


# #### Apply mask
# 
# Given the size of the data file, nside, nest, mask, $\alpha$'s, and $\delta$'s, apply the mask to the random points, and then determine the ratio of the number of random points to the number of data points.
# 
# Returned values are, in order:
# 
# * cut_ra = array of $\alpha$'s after applying mask
# * cut_dec = array of $\delta$'s after applying mask
# * ratio = ratio of number of random points after applying mask to number of data points
# * elapsed time

# In[4]:

def apply_mask(size, hp_file_path, ext=1, header=True, seed=0, ramin=0., ramax=360., dmin=-90., dmax=90.) :
	# Import modules
	import numpy as np
	import healpy as hp
	import sys

	nside, nest, mask = hp_mask_fracgood(hp_file_path, ext=ext, header=header)
	npix = hp.nside2npix(nside)

	# Generate npoints=10*size random points to begin with. We'll keep adding more points as we need them.
	factor = 10
	npoints=factor*size
	ra, dec = generate_points(npoints, seed)	
	if ramin > ramax :
		rafilt = ((ra > ramin) | (ra < ramax))
	else :
		rafilt = ((ra > ramin) & (ra < ramax))
	decfilt = ((dec > dmin) & (dec < dmax))

	# Convert to theta and phi
	phi = np.radians(ra[(rafilt & decfilt)].copy())
	theta = np.radians(90.0 - dec[(rafilt & decfilt)].copy())
	del ra, dec, rafilt, decfilt

	# Find pixel number for all angular coordinates
	pix = hp.ang2pix(nside, theta, phi, nest=nest)

	# Create filter for good pixels from mask
	good, = np.where(mask[pix] > 0)

	# Redefine theta, phi, and pix to only include points on good pixels, and extract FRACGOOD for each point
	theta = theta[good]
	phi = phi[good]
	pix = pix[good]
	fracgood = mask[pix].copy()

	# Remove 1-fracgood points from each pixel by keeping only points with random number x in (0,1) where x < fracgood
	np.random.seed(0)
	np.random.seed(seed)
	x = np.random.rand(pix.size)
	cut_theta = theta[np.where(x < fracgood)].copy()
	cut_phi = phi[np.where(x < fracgood)].copy()
	del theta, phi, x, pix, good, fracgood

	# While we still have less than size points remaining, I'll just generate 100*(size-cut_theta.size) more points 
	# to run on
	while (cut_theta.size < size) :
		# factor *= 10
		f = float(cut_theta.size)/float(size)
		remsize = size - cut_theta.size
		npoints = (factor*remsize)/f
		if npoints < 20 :
			npoints *= factor
		ra, dec = generate_points(npoints, seed)	
		if ramin > ramax :
			rafilt = ((ra > ramin) | (ra < ramax))
		else :
			rafilt = ((ra > ramin) & (ra < ramax))
		decfilt = ((dec > dmin) & (dec < dmax))
	
		phi = np.radians(ra[(rafilt & decfilt)].copy())
		theta = np.radians(90.0 - dec[(rafilt & decfilt)].copy())
		del ra, dec, rafilt, decfilt
	
		pix = hp.ang2pix(nside, theta, phi, nest=nest)
	
		good, = np.where(mask[pix] > 0)
	
		theta = theta[good]
		phi = phi[good]
		pix = pix[good]
		fracgood = mask[pix].copy()
	
		np.random.seed(0)
		np.random.seed(seed)
		x = np.random.rand(pix.size)
		cut_theta = np.append(cut_theta, theta[np.where(x < fracgood)].copy())
		cut_phi = np.append(cut_phi, phi[np.where(x < fracgood)].copy())
		del theta, phi, x, pix, good, fracgood

	# Trim to be only size number of points
	if cut_theta.size > size :
		cut_theta = cut_theta[:size]
		cut_phi = cut_phi[:size]

	# Convert back to RA and DEC
	cut_ra = np.degrees(cut_phi.copy())
	cut_dec = 90.0 - np.degrees(cut_theta.copy())

	return (cut_ra, cut_dec)

def make_points(size, ramin, ramax, dmin, dmax, seed=0) :
	# Import modules
	import numpy as np
	import sys

	# Generate npoints=10*size random points to begin with. We'll keep adding more points as we need them.
	factor = 10
	npoints=factor*size
	ra, dec = generate_points(npoints, seed)
	
	if ramin > ramax :
		rafilt = ((ra > ramin) | (ra < ramax))
	else :
		rafilt = ((ra > ramin) & (ra < ramax))
	decfilt = ((dec > dmin) & (dec < dmax))

	cut_ra = ra[(rafilt & decfilt)].copy()
	cut_dec = dec[(rafilt & decfilt)].copy()
	del ra, dec
	
	while (cut_ra.size < size) :
		# factor *= 10
		f = float(cut_ra.size)/float(size)
		remsize = size - cut_ra.size
		npoints = (factor*remsize)/f
		if npoints < 20 :
			npoints *= factor
		ra, dec = generate_points(npoints, seed)
	
		if ramin > ramax :
			rafilt = ((ra > ramin) | (ra < ramax))
		else :
			rafilt = ((ra > ramin) & (ra < ramax))
		decfilt = ((dec > dmin) & (dec < dmax))
	
		cut_ra = np.append(cut_ra, ra[(rafilt & decfilt)].copy())
		cut_dec = np.append(cut_dec, dec[(rafilt & decfilt)].copy())
		del ra, dec

	# Trim to be only size number of points
	if cut_ra.size > size :
		cut_ra = cut_ra[:size]
		cut_dec = cut_dec[:size]

	return (cut_ra, cut_dec)


# #### Find pixel densities
# 
# Given the $\alpha$'s and $\delta$'s as well as nside, nest, and possibly the mask, find the density on each pixel.
# 
# Returned objects are, in order:
# 
# * array of pixel densities
# * array of all pixel $\alpha$'s
# * array of all pixl $\delta$'s
# * elapsed time

# In[ ]:

def pixel_count(nside, nest, ra, dec, mask) :
    # Import modules
    import numpy as np
    import healpy as hp
    
    # Convert ra and dec to theta and phi
    phi = np.radians(ra.copy())
    theta = np.radians(90. - dec.copy())
    
    # Find which pixels have points
    pix = hp.ang2pix(nside, theta, phi, nest=nest)
    fracgood = mask[pix].copy()
    npix = hp.nside2npix(nside)
    # pix_area = hp.nside2pixarea(nside, degrees=True)
    
    # Find all pixels
    all_pix = np.arange(npix)
    
    # Count number of points on pixels
    pix_count = np.bincount(pix, minlength=npix, weights=1./fracgood)
    
    # If mask is supplied, correct pixel densities for fracgood
	# if not mask == [] :
	# 	for i, fg_i, pc_i in zip(all_pix, mask, pix_count) :
	# 		if not fg_i == 0 :
	# 			pc_i = pc_i/fg_i
	# 		else :
	# 			if not pc_i == 0 :
	# 				print (i, pc_i)
	# 			pc_i = 0
	# 		pix_count[i] = pc_i
    
    # Find coordinates of all pixels
    all_theta, all_phi = hp.pix2ang(nside, all_pix, nest=nest)
    all_ra = np.degrees(all_phi)
    all_dec = 90. - np.degrees(all_theta)
    
    # Return objects
    return (pix_count, all_ra, all_dec)


# #### Plot random point densities
# 
# Given $\alpha$'s and $\delta$'s, create a 2D histogram of random galaxies. Title should also be given as an input. Labels for the x and y axes can optionally be given, with defaults being $\alpha\,\cos(\delta)$ on the x-axis and $\delta$ on the y-axis in degrees.
# 
# I also want to have the min and max for both RA (xmin, xmax) and DEC (ymin, ymax) passed in. These values should be the desired min and max for the plot, not the actual min and max for the data. I will use xmin and xmax to determine the center and tick labels on the x-axis. There should be three cases:
# 
# * 0 < xmin < xmax - Axis centered between xmin and xmax if 0 is included in the range on the y-axis, or centered in the largest range of xmin$\,\cdot\,\cos$(y) if not. Tick labels are the actual angles.
# * xmin < xmax < 0 - Same as case above, but tick labels are given by x+360$^o$ instead.
# * xmin < 0 < xmax - Axis centered at 0, but tick labels for x < 0 are shifted to x+360$^o$ while those for x > 0 are left as is.
# 
# Returned objects are, in order:
# 
# * The figure object
# * elapsed time

# In[5]:

def rand_point_density(axis, ra, dec, nbins, hp_file_path, ext=1, header=True, xmin=0., xmax=360., ymin=-1., ymax=1., vmax=None) :
	# Import modules
	import numpy as np
	import matplotlib.pyplot as plt
	import healpy as hp

	nside, nest, mask = hp_mask_fracgood(hp_file_path, ext=ext, header=header)

	# Find pixel densities using pixel_density function
	pix_area = hp.nside2pixarea(nside, degrees=True)
	pix_count, all_ra, all_dec = pixel_count(nside, nest, ra=ra, dec=dec, mask=mask)
	if xmin < 0 :
		all_ra[np.where(all_ra > xmax)] -= 360

	# Find widths of each self defined pixel
	nedges=nbins+1
	deltax = (float(xmax) - float(xmin))/float(nbins)
	deltay = (float(ymax) - float(ymin))/float(nbins)

	x = all_ra.copy()
	y = np.sin(np.deg2rad(all_dec.copy()))
	# Now find the densities on each of these larger pixels
	counts = np.zeros((nbins, nbins))
	area = np.zeros_like(counts)
	pixi = (np.floor(x - xmin)/deltax).astype(int)
	pixi[np.where(pixi == nbins)] -= 1
	pixj = (np.floor(y - ymin)/deltay).astype(int)
	pixj[np.where(pixj == nbins)] -= 1
	for i, j, pc_i in zip(pixi, pixj, pix_count) :
		counts[i,j] += pc_i
		area[i,j] += pix_area
	del pix_count, all_ra, all_dec, x, y

	# Now find the densities
	pix_density = np.empty_like(counts)
	for i in range(nbins) :
		for j in range(nbins) :
			try :
				pix_density[i,j] = counts[i,j]/area[i,j]
			except ZeroDivisionError :
				pix_density[i,j] = 0.0
	del counts, area
	pix_density = np.nan_to_num(pix_density)

	# Create arrays of the locations of the edges of pixels
	xedges = np.linspace(xmin, xmax, nedges)
	yedges = np.linspace(ymin, ymax, nedges)
	yedges = np.rad2deg(np.arcsin(yedges))
	xedges, yedges = np.meshgrid(xedges, yedges, indexing='ij')
	if vmax is None :
		vmax = pix_density.max()
	
	pc = axis.pcolormesh(xedges, yedges, pix_density.T, vmin=0, vmax=vmax, cmap='Reds')

	# Return objects
	return (pc, vmax)

def mask_point_density(axis, ra, dec, nbins, xmin=0., xmax=360., ymin=-1., ymax=1., vmax=None) :
	# Import modules
	# import time
	import numpy as np
	import matplotlib.pyplot as plt

	# Find widths of each self defined pixel
	deltax = (float(xmax) - float(xmin))/float(nbins)
	deltay = (float(ymax) - float(ymin))/float(nbins)
	pix_area = deltax*deltay

	H, xedges, yedges = np.histogram2d(ra, np.sin(np.deg2rad(dec)), bins=nbins, range=[[xmin, xmax], [ymin, ymax]])
	H = H/(pix_area*ra.size)
	# yedges = np.rad2deg(np.arcsin(yedges))
	xedges, yedges = np.meshgrid(xedges, yedges, indexing='ij')
	
	if vmax is None :
		vmax = H.max()
	
	pc = axis.pcolormesh(xedges, yedges, H, vmin=0, vmax=vmax, cmap='Reds')

	# Return objects
	return (pc)