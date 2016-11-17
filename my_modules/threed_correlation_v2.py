#! /usr/bin/env python

### This function converts chord distances to great circle distances or great circle angular distances. C is the chord distance, and R is the radius if desired. If R is None (default), the angle is returned in degrees.
def great_circle_dist(C, R=None) :
	import numpy as np
	
	## The formula for the angle is delta(sigma) = 2*arcsin(C/2) according to Wikipedia
	delta_sigma = 2.*np.arcsin(C/2.)
	
	if R is not None :
		return (R*delta_sigma)
	else :
		return (np.rad2deg(delta_sigma))

#### TreeCorr catalog function
def catalog(pairing) :
	import treecorr
	import numpy as np

	cats_corr = []
	ra, dec, r = pairing
	cat = treecorr.Catalog(ra=ra, dec=dec, r=r, ra_units='deg', dec_units='deg')

	cats_corr.append(cat)
	return (cats_corr)


#### TreeCorr single autocorrelation paircount
def run_TreeCorr_single_auto(ra, dec, r, config={}) :
    import treecorr
    import numpy as np
    
    # Create the catalogs
    cat1 = treecorr.Catalog(ra=ra, dec=dec, r=r, ra_units='deg', dec_units='deg')
    
    # Set up NNCorrelation with the config dictionary
    nn = treecorr.NNCorrelation(**config)
    nn.process(cat1)
    
    # Return the results
    return (nn, cat1)


#### TreeCorr single cross correlation paircount
def run_TreeCorr_single_cross(ra1, dec1, r1, ra2, dec2, r2, config={}) :
    import treecorr
    import numpy as np
    
    # Create the catalogs
    cat1 = treecorr.Catalog(ra=ra1, dec=dec1, r=r1, ra_units='deg', dec_units='deg')
    cat2 = treecorr.Catalog(ra=ra2, dec=dec2, r=r2, ra_units='deg', dec_units='deg')
    
    # Set up NNCorrelation with the config dictionary
    nn = treecorr.NNCorrelation(**config)
    nn.process(cat1, cat2)
    
    # Return the results
    return (nn)


#### Run TreeCorr for all pairings for a set of data and randoms
def paircount_treecorr(data_ra, data_dec, data_r, rand_ra, rand_dec, rand_r, file_name, config={}, ret_time=False) :
	# Import modules
	import treecorr
	import time
	import numpy as np
	import sys

	if ret_time :  
		# Start timing
		start_time = time.time()

	# Call TreeCorr for dd, dr, and rr
	dd, dcat = run_TreeCorr_single_auto(data_ra, data_dec, data_r, config=config)
	nd = np.full_like(dd.npairs, data_ra.size)
	rr, rcat = run_TreeCorr_single_auto(rand_ra, rand_dec, rand_r, config=config)
	nr = np.full_like(rr.npairs, rand_ra.size)
	dr = run_TreeCorr_single_cross(data_ra, data_dec, data_r, rand_ra, rand_dec, rand_r, config=config)
	rd = run_TreeCorr_single_cross(rand_ra, rand_dec, rand_r, data_ra, data_dec, data_r, config=config)
	np.savetxt(file_name, np.array([dd.rnom, dd.weight, rr.weight, dr.weight, rd.weight, nd, nr]).T, fmt='%-25.18e', \
	header='{:<25s} {:<25s} {:<25s} {:<25s} {:<25s} {:<25s} {:<25s}'.format('# r (Mpc)', 'DD', 'RR', 'DR', 'RD', 'D', 'R'), \
	comments='')
	print 'Output file: {}'.format(file_name)
	sys.stdout.flush()

	if ret_time :
		elapsed_time = time.time() - start_time
		del start_time

		# Return results
		return (elapsed_time)
	else :
		return (dcat, rcat)


#### Run TreeCorr over KMeans Jackknife regions
def paircount_treecorr_regions(data_ra, data_dec, data_r, data_cent, rand_ra, rand_dec, rand_r, rand_cent, ncen, nbins, file_name, file_ext, config={}, ret_time=False) :
	# Import modules
	import treecorr
	import time
	import numpy as np
	import sys

	# Set flags: True means that there is at least one galaxy in a region, False means the region is empty
	data_reg_flags = np.bincount(data_cent, minlength=ncen) > 0
	rand_reg_flags = np.bincount(rand_cent, minlength=ncen) > 0

	if ret_time :  
		# Start timing
		start_time = time.time()

	for i in range(ncen) :  
		data_i = (data_cent == i)
		rand_i = (rand_cent == i)
		drai = data_ra[data_i].copy()
		ddeci = data_dec[data_i].copy()
		dri = data_r[data_i].copy()
		rrai = rand_ra[rand_i].copy()
		rdeci = rand_dec[rand_i].copy()
		rri = rand_r[rand_i].copy()
		ndi = np.full(nbins, drai.size)
		nri = np.full_like(ndi, rrai.size)
	#         dd = run_TreeCorr_single_auto(drai, ddeci, dri, config=config)
	#         rr = run_TreeCorr_single_auto(rrai, rdeci, rri, config=config)
	#         rr.weight = rr.weight*(rr.tot/dd.tot)
	#         dr = run_TreeCorr_single_cross(drai, ddeci, dri, rrai, rdeci, rri, config=config)
	#         dr.weight = dr.weight*(dr.tot/dd.tot)
	#         rd = run_TreeCorr_single_cross(rrai, rdeci, rri, drai, ddeci, dri, config=config)
	#         rd.weight = rd.weight*(rd.tot/dd.tot)
	#         dd.write(file_name+'_bin{}_bin{}'.format(i, i)+file_ext, rr, dr, rd)
	#         print 'Output file: {}'.format(file_name+'_bin{}_bin{}'.format(i, i)+file_ext)
	#         sys.stdout.flush()
	#         del dd, dr, rd, rr
		for j in range(i, ncen) :
			data_j = (data_cent == j)
			rand_j = (rand_cent == j)
			draj = data_ra[data_j].copy()
			ddecj = data_dec[data_j].copy()
			drj = data_r[data_j].copy()
			rraj = rand_ra[rand_j].copy()
			rdecj = rand_dec[rand_j].copy()
			rrj = rand_r[rand_j].copy()
			if i == j :
				dd = run_TreeCorr_single_auto(drai, ddeci, dri, config=config)
				rr = run_TreeCorr_single_auto(rrai, rdeci, rri, config=config)
			else :
				dd = run_TreeCorr_single_cross(drai, ddeci, dri, draj, ddecj, drj, config=config)
				rr = run_TreeCorr_single_cross(rrai, rdeci, rri, rraj, rdecj, rrj, config=config)
			ndj = np.full_like(ndi, draj.size)
			nrj = np.full_like(ndi, rraj.size)
			dr = run_TreeCorr_single_cross(drai, ddeci, dri, rraj, rdecj, rrj, config=config)
			rd = run_TreeCorr_single_cross(rrai, rdeci, rri, draj, ddecj, drj, config=config)
			np.savetxt(file_name+'_bin{}_bin{}'.format(i, j)+file_ext, np.array([dd.rnom, dd.weight, rr.weight, dr.weight, rd.weight, ndi, ndj, nri, nrj]).T, \
			fmt='%-25.18e', header='{:<25s} {:<25s} {:<25s} {:<25s} {:<25s} {:<25s} {:<25s} {:<25s} {:<25s}'.format('# r (Mpc)', 'DD', 'RR', 'DR', 'RD', 'D_i', 'D_j', 'R_i', 'R_j'), \
			comments='')
			print 'Output file: {}'.format(file_name+'_bin{}_bin{}'.format(i, j)+file_ext)
			sys.stdout.flush()
			del data_j, rand_j, draj, ddecj, drj, rraj, rdecj, rrj
			del dd, dr, rd, rr
			del ndj, nrj
		# Delete items
		del data_i, rand_i, drai, ddeci, dri, rrai, rdeci, rri. ndi, nri

	if ret_time :
		elapsed_time = time.time() - start_time
		del start_time

		# Return results
		return (elapsed_time)
	else :
		return ()

# This function finds the CF given the pair counts as inputs. These can be the average
# pair counts or pair counts from a single run for the CFs to be averaged later
def pc_corr(dd, rr, dr, rd, nd, nr) :
	import numpy as np
	# w = (dd - dr - rd + rr)/(1.*rr)
	w = ((2.*dd/(nd*(nd - 1.))) - (dr/(nd*nr)) - (rd/(nd*nr)) + (2.*rr/(nr*(nr - 1.))))/(2.*rr/(nr*(nr - 1.)))
	return w
	
# This function finds the average pair counts over a number of runs saved into files. Can
# not handle jackknifes yet! Note that file_name should not include the extension!
## Inputs :
## file_name = path without extension to pair count files. If multiple regions, this should
## not include the _run[i] portion of the file name, as this will be looped over
## file_ext = extension for pair count files (including period)
## niters = # of runs that were used
## Returns average pair counts dd, rr, dr (and rd if in file)
def ave_pc(file_name, file_ext, niters) :
	import numpy as np

	r = np.loadtxt(file_name+'_run0'+file_ext, usecols=(0,))
	dd = np.empty((niters, r.size), dtype=np.float64)
	rr = np.empty_like(dd)
	dr = np.empty_like(dd)
	rd = np.empty_like(dd)
	for i in range(niters) :
		dd[i], rr[i], dr[i], rd[i], nd, nr = read_pcs(file_name+'_run{}'.format(i), file_ext)
	dd_ave = dd.mean(axis=0)
	rr_ave = rr.mean(axis=0)
	dr_ave = dr.mean(axis=0)
	rd_ave = rd.mean(axis=0)
	if niters > 1 :
		ddi_ave = np.empty_like(dd)
		rri_ave = np.empty_like(rr)
		dri_ave = np.empty_like(dr)
		rdi_ave = np.empty_like(rd)
		for i in range(niters) :
			ddi_ave[i] = (1./(float(niters) - 1.))*(np.sum(dd, axis=0) - dd[i])
			rri_ave[i] = (1./(float(niters) - 1.))*(np.sum(rr, axis=0) - rr[i])
			dri_ave[i] = (1./(float(niters) - 1.))*(np.sum(dr, axis=0) - dr[i])
			rdi_ave[i] = (1./(float(niters) - 1.))*(np.sum(rd, axis=0) - rd[i])
		return r, dd_ave, rr_ave, dr_ave, rd_ave, nd, nr, ddi_ave, rri_ave, dri_ave, rdi_ave
	else :
		return r, dd_ave, rr_ave, dr_ave, rd_ave, nd, nr

# This function handles the calculation of the CF either from average pair counts or
# by averaging CFs
def ave_cf(file_name, file_ext, niters, flag) :
	import numpy as np
	
	if flag == 'cf' :
		r = np.loadtxt(file_name+'_run0'+file_ext, usecols=(0,))
		w = np.empty((niters, r.size), dtype=np.float64)
		for i in range(niters) :
			dd, rr, dr, rd, nd, nr = read_pcs(file_name+'_run{}'.format(i), file_ext)
			w[i,:] = pc_corr(dd.copy(), rr.copy(), dr.copy(), rd.copy(), nd, nr)
			del dd, rr, dr, rd, nd, nr
		wave = w.mean(axis=0)
		if niters > 1 :
			w_cov = (1./(float(niters) - 1.))*np.dot((w - wave).T, (w - wave))
			werr_ave = np.sqrt(np.diag(w_cov))
	
	elif flag == 'pc' :
		if niters > 1 :
			r, dd, rr, dr, rd, nd, nr, ddi, rri, dri, rdi = ave_pc(file_name, file_ext, niters)
			wave_i = np.empty((niters, r.size))
			for i in range(niters) :
				wave_i[i] = pc_corr(ddi[i], rri[i], dri[i], rdi[i], nd, nr)
			wave_iave = (1./float(niters))*np.sum(wave_i, axis=0)
			w_cov = ((float(niters) - 1.)/float(niters))*np.dot((wave_i - wave_iave).T, (wave_i - wave_iave))
			werr_ave = np.sqrt(np.diag(w_cov))
		else :
			r, dd, rr, dr, rd, nd, nr = ave_pc(file_name, file_ext, niters)
		wave = pc_corr(dd, rr, dr, rd, nd, nr)
	
	if niters > 1 :
		return r, wave, werr_ave, w_cov
	else :
		return r, wave

def read_pcs(file_name, file_ext) :
	import numpy as np
	r, dd, rr, dr, rd, nd, nr = np.loadtxt(file_name+file_ext, unpack=True)
	# dd *= (2./(nd*(nd - 1.)))
	# rr *= (2./(nr*(nr - 1.)))
	# dr *= (1./(nd*nr))
	# rd *= (1./(nd*nr))
	return dd, rr, dr, rd, nd[0], nr[0]

#### Function to calculate correlation function on full sample
def tc_ang_corr_full(results, ret_time=False) :
    # Import python modules
    import time
    import numpy as np
    import sys
    
    if ret_time : 
    	# Start timing
        start_time = time.time()
    
    ## Correlation function using full sample
    # DD
    npairs_dd = results['dd']['npairs'].copy()
    tot_dd = results['dd']['tot']
    nd = np.sqrt(2.*tot_dd)
    norm_dd = (nd*(nd - 1.))/2.
    logr = results['dd']['meanlogr'].copy()
    r = np.exp(logr.copy())
    dd = (1.0*npairs_dd.copy())/norm_dd
    ## delete logr
    del logr
    results['dd']['normed'] = dd.copy()
    
    # RR
    npairs_rr = results['rr']['npairs'].copy()
    tot_rr = results['rr']['tot']
    nr = np.sqrt(2.*tot_rr)
    norm_rr = (nr*(nr - 1.))/2.
    rr = (1.0*npairs_rr.copy())/norm_rr
    results['rr']['normed'] = rr.copy()
    
    # DR
    npairs_dr = results['dr']['npairs'].copy()
    norm_dr = (nd*nr)/2.
    dr = (1.0*npairs_dr.copy())/norm_dr
    results['dr']['normed'] = dr.copy()
    
    # RD
    npairs_rd = results['rd']['npairs'].copy()
    norm_rd = (nd*nr)/2.
    rd = (1.0*npairs_rd.copy())/norm_rd
    results['rd']['normed'] = rd.copy()
    
    # Calculate w
    # w = (dd.copy() - 2*dr.copy() + rr.copy())/rr.copy()
    w = (dd.copy()/rr.copy()) - (dr.copy()/rr.copy()) - (rd.copy()/rr.copy()) + 1.0
    results['w'] = w.copy()
    results['r'] = r.copy()
    del npairs_dd, npairs_dr, npairs_rr
    del norm_dd, norm_dr, norm_rr
    
    if ret_time :
        elapsed_time = time.time() - start_time
        del start_time

        # Return results
        return (results, r, w, elapsed_time)
    else :
        return (results, r, w)


#### Calculate correlation function using jackknife regions
def tc_ang_corr_reg(nbins, ncen, results, ret_time=False) :
    # Import python modules
    import time
    import numpy as np
    import sys
    
    if ret_time : 
    	# Start timing
        start_time = time.time()
    
    ## Correlation function using all regions
    ### Find first non-empty logr to read r from
    class Found(Exception) : pass
    try :
    	for i in range(ncen) :
    		for j in range(i, ncen) :
    			if 'logr' in results['dd'][i][j] :
    				raise Found
    except Found :
    	r = np.exp(results['dd'][i][j]['logr'].copy())
    results['r'] = r.copy()
    npairs_dd = np.zeros(nbins)
    tot_dd = 0
    npairs_dr = np.zeros(nbins)
    tot_dr = 0
    npairs_rd = np.zeros(nbins)
    tot_rd = 0
    npairs_rr = np.zeros(nbins)
    tot_rr = 0
    results['dd']['summed'] = {}
    results['dd']['summed']['full'] = {}
    results['dr']['summed'] = {}
    results['dr']['summed']['full'] = {}
    results['rd']['summed'] = {}
    results['rd']['summed']['full'] = {}
    results['rr']['summed'] = {}
    results['rr']['summed']['full'] = {}
    results['w'] = {}
    for i in range(ncen) :
        for j in range(i, ncen) :
            npairs_dd += results['dd'][i][j]['npairs'].copy()
            tot_dd += results['dd'][i][j]['tot']
            npairs_dr += results['dr'][i][j]['npairs'].copy()
            tot_dr += results['dr'][i][j]['tot']
            npairs_rd += results['rd'][i][j]['npairs'].copy()
            tot_rd += results['rd'][i][j]['tot']
            npairs_rr += results['rr'][i][j]['npairs'].copy()
            tot_rr += results['rr'][i][j]['tot']
    nd = np.sqrt(2.*tot_dd)
    nr = np.sqrt(2.*tot_rr)
    norm_dd = (nd*(nd - 1.))/2.
    norm_dr = (nd*nr)/2.
    norm_rd = (nd*nr)/2.
    norm_rr = (nr*(nr - 1.))/2.
    dd = (1.0*npairs_dd.copy())/norm_dd
    dr = (1.0*npairs_dr.copy())/norm_dr
    rd = (1.0*npairs_rd.copy())/norm_rd
    rr = (1.0*npairs_rr.copy())/norm_rr
    results['dd']['summed']['full']['npairs'] = 1.0*npairs_dd.copy()
    results['dd']['summed']['full']['tot'] = nd
    results['dd']['summed']['full']['norm'] = norm_dd
    results['dr']['summed']['full']['npairs'] = 1.0*npairs_dr.copy()
    results['dr']['summed']['full']['norm'] = norm_dr
    results['rd']['summed']['full']['npairs'] = 1.0*npairs_rd.copy()
    results['rd']['summed']['full']['norm'] = norm_rd
    results['rr']['summed']['full']['npairs'] = 1.0*npairs_rr.copy()
    results['rr']['summed']['full']['tot'] = nr
    results['rr']['summed']['full']['norm'] = norm_rr
    w_reg = (dd.copy()/rr.copy()) - (dr.copy()/rr.copy()) - (rd.copy()/rr.copy()) + 1.0
    results['w']['all'] = w_reg.copy()
    del npairs_dd, npairs_dr, npairs_rd, npairs_rr
    del tot_dd, tot_dr, tot_rd, tot_rr
    del norm_dd, norm_dr, norm_rd, norm_rr
    del nd, nr
    del dd, dr, rd, rr
    
    ## Now find the errors
    w_i = np.empty((ncen, nbins))
    for i in range(ncen) :
    	results['dd']['summed'][i] = {}
    	results['dr']['summed'][i] = {}
    	results['rd']['summed'][i] = {}
    	results['rr']['summed'][i] = {}
        npairs_dd = np.zeros_like(w_reg)
        tot_dd = 0
        npairs_dr = np.zeros_like(w_reg)
        tot_dr = 0
        npairs_rd = np.zeros_like(w_reg)
        tot_rd = 0
        npairs_rr = np.zeros_like(w_reg)
        tot_rr = 0
        for j in range(ncen) :
            if not j == i :
                for k in range(j, ncen) :
                    if not k == i :
                        npairs_dd += results['dd'][j][k]['npairs'].copy()
                        tot_dd += results['dd'][j][k]['tot']
                        npairs_dr += results['dr'][j][k]['npairs'].copy()
                        tot_dr += results['dr'][j][k]['tot']
                        npairs_rd += results['rd'][j][k]['npairs'].copy()
                        tot_rd += results['rd'][j][k]['tot']
                        npairs_rr += results['rr'][j][k]['npairs'].copy()
                        tot_rr += results['rr'][j][k]['tot']
        nd = np.sqrt(2.*tot_dd)
        nr = np.sqrt(2.*tot_rr)
        norm_dd = (nd*(nd - 1.))/2.
        norm_dr = (nd*nr)/2.
        norm_rd = (nd*nr)/2.
        norm_rr = (nr*(nr - 1.))/2.
        dd = (1.0*npairs_dd.copy())/norm_dd
        dr = (1.0*npairs_dr.copy())/norm_dr
        rd = (1.0*npairs_rd.copy())/norm_rd
        rr = (1.0*npairs_rr.copy())/norm_rr
        results['dd']['summed'][i]['npairs'] = 1.0*npairs_dd.copy()
        results['dd']['summed'][i]['tot'] = nd
        results['dd']['summed'][i]['norm'] = norm_dd
        results['dr']['summed'][i]['npairs'] = 1.0*npairs_dr.copy()
        results['dr']['summed'][i]['norm'] = norm_dr
        results['rd']['summed'][i]['npairs'] = 1.0*npairs_rd.copy()
        results['rd']['summed'][i]['norm'] = norm_rd
        results['rr']['summed'][i]['npairs'] = 1.0*npairs_rr.copy()
        results['rr']['summed'][i]['tot'] = nr
        results['rr']['summed'][i]['norm'] = norm_rr
        ## delete 
        del npairs_dd, npairs_dr, npairs_rd, npairs_rr
        del tot_dd, tot_dr, tot_rd, tot_rr
        w_i[i] = (dd.copy()/rr.copy()) - (dr.copy()/rr.copy()) - (rd.copy()/rr.copy()) + 1.0
        results['w'][i] = w_i[i].copy()
        del dd, dr, rd, rr
        del norm_dd, norm_dr, norm_rd, norm_rr
        del nd, nr
    w_ave = (1./float(ncen))*np.sum(w_i, axis=0)
    sig_sq = (float(ncen-1)/float(ncen))*np.dot((w_i.copy() - w_ave.copy()).T, (w_i.copy() - w_ave.copy()))
    results['cov'] = sig_sq.copy()
    
    if ret_time :
        elapsed_time = time.time() - start_time
        del start_time

        # Return results
        return (results, r, w_reg, sig_sq, elapsed_time)
    else :
        return (results, r, w_reg, sig_sq)


#### TreeCorr Xi calculation
def correlation_TreeCorr(data_ra, data_dec, data_r, rand_ra, rand_dec, rand_r, config) :
	import time
	import numpy as np
	import treecorr
	import sys

	# Begin timing
	start = time.time()

	# Make sure arrays match
	assert data_ra.size == data_dec.size, "Data must have both RA and DEC"
	assert rand_ra.size == rand_dec.size, "Randoms must have both RA and DEC"

	# Create TreeCorr catalog objects
	dcat = treecorr.Catalog(ra=data_ra, dec=data_dec, r=data_r, ra_units='deg', dec_units='deg')
	rcat = treecorr.Catalog(ra=rand_ra, dec=rand_dec, r=rand_r, ra_units='deg', dec_units='deg')
	print ('TreeCorr catalogs created')
	sys.stdout.flush()

	# Run TreeCorr processes for DD, DR, RD, and RR
	dd = treecorr.NNCorrelation(config)
	dr = treecorr.NNCorrelation(config)
	# rd = treecorr.NNCorrelation(config)
	rr = treecorr.NNCorrelation(config)
	dd.process(dcat)
	print ('DD done')
	sys.stdout.flush()
	# I also need to get the bin locations for plotting
	logr = dd.logr
	dr.process(dcat, rcat)
	print ('DR done')
	sys.stdout.flush()
	# rd.process(rcat, dcat)
	# print ('RD done')
	# sys.stdout.flush()
	rr.process(rcat)
	print ('RR done')
	sys.stdout.flush()

	# Find the correlation function and errors
	# xi, varxi = dd.calculateXi(rr, dr, rd)
	xi, varxi = dd.calculateXi(rr, dr)
	print ('Correlation function and errors calculated')
	sys.stdout.flush()

	# Find elapsed time
	runtime = time.time() - start
	del start
	## Print the time it took
	h = int(np.floor(runtime/(60.0*60.0)))
	m = int(np.floor((runtime - (60.0*60.0*h))/60.0))
	s = runtime - 60.0*60.0*h - 60.0*m
	print ('Elapsed time: {:>02d}:{:>02d}:{:>05.2f}'.format(h, m, s))
	sys.stdout.flush()
	del runtime, h, m, s

	# Return xi, varxi, and bin locations
	return (xi, varxi, logr)