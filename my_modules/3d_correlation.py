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
def run_TreeCorr_single_auto(pairing, config={}) :
    import treecorr
    import numpy as np
    
    # Create the catalogs
    cats_corr = catalog(pairing, perp=perp)
    
    # Set up NNCorrelation with the config dictionary
    nn = treecorr.NNCorrelation(**config)
    nn.process(cats_corr)
    # Save the results
    results = {'npairs': nn.npairs,
               'tot': nn.tot,
               'logr': nn.logr,
               'meanlogr': nn.meanlogr,
               'b': nn.b,
               'bin_size': nn.bin_size,
               'bin_slop': nn.bin_slop,
               'config': nn.config,
               'corr': nn.corr,
               'log_sep_units': nn.log_sep_units,
               'max_sep': nn.max_sep,
               'min_sep': nn.min_sep,
               'sep_units': nn.sep_units}
    # Return the results
    return (results)


#### TreeCorr single cross correlation paircount
def run_TreeCorr_single_cross(pairing1, pairing2, config={}) :
    import treecorr
    import numpy as np
    
    # Create the catalogs
    cats_corr1 = catalog(pairing1)
    cats_corr2 = catalog(pairing2)
    
    # Set up NNCorrelation with the config dictionary
    nn = treecorr.NNCorrelation(**config)
    nn.process(cats_corr1, cats_corr2)
    # Save the results
    results = {'npairs': nn.npairs,
               'tot': nn.tot,
               'logr': nn.logr,
               'meanlogr': nn.meanlogr,
               'b': nn.b,
               'bin_size': nn.bin_size,
               'bin_slop': nn.bin_slop,
               'config': nn.config,
               'corr': nn.corr,
               'log_sep_units': nn.log_sep_units,
               'max_sep': nn.max_sep,
               'min_sep': nn.min_sep,
               'sep_units': nn.sep_units}
    
    # Return the results
    return (results)


#### Run TreeCorr for all pairings for a set of data and randoms
def paircount_treecorr(data_ra, data_dec, data_r, rand_ra, rand_dec, rand_r, config={}, ret_time=False) :
    # Import modules
    import time
    import numpy as np
    import sys
    
    if ret_time :
		# Start timing
		start_time = time.time()
    
    # Store pairings and labels
	pairings = [[[data_ra.copy(), data_dec.copy(), data_r.copy()], [data_ra.copy(), data_dec.copy(), data_r.copy()]],
				[[data_ra.copy(), data_dec.copy(), data_r.copy()], [rand_ra.copy(), rand_dec.copy(), rand_r.copy()]],
				[[rand_ra.copy(), rand_dec.copy(), rand_r.copy()], [data_ra.copy(), data_dec.copy(), data_r.copy()]],
				[[rand_ra.copy(), rand_dec.copy(), rand_r.copy()], [rand_ra.copy(), rand_dec.copy(), rand_r.copy()]]]
    labels = ['dd', 'dr', 'rd', 'rr']
    
    # Call TreeCorr for all pairings
    results = {}
    for pairing, label in zip(pairings, labels) :
        if len(pairing)==1 :
            results[label] = run_TreeCorr_single_auto(pairing[0], config=config)
        elif len(pairing)==2 :
            results[label] = run_TreeCorr_single_cross(pairing[0], pairing[1])
            results[label]['npairs'] = results[label]['npairs']/2.0
            results[label]['tot'] = results[label]['tot']/2.0
    
    if ret_time :
        elapsed_time = time.time() - start_time
        del start_time

        # Return results
        return (results, elapsed_time)
    else :
        return (results)


#### Run TreeCorr over KMeans Jackknife regions
def paircount_treecorr_regions(data_ra, data_dec, data_r, data_cent, rand_ra, rand_dec, rand_r, rand_cent, ncen, nbins, config={}, ret_time=False) :
    # Import modules
    import time
    import numpy as np
    import sys
    
    results = {}
    labels = ['dd', 'dr', 'rd', 'rr']
    for label in labels :
        results[label] = {}
    
    # Set flags: True means that there is at least one galaxy in a region, False means the region is empty
    data_reg_flags = np.bincount(data_cent, minlength=ncen) > 0
    rand_reg_flags = np.bincount(rand_cent, minlength=ncen) > 0
    
    if ret_time :
		# Start timing
		start_time = time.time()
	
    for i in range(ncen) :
        ## Initialize the i'th region in results[labels]
        for label in labels :
            results[label][i] = {}
        data_i = (data_cent == i)
        rand_i = (rand_cent == i)
        drai = data_ra[data_i].copy()
        ddeci = data_dec[data_i].copy()
        rrai = rand_ra[rand_i].copy()
        rdeci = rand_dec[rand_i].copy()
		dri = data_r[data_i].copy()
		rri = rand_r[rand_i].copy()
        for j in range(i, ncen) :
            for label in labels :
                results[label][i][j] = {}
            pc_flags = [data_reg_flags[i] and data_reg_flags[j],
                        data_reg_flags[i] and rand_reg_flags[j],
                        rand_reg_flags[i] and data_reg_flags[j],
                        rand_reg_flags[i] and rand_reg_flags[j]]
            data_j = (data_cent == j)
            rand_j = (rand_cent == j)
            draj = data_ra[data_j].copy()
            ddecj = data_dec[data_j].copy()
            rraj = rand_ra[rand_j].copy()
            rdecj = rand_dec[rand_j].copy()
			drj = data_r[data_j].copy()
			rrj = rand_r[rand_j].copy()
			pairings = [[[drai.copy(), ddeci.copy(), dri.copy()], [draj.copy(), ddecj.copy(), drj.copy()]],
						[[drai.copy(), ddeci.copy(), dri.copy()], [rraj.copy(), rdecj.copy(), rrj.copy()]],
						[[rrai.copy(), rdeci.copy(), rri.copy()], [draj.copy(), ddecj.copy(), drj.copy()]],
						[[rrai.copy(), rdeci.copy(), rri.copy()], [rraj.copy(), rdecj.copy(), rrj.copy()]]]
            for pc_flag, pairing, label in zip(pc_flags, pairings, labels) :
                if pc_flag :
                    results[label][i][j] = run_TreeCorr_single_cross(pairing[0], pairing[1], config=config)
                    if i == j :
	                    results[label][i][j]['npairs'] = results[label][i][j]['npairs']/2.0
	                    results[label][i][j]['tot'] = results[label][i][j]['tot']/2.0
                else :
                    results[label][i][j]['npairs'] = np.zeros(nbins)
                    results[label][i][j]['tot'] = 0
            del pc_flags, data_j, rand_j, draj, ddecj, rraj, rdecj, pairings
            del drj, rrj
        # Delete items
        del data_i, rand_i, drai, ddeci, rrai, rdeci
        del dri, rri
    
    # Delete items
    del labels
    
    if ret_time :
        elapsed_time = time.time() - start_time
        del start_time

        # Return results
        return (results, elapsed_time)
    else :
        return (results)


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