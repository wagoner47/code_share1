
# coding: utf-8

# # Module jackknife.py
# ## Version 2.3
# ## Last Update: August 25, 2015
# ### [Erika Wagoner](mailto:wagoner47@email.arizona.edu "send an email")
# Module defining functions relating to jackknife region creation and use.
# 
# Available functions:
# 
# - find_kmeans_centers: Run KMeans to get the centers for a set of data
# - find_kmeans_labels: Get the labels for points in a set of data given centers
# - equal_areas_centers_test: Function that runs many iterations of the method for equal_areas_centers and checks the $\chi^2$ value for each iteration with the assumed uniformity to see if it is possible to get equal areas KMeans regions.
# - equal_areas_centers: Method that attempts to use an initial guess that will produce equal area KMeans regions. Not currently available.
# - region_map: Function that plots the KMeans regions color coded with the centers as dots given the points, labels, and centers, as well as a few other parameters. Returns the figure object to be saved if desired.
# 
# In version 2.0, I've allowed for the possibility of using either north and south data or having the code split the data into north and south. I've also added a check to make sure the number of regions isn't too great to be able to have the same area to within about 5%. Worked out, this means that we need to have a maximum of 400 points/region, so there is a check to make sure that npoints/ncen $\le$ 400. This failing, the code defaults to the maximum number of centers allowable such that there are at least 400 points/region. But if there are more points than needed to have 400 points/region, the code also down samples to only use 400$\cdot$ncen points when finding the centers with KMeans.
# 
# ##### Update on August 13, 2015
# Instead of choosing a new subsample for each iteration in equal_areas_centers_test, I am trying setting the down sample once and using that to find the centers every time. This way, the points used when finding the centers are consistent, and not causing additional differences in the chosen centers.
# 
# ##### Update on August 17, 2015
# I've added an additional function (see above) called region_map to plot the KMeans regions for points. This doesn't account for splitting the data between RA and DEC, so I'll need to append the two for it to work. I'll also need to add the ability to get out the actual ncen used if changed and the locations of the centers.
# 
# ##### Update on August 18, 2015
# I think it is more important to check the ability to handle different numbers of regions first, and then check each ncen for different numbers of iterations. This is because I'm not entirely convinced all is working, even when using the subsample, when I begin increasing the number of regions. Thus, I'm modifying what I output to be easier to use the code in a for loop. I will now output the average time/iteration (rather than an array of times), the $\chi^2$ values from just the subsample, the subsample used, and the centers and labels for the iteration with the lowest $\chi^2$. I want to be able to plot the regions for NGC and SGC separately, as they are too difficult to see in the same plot, so I will either need to output the labels for NGC and SGC separately, or else I will need to output the number of regions used for NGC to allow the labels to be separated by the user later. I have not yet decided which of these is the better idea, but I feel like outputting the labels and centers separately might be easier to work with.
# 
# ##### Update on August 25, 2015
# I need to allow for the case that I only have NGC or SGC when running KMeans, as I am to test them separately for now. Thus, I added some flags to check if I have NGC points, SGC points, or both.

# #### find_kmeans_centers
# Inputs: 
# - ncen = # of regions
# - maxiter = maximum # of KMeans iterations
# - tol = KMeans tolerance limit
# - points = array of shape [2, nsamp] or [nsamp, 2] containing the RA and DEC of the sample points
# - centers (optional, default=None) = array of shape [2, ncen] or [ncen, 2] containing the RA and DEC of the guessed centers
# - verbose (optional, default=1) = verbose option for KMeans
# 
# Process: If centers is given, use kmeans_radec.kmeans with the sample (points) and guessed centers (centers), also feeding the tolerance (tol), maximum iterations (maxiter), and verbose options. Otherwise, use kmeans_radec.kmeans_sample on the points with specified ncen and the other options. The final center locations can be obtained using the centers attribute of the KMeans object after running in either case.
# 
# Output: Array of shape [ncen, 2] containing the RA and DEC of the centers from KMeans.

# In[ ]:

def find_kmeans_centers(ncen, maxiter, tol, points, centers=None, verbose=1) :
    import numpy as np
    import kmeans_radec as kmrd
    import sys
    
    np.random.seed(0)
    # Make sure pra and pdec are both given
    assert points.shape[0] == 2 or points.shape[1] == 2, "Must have RA and Dec"
    if points.shape[0] == 2 :
        points.T
    
    # Check if centers are guessed, and if so, make sure both ra and dec are given
    if centers is not None :
        assert centers.shape[0] == 2 or centers.shape[1] == 2, "Centers must have RA and Dec if given"
        if centers.shape[0] == 2 :
            centers.T
        km = kmrd.kmeans(points, centers, tol=tol, maxiter=maxiter, verbose=verbose)
    else :
        km = kmrd.kmeans_sample(points, ncen, maxiter=maxiter, tol=tol, verbose=verbose)
    
    return (km.centers)


# #### find_kmeans_labels
# Inputs: 
# - points = array of shape [npoints, 2] or [2, npoints] containing the RA and DEC of the points.
# - centers = array of shape [ncen, 2] containing the RA and DEC of the centers.
# 
# Process: Use kmeans_radec.find_nearest to find the labels of all points given the positions of the centers.
# 
# Output: Array of labels giving the region for each point in points.

# In[ ]:

def find_kmeans_labels(points, centers) :
    import numpy as np
    import kmeans_radec as kmrd
    import sys
    
    # Make sure ra, dec, and centers have correct shapes
    assert points.shape[0] == 2 or points.shape[1] == 2, "Points must have RA and Dec"
    if points.shape[0] == 2 :
        points.T
    assert centers.shape[0] == 2 or centers.shape[1] == 2, "Centers must have RA and Dec"
    if centers.shape[0] == 2 :
        centers.T
    
    # Run find_nearest
    return (kmrd.find_nearest(points, centers))


# #### equal_areas_centers_test
# Inputs: 
# - coords = array of shape [npoints, 2] or [2, npoints] containing the RA and DEC of points.
# - ncen = # of regions to use. If npoints/ncen $\le$ 400, default to ncen = floor(npoints/400).
# - maxiter = maximum # of iterations for KMeans.
# - tol = tolerance limit for KMeans.
# - niters (optional, default=100) = # of times to guess new centers and run to see if it is possible to reach equal area regions within the 5% tolerance limit.
# - verbose (optional, default=1) = the verbose option for KMeans.
# - scoords (optional, default=None) = array of shape [npoints_sgc, 2] or [2, npoints_sgc] containing RA and DEC of SGC points. If given, ncen is split between north and south according to the fraction of the total number of points contained in each, and the tolerance limit must be met in both regions. If not given, code will try to split coords into north and south and assign ncen to each of those as appropriate.
# 
# Process: First split ncen between north and south based on the fraction of the total points in each. If the data is truly only in north or south, the other will not have any regions and will not be run through KMeans. Next, the code tests to make sure the number of centers given allows for at least 400 points per region in both north and south. If not, the number of centers defaults to the maximum number of centers such that this is the case.
# 
# If the number of points in north/south is greater than 400 times the number of centers in the same area, the code will randomly down sample in this area when calling the function that finds the centers. But the guess for the centers will be drawn as ncen_ngc/ncen_sgc random points within north/south before being down sampled. These random points and the sample (with enough for 400 points/region) will be fed to the function find_kmeans_centers to get the guessed centers. Then, the same sample will be given with the final centers to find_kmeans_labels to get the regions for each point in the sample.
# 
# The $\chi^2$ will be computed with a model of constant 400 points/region on just the subsample, or the full sample if it works out to be exactly 400 points/region.
# 
# Output:
# - Average time/iteration
# - subsample used, with NGC and SGC separate and SGC shifted
# - $\chi^2$ for just the down sample
# - centers for the iteration with the lowest $\chi^2$, with NGC and SGC separate
# - labels for the iteration with the lowest $\chi^2$, with NGC and SGC separate

# In[ ]:

def equal_areas_centers_test(coords, ncen, maxiter, tol, niters=100, verbose=1, scoords=None) :
    import numpy as np
    import sys
    import time
    
    nexp = 400
    # Make sure things were input correctly
    assert coords.shape[0] == 2 or coords.shape[1] == 2, "Coordinates must have RA and Dec"
    if coords.shape[0] == 2 :
        coords.T
    if scoords is not None :
        assert scoords.shape[0] == 2 or scoords.shape[1] == 2, "Coordinates must have RA and Dec"
        if scoords.shape[0] == 2 :
            scoords.T
        npoints_ngc = coords.shape[0]
        npoints_sgc = scoords.shape[0]
        ncoords = coords.copy()
    else :
        ngc_filt = (coords[:,0] > 60) & (coords[:,0] < 300)
        sgc_filt = (coords[:,0] < 60) | (coords[:,0] > 300)
        ncoords = coords[ngc_filt,:].copy()
        scoords = coords[sgc_filt,:].copy()
        npoints_ngc = ncoords.shape[0]
        npoints_sgc = scoords.shape[0]
    scoords[np.where(scoords[:,0] > 180),0] -= 360
    npoints = npoints_ngc + npoints_sgc
    # As long as the minimum points/region holds for the full set and the number of regions is split correctly
    # between NGC and SGC, the points/region will hold for NGC and SGC.
    if float(npoints)/float(ncen) < nexp :
        print ('Not enough points to support {} regions.         Defaulting to maximum number of regions for accuracy.'.format(ncen))
        sys.stdout.flush()
        ncen = int(np.floor(npoints/float(nexp)))
        print ('Number of regions set to {}'.format(ncen))
        sys.stdout.flush()
    if ncen == 2 and npoints_ngc > 0 and npoints_sgc > 0 :
        ncen_ngc = 1
        ncen_sgc = 1
    else :
        frac_ngc = float(npoints_ngc)/float(npoints)
        ncen_ngc = int(np.rint(frac_ngc*ncen))
        ncen_sgc = ncen - ncen_ngc
    sample_size_ngc = nexp*ncen_ngc
    sample_size_sgc = nexp*ncen_sgc
    if npoints_ngc > 0 and not sample_size_ngc > 0 :
        print ('NGC points exist but NGC sample size is 0')
        sys.stdout.flush()
    if npoints_sgc > 0 and not sample_size_sgc > 0 :
        print ('SGC points exist but SGC sample size is 0')
        sys.stdout.flush()
    
    np.random.seed(0)
    if sample_size_ngc > 0 :
        sample_ngc = np.random.choice(npoints_ngc, sample_size_ngc, replace=False)
    if sample_size_sgc > 0 :
        sample_sgc = np.random.choice(npoints_sgc, sample_size_sgc, replace=False)
    
    # Create a flag to more easily determine which regions to run KMeans on and output results for
    ngc_go = ncen_ngc > 0
    sgc_go = ncen_sgc > 0
    both_go = ngc_go and sgc_go
    
    # Run
    ave_time = 0.0
    chisq = np.empty(niters)
    if ngc_go :
        centers_ngc_best = np.empty((ncen_ngc,2))
        labels_ngc_best = np.empty(ncen_ngc)
    if sgc_go :
        centers_sgc_best = np.empty((ncen_sgc,2))
        labels_sgc_best = np.empty(ncen_sgc)
    for i in range(niters) :
        if i%(niters//10) == 0 :
            print i
            sys.stdout.flush()
        start = time.time()
        np.random.seed(0)
        # Because sub samples are taken with seed 0, start the seed at 1 here.
        np.random.seed(i+1)
        if ngc_go or both_go :
            choice_ngc = np.random.choice(npoints_ngc, ncen_ngc, replace=False)
            centers_ngc = find_kmeans_centers(ncen_ngc, maxiter, tol, points=ncoords[sample_ngc,:], 
                                              centers=ncoords[choice_ngc,:], verbose=verbose)
            labels_ngc = find_kmeans_labels(ncoords[sample_ngc,:], centers_ngc)
        if sgc_go or both_go :
            choice_sgc = np.random.choice(npoints_sgc, ncen_sgc, replace=False)
            centers_sgc = find_kmeans_centers(ncen_sgc, maxiter, tol, points=scoords[sample_sgc,:], 
                                              centers=scoords[choice_sgc,:], verbose=verbose)
            labels_sgc = find_kmeans_labels(scoords[sample_sgc,:], centers_sgc)
        if both_go :
            ni = np.append(np.bincount(labels_ngc,minlength=ncen_ngc),np.bincount(labels_sgc,minlength=ncen_sgc))
        else :
            if ngc_go :
                ni = np.bincount(labels_ngc, minlength=ncen_ngc)
            if sgc_go :
                ni = np.bincount(labels_sgc, minlength=ncen_sgc)
        chisq[i] = (((ni.astype(float) - float(nexp))**2)/float(nexp)).sum()
        if i == 0 or chisq[:i+1].argmin() == i :
            if ngc_go :
                centers_ngc_best = centers_ngc.copy()
                labels_ngc_best = labels_ngc.copy()
            if sgc_go :
                centers_sgc_best = centers_sgc.copy()
                labels_sgc_best = labels_sgc.copy()
        if sgc_go :
            del choice_sgc, centers_sgc, labels_sgc
        if ngc_go :
            del choice_ngc, centers_ngc, labels_ngc
        del ni
        ave_time += time.time() - start
        del start
    #return (chisq_ds, chisq, times, nexp_ngc, nexp_sgc, ncen, results)
    #scoords[np.where(scoords[:,0] < 0),0] += 360
    ave_time /= float(niters)
    if both_go :
        return (ave_time,ncoords[sample_ngc,:].T,scoords[sample_sgc,:].T,chisq,centers_ngc_best.T,centers_sgc_best.T,
                labels_ngc_best,labels_sgc_best)
    else :
        if ngc_go :
            return (ave_time,ncoords[sample_ngc,:].T,chisq,centers_ngc_best.T,labels_ngc_best)
        if sgc_go :
            return (ave_time,scoords[sample_sgc,:].T,chisq,centers_sgc_best.T,labels_sgc_best)


# #### region_map
# I need to make sure that the $\chi^2$ values are being calculated correctly, so I'm going to add an optional function to plot the points colored by region with the centers on top.
# 
# Inputs:
# - data_ra, data_dec: arrays of the RA and DEC for the points
# - labels: array of the labels of the regions for the points
# - center_ra, center_dec: arrays of the RA and DEC for the centers
# - area: string giving the area of the sky --- NGC, SGC, or full sky
# - niters: # of iterations run, used in title
# - xlim (optional, default=[-180,180]): specify the x limits for the plot
# - ylim (optional, defualt=[-90,90]): specify the y limits for the plot
# - ds (optional bool, default=True): specify if the data is only the downsample or the full catalog
# - ncen (optional, default=None): # of regions for the points, used in title, for setting up colors and markers, and for plotting regions. If not given, will default to the number of centers in the given arrays.
# - sphere (optional bool, default=False): Whether the points are uniform spherical data, used in title. If false, the title will specify that it is CMASS points; if true, will specify sphere in title.
# 
# Process: Plot the points in ($\alpha, \delta$) space. The label will determine the color and shape of the points. Once these are plotted, the centers will be plotted as black dots on top.
# 
# Output: Figure object. In case we want to save the plots, I need some way to access them from outside, as I'm not going to include the save option here.

# In[ ]:

def region_map(data_ra, data_dec, labels, center_ra, center_dec, area, niters, xlim=[-180,180], ylim=[-90,90], ds=True, ncen=None, sphere=False) :
    import numpy as np
    import sys
    import matplotlib.pyplot as plt
    from matplotlib import rc
    import itertools
    get_ipython().magic(u'matplotlib inline')
    
    assert data_ra.size == data_dec.size, 'Must have same number of points in RA and DEC'
    assert data_ra.size == labels.size, 'Must have same number of points as labels'
    assert center_ra.size == center_dec.size, 'Must have same number of centers in RA and DEC'
    
    rc('text', usetex=True)
    rc('font', **{'family':'serif', 'serif':'cm', 'size':12})
    
    if ncen is None :
        ncen = center_ra.size
    
    # Set iterative colors and markers
    colors = np.arange(ncen)/float(ncen)
    markers_iter = itertools.cycle((',', '+', '.', '*'))
    markers = [markers_iter.next() for i in range(ncen)]
    
    # Set up the figure
    fig = plt.figure(facecolor='w', figsize=(10,6))
    if ds :
        alpha = 0.6
        if sphere :
            fig_title = r'Sphere {:s} subsample, {} regions, {} iterations'.format(area, ncen, niters)
        else :
            fig_title = r'CMASS {:s} subsample, {} regions, {} iterations'.format(area, ncen, niters)
    else :
        alpha = 0.3
        if sphere :
            fig_title = r'Sphere {:s}, {} regions, {} iterations'.format(area, ncen, niters)
        else :
            fig_title = r'CMASS {:s}, {} regions, {} iterations'.format(area, ncen, niters)
    plt.title(fig_title)
    plt.xlabel(r'$\alpha$ [Degrees]')
    plt.ylabel(r'$\delta$ [Degrees]')
    plt.xlim(xlim)
    plt.ylim(ylim)
    
    # Plot the points
    for i in range(ncen) :
        c = colors[i]
        m = markers[i]
        plt.scatter(data_ra[(labels == i)], data_dec[(labels == i)], alpha=alpha, marker=m, color=plt.cm.rainbow(c))
    
    # And now plot the centers
    plt.plot(center_ra, center_dec, 'ko')
    
    # Return the figure object
    return (fig)

