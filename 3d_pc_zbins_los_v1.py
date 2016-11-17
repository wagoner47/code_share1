#!/usr/bin/env python

## Imports
import time, os, sys
import numpy as np
import argparse
from astropy.io import fits
### Change this depending on where the module is
# ac_path = os.path.expanduser('~/autocorrelation_code/my_modules')
ac_path = os.path.abspath('my_modules')
if ac_path not in sys.path :
    sys.path.insert(0,ac_path)
import angular_correlation as corr
# import kmeans_radec as kmrd
import jackknife as jk
from datetime import datetime
import ConfigParser
from astropy.cosmology import w0waCDM, FlatLambdaCDM
from scipy.interpolate import UnivariateSpline

parser = argparse.ArgumentParser(description='Complete projected correlation calculation for data with perturbations along the los (redshift uncertainties)')
## File locations
parser.add_argument('data', help='Full path to data FITS file', metavar='/path/to/data_file.fit')
parser.add_argument('rand', help='Full path to random FITS file', metavar='/path/to/rand_file.fit')
parser.add_argument('save', help='Path to save the files. Should be a directory rather than a file!', metavar='/path/to/save/files')
parser.add_argument('cosmology', metavar='/path/to/cosmology_file.ini', help='Full path to cosmology ini file')
## The tags to use for accessing redshifts
parser.add_argument('--dztag', default='ZPHOT', help='The tag to use to access the redshifts in the data catalog (default: %(default)s)')
# parser.add_argument('--dzetag', help='The tag to use to access the redshift errors in the data catalog. (default: [dztag]_E)')
parser.add_argument('--dzbtag', default='ZPHOT', help='The tag to use for redshift binning in the data (default: %(default)s)')
parser.add_argument('--rztag', default='ZPHOT', help='The tag to use to access the redshifts in the random catalog (default: %(default)s)')
# parser.add_argument('--rzetag', help='The tag to use to access the redshift errors in the random catalog. (default: [rztag]_E)')
parser.add_argument('--rzbtag', default='ZPHOT', help='The tag to use for redshift binning in the randoms (default: %(default)s)')
parser.add_argument('--zbflag', choices=['phot', 'spec'], default='phot', help='Flag for binning in photometric or spectroscopic redshift for file naming (choices: %(choices)s; default: %(default)s)')
## Options for redshift bins and range, separation range and bins, KMeans regions
parser.add_argument('--zmin', type=float, default=0.6, help='Minimum redshift to use (default: %(default)s)')
parser.add_argument('--zmax', type=float, default=0.8, help='Maximum redshift to use (default: %(default)s)')
parser.add_argument('--nz', type=int, default=4, help='Number of redshift bins to use (default: %(default)s)')
parser.add_argument('--min_sep', type=float, default=60.0, help='Minimum separation to use in units of Mpc (default: %(default)s)')
parser.add_argument('--max_sep', type=float, default=200.0, help='Maximum separation to use in units of Mpc (default: %(default)s)')
parser.add_argument('--nbins', type=int, default=20, help='Number of separation bins to use (default: %(default)s)')
parser.add_argument('--ncen', type=int, default=100, help='Number of KMeans regions to use (default: %(default)s)')
## Use this argument if the data is not continuous in RA
parser.add_argument('--wrap', type=float, default=None, help='Where to wrap RA to make contiguous data set (default: %(default)s')
## An option for timing
parser.add_argument('-t', action='store_true', dest='timing', help='If flag given, return timing info')
## Bug testing: this factor to multiply the errors should not be used normally
parser.add_argument('--sigmaz', default=0.0, type=float, help='Fractional error on redshift, where 0 means the redshifts are exact and the individual redshift errors are given by sigmaz*(1+zspec)')
## Notice email options
parser.add_argument('--mail_options', metavar='/path/to/mail_options.ini', help='Path to mail options file for sending notice email. If None, no email will be sent (default: %(default)s)')
parser.add_argument('--nohup', help='nohup file used (default: %(default)s)')
args = parser.parse_args()

## For timing purposes
start = datetime.now()
## Define input file paths
assert os.path.splitext(args.data)[1] == '.fit' or os.path.splitext(args.data)[1] == '.fits', 'Invalid file extension for data: {}. Data should be stored in a FITS file with extension .fit or .fits'.format(os.path.splitext(args.data)[1])
dpath = args.data
assert os.path.splitext(args.rand)[1] == '.fit' or os.path.splitext(args.rand)[1] == '.fits', 'Invalid file extension for randoms: {}. Randoms should be stored in a FITS file with extension .fit or .fits'.format(os.path.splitext(args.rand)[1])
rpath = args.rand
assert not os.path.isfile(args.save), 'Invalid save location: {}. Path to save should be a directory, not a file'.format(args.save)
save_base = args.save
assert os.path.splitext(args.cosmology)[1] == '.ini', 'Invalid file extension for cosmology: {}. Cosmology file must have extension .ini'.format(os.path.splitext(args.cosmology)[1])
cosmology = args.cosmology
dztag = args.dztag
# dzetag = args.dzetag if args.dzetag is not None else dztag + '_E''
dzbtag = args.dzbtag
rztag = args.rztag
# rzetag = args.rzetag if args.rzetag is not None else rztag + '_E'
rzbtag = args.rzbtag
wrap = args.wrap
zbf = 'p' if args.zbflag == 'phot' else 's'
timing = args.timing

# Setup
## Redshift range
zmin = args.zmin
zmax = args.zmax
nz = args.nz
deltaz = (zmax - zmin)/float(2*nz)
zbins = np.linspace(zmin+deltaz, zmax-deltaz, num=nz)
## KMeans
ncen = args.ncen
maxiter = 100
tol = 1e-5
size = int(3e5)
## TreeCorr
min_sep = args.min_sep
max_sep = args.max_sep
nbins = args.nbins
b = 0.05
bin_slop = min((b*float(nbins))/np.log(float(max_sep)/float(min_sep)), 1.)
config = {'min_sep':min_sep,
		  'max_sep':max_sep,
		  'nbins':nbins,
		  'bin_slop':bin_slop}
### The r_parallel grid that we will want to use
r_par = np.logspace(np.log(min_sep), np.log(max_sep), num=nbins+1, base=np.e)
## Redshift-distance conversion spline
### First, get the cosmology
config1 = ConfigParser.RawConfigParser()
config1.read(cosmology)
cos_param_sec = 'cosmological_parameters'
h0 = config1.getfloat(cos_param_sec, 'h0')
omega_m = config1.getfloat(cos_param_sec, 'omega_m')
omega_k = config1.getfloat(cos_param_sec, 'omega_k')
omega_L = 1.0 - omega_m - omega_k
w0 = config1.getfloat(cos_param_sec, 'w')
wa = config1.getfloat(cos_param_sec, 'wa')
### Now set up the astropy cosmology object
if w0 != -1.0 or wa != 0.0 :
	cosmo = w0waCDM(100*h0, omega_m, omega_L, w0=w0, wa=wa)
else :
	cosmo = FlatLambdaCDM(100*h0, omega_m)
### Finally, get distances at a grid of redshfits and use a spline to interpolate
table_z, delta_z = np.linspace(0.0, 2.0, num=101, retstep=True)
table_r = cosmo.comoving_distance(table_z).value
dist = UnivariateSpline(table_z, table_r, s=0)
invdist = UnivariateSpline(table_r, table_z, s=0)

# Define output file paths
save_dir = os.path.join(save_base, 'pc_nbins{}_smin{}_smax{}_b{}_ncen{}_z{}min{}_z{}max{}_nz{}_scatter{}'.format(nbins, min_sep, max_sep, b, ncen, zbf, zmin, zbf, zmax, nz, args.sigmaz))
if not os.path.exists(save_dir) :
	os.makedirs(save_dir)
fbase = 'pc'
nbase = 'nd_nr'
np.savetxt(os.path.join(save_dir, 'r_parallel.dat'), r_par, fmt='%-25.18f')
	
# Read in from catalogs
start_read = datetime.now()
## Data
data = fits.getdata(dpath)
dfilt = ((data[dzbtag] >= zmin) & (data[dzbtag] <= zmax))
dra = data['RA'][dfilt].copy()
ddec = data['DEC'][dfilt].copy()
odz = data[dztag][dfilt].copy()
# try :
# 	odze = data[dzetag][dfilt].copy()
# except KeyError :
# 	pass
dzb = data[dzbtag][dfilt].copy()
## Randoms
rand = fits.getdata(rpath)
rfilt = ((rand[rzbtag] >= zmin) & (rand[rzbtag] <= zmax))
rra = rand['RA'][rfilt].copy()
rdec = rand['DEC'][rfilt].copy()
orz = rand[rztag][rfilt].copy()
# try :
# 	orze = rand[rzetag][rfilt].copy()
# except KeyError :
# 	pass
rzb = rand[rzbtag][rfilt].copy()
## Alert user
print 'Data and randoms read in'
if timing :
	print 'Time to read files = {}'.format(datetime.now() - start_read)
sys.stdout.flush()
del start_read

# If more than one jackknife region wanted, find centers and labels
if ncen > 1 :
	start_reg = datetime.now()
	## Find centers
	np.random.seed(0)
	idx = np.random.randint(rra.size, size=size)
	x_samp = np.array([rra[idx].copy(), rdec[idx].copy()]).T
	if wrap is not None :
		x_samp[np.where(x_samp[:,0] > wrap),0] -= 360
	cent = jk.find_kmeans_centers(ncen, maxiter, tol, x_samp)
	del x_samp
	## Data labels
	x = np.array([dra.copy(), ddec.copy()]).T
	if wrap is not None :
		x[np.where(x[:,0] > wrap),0] -= 360
	dlabels = jk.find_kmeans_labels(x, cent)
	del x
	## Random labels
	x = np.array([rra.copy(), rdec.copy()]).T
	if wrap is not None :
		x[np.where(x[:,0] > wrap),0] -= 360
	rlabels = jk.find_kmeans_labels(x, cent)
	del x
	## Alert user
	print 'KMeans regions found'
	if timing :
		print 'Time to find KMeans regions = {}'.format(datetime.now() - start_reg)
	sys.stdout.flush()
	del start_reg

# Start pair counting
start_pc = datetime.now()
for j in range(nbins) :
	start_iter = datetime.now()
	if nz > 1 :
		for i, zbini in zip(range(nz), zbins) :
			start_zbin = datetime.now()
			fname = os.path.join(save_dir, fbase+'_zbin{}_rlbin{}'.format(i,j))
			nfname = os.path.join(save_dir, nbase+'_zbin{}_rlbin{}.dat'.format(i,j))
			zmini = zbini - deltaz
			zmaxi = zbini + deltaz
			dzfilt = ((dzb >= zmini) & (dzb <= zmaxi))
			rzfilt = ((rzb >= zmini) & (rzb <= zmaxi))
			## Do pair count over regions if desired, or single count if not
			if ncen > 1 :
				corr.paircount_treecorr_regions(dra[dzfilt], ddec[dzfilt], dlabels[dzfilt], rra[rzfilt], rdec[rzfilt], rlabels[rzfilt], ncen, nbins, config, fname, '.dat', nfname, data_r=dist(odz[dzfilt]), rand_r=dist(orz[rzfilt]), min_rpar=r_par[j], max_rpar=r_par[j+1])
			else :
				corr.paircount_treecorr(dra[dzfilt], ddec[dzfilt], rra[rzfilt], rdec[rzfilt], config, fname+'.dat', nfname, data_r=dist(odz[dzfilt]), rand_r=dist(orz[rzfilt]), min_rpar=r_par[j], max_rpar=r_par[j+1])
			print 'Pair count finished: zbin {}, r_parallel bin {}'.format(zbini, j)
			if timing :
				print 'Time for pair counting single iteration and single redshift bin = {}'.format(datetime.now() - start_zbin)
			sys.stdout.flush()
			del start_zbin
	else :
		fname = os.path.join(save_dir, fbase+'_rlbin{}'.format(j))
		nfname = os.path.join(save_dir, nbase+'_rlbin{}.dat'.format(j))
		if ncen > 1 :
			corr.paircount_treecorr_regions(dra, ddec, dlabels, rra, rdec, rlabels, ncen, nbins, config, fname, '.dat', nfname, data_r=dist(odz), rand_r=dist(orz), min_rpar=r_par[j], max_rpar=r_par[j+1])
		else :
			corr.paircount_treecorr(dra, ddec, rra, rdec, config, fname+'.dat', nfname, data_r=dist(odz), rand_r=dist(orz), min_rpar=r_par[j], max_rpar=r_par[j+1])
	if timing :
		print 'Time for pair counting single iteration = {}'.format(datetime.now() - start_iter)
	sys.stdout.flush()
	del start_iter
## Alert user
if timing :
	print 'Time for pair counting {} iterations = {}'.format(nbins, datetime.now() - start_pc)
sys.stdout.flush()
del start_pc

if args.mail_options is not None :
	import noticeEMail as mail
	config2 = ConfigParser.RawConfigParser()
	config2.read(args.mail_options)
	mail_sec = 'mail_options'
	usr = config2.get(mail_sec, 'usr')
	psw = config2.get(mail_sec, 'psw')
	fromaddr = config2.get(mail_sec, 'fromaddr')
	toaddr = config2.get(mail_sec, 'toaddr')
	mail.noticeEMail(start, usr, psw, fromaddr, toaddr, nohup=args.nohup, plots=plot_list)

else :
	print 'Time elapsed: {}'.format(datetime.now() - start)
	sys.stdout.flush()
