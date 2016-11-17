#!/usr/bin/env python

## Imports
import os, sys
import numpy as np
import argparse
import kmeans_radec as kmrd
from astropy.io import fits
### Change this depending on where the module is
# ac_path = os.path.expanduser('~/autocorrelation_code/my_modules')
ac_path = os.path.abspath('my_modules')
if ac_path not in sys.path :
    sys.path.insert(0, ac_path)
import angular_correlation as corr
from datetime import datetime
import ConfigParser

# Command line arguments
parser = argparse.ArgumentParser(description='Compute correlation function from pair counts in files')
## File locations
parser.add_argument('data', help='Path to the data files. Should be a directory rather than a file!')
parser.add_argument('save', help='Path to save the files. Should be a directory rather than a file!')
## Options for redshift bins and range, separation range and bins, KMeans regions
parser.add_argument('--zmin', type=float, default=0.6, help='Minimum redshift to use (default: %(default)s)')
parser.add_argument('--zmax', type=float, default=0.8, help='Maximum redshift to use (default: %(default)s)')
parser.add_argument('--nz', type=int, default=4, help='Number of redshift bins to use (default: %(default)s)')
parser.add_argument('--min_sep', type=float, default=60.0, help='Minimum separation to use in units of Mpc (default: %(default)s)')
parser.add_argument('--max_sep', type=float, default=200.0, help='Maximum separation to use in units of Mpc (default: %(default)s)')
parser.add_argument('--nbins', type=int, default=20, help='Number of separation bins to use (default: %(default)s)')
parser.add_argument('--ncen', type=int, default=100, help='Number of KMeans regions to use (default: %(default)s)')
parser.add_argument('--zbflag', choices=['phot', 'spec'], default='phot', help='Flag for binning in photometric or spectroscopic redshift for file naming (choices: %(choices)s; default: %(default)s)')
## The scatter on the redshifts, translating to the perturbations along the los
parser.add_argument('--sigmaz', default=0.0, type=float, help='Fractional error on redshift, where 0 means the redshifts are exact and the individual redshift errors are given by sigmaz*(1+zspec)')
## An option for timing
parser.add_argument('-t', action='store_true', dest='timing', help='If flag given, return timing info')
## Notice email options
parser.add_argument('--mail_options', metavar='/path/to/mail_options.ini', help='Path to mail options file for sending notice email. If None, no email will be sent (default: %(default)s)')
parser.add_argument('--nohup', help='nohup file used (default: %(default)s)')
args = parser.parse_args()

## Start timing
start = datetime.now()

zbf = 'p' if args.zbflag == 'phot' else 's'
dpath = args.data
spath = args.save
nz = args.nz
ncen = args.ncen
timing = args.timing
spath = os.path.join(spath, 'cf_nbins{}_smin{}_smax{}_ncen{}_z{}min{}_z{}max{}_nz{}_scatter{}'.format(args.nbins, args.min_sep, args.max_sep, args.ncen, zbf, args.zmin, zbf, args.zmax, nz, args.sigmaz))
if not os.path.exists(spath) :
	os.makedirs(spath)
fbase1 = 'pc'
nfbase = 'nd_nr'
fbase2 = 'cf'
file_ext = '.dat'
hbs = (np.log(args.max_sep) - np.log(args.min_sep))/(2.*float(args.nbins))
r_par = np.logspace(np.log(args.min_sep)+hbs, np.log(args.max_sep)-hbs, num=args.nbins, base=np.e)
np.savetxt(os.path.join(spath, 'r_parallel.dat'), r_par, fmt='%-25.18f')

start_cf = datetime.now()
if nz > 1 :
	r_perp = np.loadtxt(os.path.join(dpath, '{}_zbin0_rlbin0_bin0_bin0{}'.format(fbase1, file_ext)), usecols=(0,))
	np.savetxt(os.path.join(spath, 'r_perp.dat'), r_perp, fmt='%-25.18f')
	for i in range(nz) :
		start_cf_i = datetime.now()
		xi = np.empty((r_perp.size, r_par.size), dtype=float)
		if ncen > 1 :
			err = np.empty_like(xi)
		for j in range(r_par.size) :
			start_cf_j = datetime.now()
			nd, nr = np.loadtxt(os.path.join(dpath, '{}_zbin{}_rlbin{}{}'.format(nfbase, i, j, file_ext)), unpack=True)
			if ncen == 1 :
				xi[:,j] = corr.tc_proj_corr(args.nbins, os.path.join(dpath, '{}_zbin{}_rlbin{}{}'.format(fbase1, i, j, file_ext)), nd, nr)[1]
			else :
				xi[:,j], err[:,j] = corr.tc_proj_corr_reg(args.nbins, ncen, os.path.join(dpath, '{}_zbin{}_rlbin{}'.format(fbase1, i, j)), file_ext, nd, nr)[1:3]
			del nd, nr
			if timing :
				print 'Time to compute one bin in r_parallel and z = {}'.format(datetime.now() - start_cf_j)
				sys.stdout.flush()
			del start_cf_j
		np.savetxt(os.path.join(spath, '{}_xi_zbin{}_rperp_rpar{}'.format(fbase2, i, file_ext)), xi)
		if ncen > 1 :
			np.savetxt(os.path.join(spath, '{}_err_zbin{}_rperp_rpar{}'.format(fbase2, i, file_ext)), err)
			del err
		del xi
		if timing :
			print 'Time to compute one bin z = {}'.format(datetime.now() - start_cf_i)
			sys.stdout.flush()
		del start_cf_i
else :
	r_perp = np.loadtxt(os.path.join(dpath, '{}_rlbin0_bin0_bin0{}'.format(fbase1, file_ext)), usecols=(0,))
	np.savetxt(os.path.join(spath, 'r_perp.dat'), r_perp, fmt='%-25.18f')
	xi = np.empty((r_perp.size, r_par.size), dtype=float)
	if ncen > 1 :
		err = np.empty_like(xi)
	for j in range(r_par.size) :
		start_cf_j = datetime.now()
		nd, nr = np.loadtxt(os.path.join(dpath, '{}_rlbin{}{}'.format(nfbase, j, file_ext)), unpack=True)
		if ncen == 1 :
			xi[:,j] = corr.tc_proj_corr(args.nbins, os.path.join(dpath, '{}_rlbin{}{}'.format(fbase1, j, file_ext)), nd, nr)[1]
		else :
			xi[:,j], err[:,j] = corr.tc_proj_corr_reg(args.nbins, ncen, os.path.join(dpath, '{}_rlbin{}'.format(fbase1, j)), file_ext, nd, nr)[1:3]
		if timing :
			print 'Time to compute one bin in r_parallel = {}'.format(datetime.now() - start_cf_j)
			sys.stdout.flush()
		del start_cf_j
	np.savetxt(os.path.join(spath, '{}_xi_rperp_rpar{}'.format(fbase2, file_ext)), xi)
	if ncen > 1 :
		np.savetxt(os.path.join(spath, '{}_err_rperp_rpar{}'.format(fbase2, file_ext)), err)
		del err
	del xi
if timing :
	print 'Time to compute in all bins = {}'.format(datetime.now() - start_cf)
	sys.stdout.flush()
del start_cf
	

if args.mail_options is not None :
	import noticeEMail as mail
	config2 = ConfigParser.RawConfigParser()
	config2.read(args.mail_options)
	mail_sec = 'mail_options'
	usr = config2.get(mail_sec, 'usr')
	psw = config2.get(mail_sec, 'psw')
	fromaddr = config2.get(mail_sec, 'fromaddr')
	toaddr = config2.get(mail_sec, 'toaddr')
	mail.noticeEMail(start, usr, psw, fromaddr, toaddr, nohup=args.nohup, plots=[])

else :
	print 'Time elapsed: {}'.format(datetime.now() - start)
	sys.stdout.flush()
