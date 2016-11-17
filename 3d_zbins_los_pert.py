#!/usr/bin/env python

# Imports
import os, sys
import numpy as np
import argparse
import shlex, subprocess
### Change this depending on where the module is
# ac_path = os.path.expanduser('~/autocorrelation_code/my_modules')
ac_path = os.path.abspath('my_modules')
if ac_path not in sys.path :
    sys.path.insert(0,ac_path)
from noticeEMail import noticeEMail as mail
from datetime import datetime
from ConfigParser import ConfigParser

# Command line arguments
parser = argparse.ArgumentParser(description='Code to call the scripts for calculating \
the projected correlation function with los perturbations')
## File locations
parser.add_argument('data', help='Full path to data FITS file')
parser.add_argument('rand', help='Full path to random FITS file')
parser.add_argument('save', help='Path to save the files. Should be a directory rather \
than a file!')
parser.add_argument('cosmology', metavar='/path/to/cosmology_file.ini', help='Full path \
to cosmology ini file')
## The scatter on the redshifts, translating to the perturbations along the los
parser.add_argument('--sigmaz', default=0.0, type=float, help='Fractional error on \
redshift, where 0 means the redshifts are exact and the individual redshift errors are \
given by sigmaz*(1+zspec)')
## The tags to use for accessing redshifts
parser.add_argument('--dztag', default='ZPHOT', help='The tag to use to access the \
redshifts in the data catalog (default: %(default)s)')
parser.add_argument('--dzbtag', default='ZPHOT', help='The tag to use for redshift \
binning in the data (default: %(default)s)')
parser.add_argument('--rztag', default='ZPHOT', help='The tag to use to access the \
redshifts in the random catalog (default: %(default)s)')
parser.add_argument('--rzbtag', default='ZPHOT', help='The tag to use for redshift \
binning in the randoms (default: %(default)s)')
parser.add_argument('--zbflag', choices=['phot', 'spec'], default='phot', help='Flag for \
binning in photometric or spectroscopic redshift for file naming (choices: %(choices)s; \
default: %(default)s)')
## Options for redshift bins and range, separation range and bins, KMeans regions
parser.add_argument('--zmin', type=float, default=0.6, help='Minimum redshift to use \
(default: %(default)s)')
parser.add_argument('--zmax', type=float, default=0.8, help='Maximum redshift to use \
(default: %(default)s)')
parser.add_argument('--nz', type=int, default=1, help='Number of redshift bins to use \
(default: %(default)s)')
parser.add_argument('--min_sep', type=float, default=60.0, help='Minimum separation to \
use in units of Mpc (default: %(default)s)')
parser.add_argument('--max_sep', type=float, default=200.0, help='Maximum separation to \
use in units of Mpc (default: %(default)s)')
parser.add_argument('--nbins', type=int, default=20, help='Number of separation bins to \
use (default: %(default)s)')
parser.add_argument('--ncen', type=int, default=100, help='Number of KMeans regions to \
use (default: %(default)s)')
## Use this argument if the data is not continuous in RA
parser.add_argument('--wrap', type=float, default=None, help='Where to wrap RA to make \
contiguous data set (default: %(default)s')
## An option for timing
parser.add_argument('-t', action='store_true', dest='timing', help='If flag given, \
return timing info')
## Options to skip things that have already been done
parser.add_argument('-p', action='store_true', help='If flag specified, skip pair \
counting (assuming pair counting already completed at the correct file location)')
parser.add_argument('-c', action='store_true', help='If flag specified, skip correlation \
function calculation (assuming calculation already completed at the correct file location)')
## Option to make plot: if anything is given for a plot directory, make plot.
## Else, skip this.
parser.add_argument('--plot', default=None, help='Directory in which to save plot. If \
None, plot will not be made (default: %(default)s)')
## Notice email options
parser.add_argument('--mail_options', metavar='/path/to/mail_options.ini', help='Path to \
mail options file for sending notice email. If None, no email will be sent (default: \
%(default)s)')
parser.add_argument('--nohup', help='nohup file used (default: %(default)s)')

## Parse command line arguments
args = parser.parse_args()

# Start timing and get a few options from the command line
start_time = datetime.now()
timing = args.timing
zbf = 'p' if args.zbflag == 'phot' else 's'
if args.mail_options is not None :
	mail_config = ConfigParser()
	mail_config.read(args.mail_options)
	mail_ops = mail_config._sections['mail_options']

# Just so we have this available
pc_save = os.path.join(args.save, 'pc_nbins{}_smin{}_smax{}_b{}_ncen{}_z{}min{}_z{}max{}\
_nz{}_scatter{}'.format(args.nbins, args.min_sep, args.max_sep, 0.05, args.ncen, zbf, 
args.zmin, zbf, args.zmax, args.nz, args.sigmaz))
cf_save = os.path.join(args.save, 'cf_nbins{}_smin{}_smax{}_ncen{}_z{}min{}_z{}max{}\
_nz{}_scatter{}'.format(args.nbins, args.min_sep, args.max_sep, args.ncen, zbf, 
args.zmin, zbf, args.zmax, args.nz, args.sigmaz))

# Run pair count code
if not args.p :
	start_pc = datetime.now()
	## Get the command to run the script
	pc_args = shlex.split('nice -n 19 python 3d_pc_zbins_los_v1.py {} {} {} {} '\
	'--dztag {} --dzbtag {} --rztag {} --rzbtag {} --zbflag {} --zmin {} --zmax {} --nz '\
	'{} --min_sep {} --max_sep {} --nbins {} --ncen {}{} --sigmaz {}{}'.format(args.data, \
	args.rand, args.save, args.cosmology, args.dztag, args.dzbtag, \
	args.rztag, args.rzbtag, args.zbflag, args.zmin, args.zmax, args.nz, args.min_sep, args.max_sep, \
	args.nbins, args.ncen, ' --wrap '+args.wrap if args.wrap is not None else '', \
	args.sigmaz, ' -t' if timing else ''))
	## Run the command and pipe output to a file
	with open('pc_proj_los.txt', 'wb') as out :
		p_pc = subprocess.Popen(pc_args, stdout=out, stderr=out)
		### This communicates with the subprocess so we wait until it finishes
		p_pc.communicate()
	## Print timing if desired
	if timing :
		print 'Time to run pair counting = {}'.format(datetime.now() - start_pc)
		sys.stdout.flush()
	## Alert user if desired
	if args.mail_options is not None :
		mail(start_pc, mail_ops['usr'], mail_ops['psw'], mail_ops['fromaddr'], \
		mail_ops['toaddr'], exit=p_pc.returncode, nohup='pc_proj_los.txt')
	del start_pc, pc_args, p_pc
else :
	print 'Pair counting skipped'
	sys.stdout.flush()

# Run correlation function code
if not args.c :
	start_cf = datetime.now()
	## Get the command to run the script
	cf_args = shlex.split('nice -n 19 python 3d_cf_from_pc_los_v1.py {} {} --zmin {} '\
	'--zmax {} --nz {} --min_sep {} --max_sep {} --nbins {} --ncen {}  --sigmaz {}'\
	'{}'.format(pc_save, args.save, args.zmin, args.zmax, args.nz, args.min_sep, \
	args.max_sep, args.nbins, args.ncen, args.sigmaz, ' -t' if timing else ''))
	## Run the command and pipe output to a file
	with open('cf_proj_los.txt', 'wb') as out :
		p_cf = subprocess.Popen(cf_args, stdout=out, stderr=out)
		### This communicates with the subprocess so we wait until it finishes
		p_cf.communicate()
	## Print timing if desired
	if timing :
		print 'Time to run correlation function = {}'.format(datetime.now() - start_cf)
		sys.stdout.flush()
	## Alert user if desired
	if args.mail_options is not None :
		mail(start_cf, mail_ops['usr'], mail_ops['psw'], mail_ops['fromaddr'], \
		mail_ops['toaddr'], exit=p_cf.returncode, nohup='cf_proj_los.txt')
	del start_cf, cf_args, p_cf
else :
	print 'Correlation function skipped'
	sys.stdout.flush()

## Make plot if desired
if args.plot is not None :
	start_plot = datetime.now()
	## Get the command to run the script
	plot_args = shlex.split('python plot_xi_comp.py {} cf_xi_rperp_rpar.dat {} --sperp '\
	'r_perp.dat --spar r_parallel.dat --zmin {} --zmax {} --nz {} --min_sep {} '\
	'--max_sep {} --nbins {} --ncen {}'.format(cf_save, args.plot, args.zmin, args.zmax, \
	args.nz, args.min_sep, args.max_sep, args.nbins, args.ncen))
	## Run the command and pipe output to a file
	with open('plot_proj_los.txt', 'wb') as out :
		p_plot = subprocess.Popen(plot_args, stdout=out, stderr=out)
		### This communicates with the subprocess so we wait until it finishes
		p_plot.communicate()
	## Print timing if desired
	if timing :
		print 'Time to make plot = {}'.format(datetime.now() - start_plot)
		sys.stdout.flush()
	## Alert user if desired
	if args.mail_options is not None :
		### Let's get the list of plots that we made to email those also
		plot_list = [os.path.join(args.plot, 'smin{}_smax{}_nbins{}_ncen{}_xi_zbin{}.'\
		'png'.format(args.min_sep, args.max_sep, args.nbins, args.ncen, i)) for i in \
		range(args.nz)]
		mail(start_plot, mail_ops['usr'], mail_ops['psw'], mail_ops['fromaddr'], \
		mail_ops['toaddr'], exit=p_plot.returncode, nohup='plot_proj_los.txt', \
		plots=plot_list)
		del plot_list
	del start_plot, plot_args, p_plot
else :
	print 'No plot made'
	sys.stdout.flush()
	

# End of program: alert user if desired. Otherwise, print total time elapsed
if args.mail_options is not None :
	mail(start_time, mail_ops['usr'], mail_ops['psw'], mail_ops['fromaddr'], \
	mail_ops['toaddr'], exit=0, nohup=args.nohup)
else :
	print 'Elapsed time = {}'.format(datetime.now() - start_time)
