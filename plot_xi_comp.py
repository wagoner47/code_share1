#!/usr/bin/env python

import os, sys
import numpy as np
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc
import corner
# ac_path = os.path.expanduser('~/autocorrelation_code/my_modules')
ac_path = os.path.abspath('my_modules')
if ac_path not in sys.path :
    sys.path.insert(0,ac_path)
import ConfigParser
from scipy.interpolate import UnivariateSpline
import emcee
from datetime import datetime

parser = argparse.ArgumentParser(description='Plot correlation functions in terms of components of separation')
## File paths and catalog type
parser.add_argument('load_dir', help='Directory in which data, including separation bin files, are located')
parser.add_argument('data', help='Data file name. If only one redshift bin, may be full name.')
parser.add_argument('save_dir', help='Full path to directory in which to save the output plot files', metavar='/path/to/save/plots')
group1 = parser.add_mutually_exclusive_group(required=True)
group1.add_argument('--s', help='Name of file containing bins in 3D separation, s, in Mpc')
group1.add_argument('--sperp', help='Name of file containing bins in perpendicular separation, s_perp, in Mpc')
group2 = parser.add_mutually_exclusive_group(required=True)
group2.add_argument('--mu', help='Name of file containing bins in mu')
group2.add_argument('--spar', help='Name of file containing bins in parallel separation, s_par, in Mpc')
## Details for data: down sample, snapshot, min sep, max sep, nbins, ncen, flags (probably not needed anymore)
parser.add_argument('--zmin', type=float, default=0.6, help='Minimum redshift (default: %(default)s)')
parser.add_argument('--zmax', type=float, default=0.8, help='Maximum redshift (default: %(default)s)')
parser.add_argument('--nz', type=int, default=1, help='Number of redshift bins (default: %(default)s)')
parser.add_argument('--min_sep', default=60.0, type=float, help='Minimum separation in Mpc (default: %(default)s)')
parser.add_argument('--max_sep', default=200.0, type=float, help='Maximum separation in Mpc (default: %(default)s)')
parser.add_argument('--nbins', type=int, default=20, help='Number of log-spaced separation bins (default: %(default)s)')
parser.add_argument('--ncen', type=int, default=100, help='Number of jackknife regions to use. (default: %(default)s)')
## Alert options
parser.add_argument('--mail_options', metavar='/path/to/mail_options.ini', help='Path to mail options file for sending notice email. If None, no email will be sent (default: %(default)s)')
parser.add_argument('--nohup', help='nohup file used (default: %(default)s)')
args = parser.parse_args()
nz = args.nz
zmin = args.zmin
zmax = args.zmax
min_sep = args.min_sep
max_sep = args.max_sep
assert os.path.exists(args.load_dir), 'Data directory not found: {}'.format(args.load_dir)
assert os.path.isdir(args.load_dir), 'Invalid path specified for data directory: {}'.format(args.load_dir)
load_dir = args.load_dir
if os.path.basename(args.data) == args.data :
	dpath = os.path.join(load_dir, args.data)
else :
	dpath = args.data
save_dir = args.save_dir
assert (args.s is not None and args.mu is not None) or (args.sperp is not None and args.spar is not None), 'Invalid separation bin combination. Valid cominations are s and mu or s_perp and s_par'

start = datetime.now()

## Get the rest of the file path details set up
xi = np.loadtxt(dpath)
if args.s is not None :
	if os.path.basename(args.s) == args.s :
		s = np.loadtxt(os.path.join(load_dir, args.s))
	else :
		s = np.loadtxt(args.s)
	if os.path.basename(args.mu) == args.mu :
		mu = np.loadtxt(os.path.join(load_dir, args.mu))
	else :
		mu = np.loadtxt(args.mu)
	if not xi.shape == (s.size, mu.size) :
		xi = xi.T
	S, MU = np.meshgrid(s, mu, indexing='ij', sparse=True)
	s_par = S.copy()*MU.copy()
	s_perp = S.copy()*np.sqrt(1. - np.power(MU.copy(), 2))
	del s, mu, S, MU
else :
	if os.path.basename(args.sperp) == args.sperp :
		S_perp = np.loadtxt(os.path.join(load_dir, args.sperp))
	else :
		S_perp = np.loadtxt(args.sperp)
	if os.path.basename(args.spar) == args.spar :
		S_par = np.loadtxt(os.path.join(load_dir, args.spar))
	else :
		S_par = np.loadtxt(args.spar)
	if not xi.shape == (S_perp.size, S_par.size) :
		xi = xi.T
	s_perp, s_par = np.meshgrid(S_perp, S_par, indexing='ij')
	del S_perp, S_par
if not os.path.exists(save_dir) :
	os.makedirs(save_dir)
plot_list = []

rc('text', usetex=True)
font_dict = {'family':'serif', 'serif':'cm', 'size':16}
rc('font', **font_dict)

## Make plot
if nz > 1 :
	for i in range(nz) :
		pname = 'smin{}_smax{}_nbins{}_ncen{}_xi_zbin{}.png'.format(min_sep, max_sep, args.nbins, args.ncen, i)
		plt.figure(figsize=(10,6), facecolor='w')
		plt.xlabel(r'$r_p$ (Mpc)')
		plt.ylabel(r'$\pi$ (Mpc)')
		plt.contour(s_perp, s_par, xi, cmap='YlGnBu')
		plt.contourf(s_perp, s_par, xi, cmap='YlGnBu')
		cbl = plt.colorbar()
		cbl.set_label(r'$\xi(r_p, \pi)$')
		plt.tight_layout()
		plt.savefig(os.path.join(save_dir, pname))
		plot_list.append(os.path.join(save_dir, pname))
		plt.close()
		del pname, cbl
else :
	pname = 'smin{}_smax{}_nbins{}_ncen{}_xi_zbin{}.png'.format(min_sep, max_sep, args.nbins, args.ncen, 0)
	plt.figure(figsize=(10,6), facecolor='w')
	plt.xlabel(r'$r_p$ (Mpc)')
	plt.ylabel(r'$\pi$ (Mpc)')
	plt.contour(s_perp, s_par, xi, cmap='YlGnBu')
	plt.contourf(s_perp, s_par, xi, cmap='YlGnBu')
	cbl = plt.colorbar()
	cbl.set_label(r'$\xi(r_p, \pi)$')
	plt.tight_layout()
	plt.savefig(os.path.join(save_dir, pname))
	plot_list.append(os.path.join(save_dir, pname))
	plt.close()
	del pname, cbl

## Alert user
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
