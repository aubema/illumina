#!/usr/bin/env python2
#
# Preprocessing for Illumina
#
# Author : Alexandre Simoneau
#
# January 2019

import numpy as np
import shutil, re, os, yaml
from glob import glob
import pytools as pt, hdftools as hdf
import MultiScaleData as MSD
from scipy.interpolate import interp1d as interp
from itertools import izip

dir_name = "Inputs/"
shutil.rmtree(dir_name,True)
os.makedirs(dir_name)

with open("inputs_params.in") as f:
	params = yaml.safe_load(f)

# TODO: Check if lamp is in non-zero zone
if params['zones_inventory'] is not None and \
	params['lamps_inventory'] is not None:
	lamps = np.loadtxt(params['lamps_inventory'],usecols=[0,1])
	zones = np.loadtxt(params['zones_inventory'],usecols=[0,1,2])
	zonData = pt.parse_inventory(params['zones_inventory'],7)

	hasLights = [ sum( x[0] for x in z ) != 0 for z in zonData ]

	circles = hdf.from_domain("domain.ini")
	for dat,b in izip(zones,hasLights):
		circles.set_circle((dat[0],dat[1]),dat[2]*1000,b)

	zones_ind = hdf.from_domain("domain.ini")
	for i,dat in enumerate(zones,1):
		zones_ind.set_circle((dat[0],dat[1]),dat[2]*1000,i)

	for lat,lon in lamps:
		for i in xrange(len(circles)):
			try:
				col,row = circles._get_col_row((lat,lon),i)
				if circles[i][row,col] and col >= 0 and row >= 0:
					zon_ind = zones_ind[i][row,col]
					raise ValueError("WARNING: Point source at (%g,%g) falls within non-null zone #%d." % (lat,lon,zon_ind))
			except IndexError:
				continue

# Angular distribution (normalised to 1)
lop_files = glob("Lights/*.lop")
lop = {
	os.path.basename(s).split('_',1)[0] : \
	pt.load_lop(angles,s) for s in lop_files }

# Spectral distribution (normalised with scotopric vision to 1)
wav, viirs = np.loadtxt( "Lights/viirs.dat", skiprows=1 ).T
viirs = pt.spct_norm(wav,viirs)
norm_spectrum = pt.load_spct(
	wav,
	np.ones(wav.shape),
	"Lights/photopic.dat",
	1 )

spct_files = glob("Lights/*.spct")
spct = {
	os.path.basename(s).split('_',1)[0] : \
	pt.load_spct(wav,norm_spectrum,s) \
	for s in spct_files }

print "Splitting in a few wavelengths."
n_bins = params['nb_bins']
lmin = params['lambda_min']
lmax = params['lambda_max']

bool_array = (lmin<=wav)*(wav<lmax)

limits = np.array(map(
	np.min, np.array_split(
		wav[bool_array],
		n_bins,
		-1 ) ) + [lmax] )

lim_file = dir_name + "integration_limits.dat"
with open(lim_file,'w') as f:
	f.write("%d\n"%n_bins)
with open(lim_file,'ab') as f:
	np.savetxt(f,limits[:,np.newaxis])

x = np.mean([limits[1:],limits[:-1]],0)

print "Linking mie files."

aero_profile = params['aerosol_profile']
RH = params['relative_humidity']
mie_pre = aero_profile+"_RH%02d" % RH

ppath = os.environ['PATH'].split(os.pathsep)
illumpath = filter(
	lambda s: "illumina" in s and "bin" not in s,
	ppath )[0]

mie_path = illumpath + "/Aerosol_optical_prop/"
mie_files = glob(mie_path+mie_pre+"*.mie.out")
mie_files = { int(s.split('.')[-3][:3]):s for s in mie_files }
mie_wl = np.asarray(sorted(mie_files.keys()))
wl2mie = np.asarray([min(mie_wl, key=lambda i: abs(i-j)) for j in x])

for i in xrange(len(wl2mie)):
	name = dir_name+mie_pre.strip('_')+"_0.%03d0um.mie.out"%x[i]
	try:
		os.symlink(os.path.abspath(mie_files[wl2mie[i]]),name)
	except OSError as e:
		if e[0] != 17:
			raise

with open(dir_name+"/wav.lst",'w') as zfile:
	zfile.write('\n'.join( map(lambda n:"%03d"%n, x ))+'\n')

# TODO: Execution logic

# TODO: Change to lamp instead of zones
execfile(os.path.join(illumpath,"make_zones.py"))

# TODO: Make from lamp inventory

# TODO: Add both inputs

print "Done."
