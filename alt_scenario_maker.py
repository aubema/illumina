#!/usr/bin/env python

import numpy as np
from pytools import *
from glob import glob
import os

new_name = raw_input("New name : ")
old_invp = raw_input("Old inventory path : ")
inv_path = raw_input("New inventory path : ")

print "\nLoading data..."

# Angular distribution (normalised to 1)
lop_files = glob("Lights/*.lop")
angles = np.arange(181,dtype=float)
lop = { os.path.basename(s).split('_',1)[0]:load_lop(angles,s) for s in lop_files }

# Spectral distribution (normalised with scotopric vision to 1)
wav, viirs = np.loadtxt("Lights/viirs.dat", skiprows=1).T
viirs = spct_norm(wav,viirs)
scotopic = load_spct(wav, np.ones(wav.shape), "Lights/scotopic.dat", 1)
photopic = load_spct(wav, np.ones(wav.shape), "Lights/photopic.dat", 1)

ratio_ps = float(raw_input("    photopic/scotopic ratio ? (0 <= p/(s+p) <= 1) : "))
norm_spectrum = ratio_ps*photopic + (1-ratio_ps)*scotopic

spct_files = glob("Lights/*.spct")
spct = { os.path.basename(s).split('_',1)[0]:load_spct(wav,norm_spectrum,s) for s in spct_files }

zonData = parse_inventory(inv_path)
oldZonData = parse_inventory(old_invp,6)

print "\nCalculating the generalized lamps..."

# Calculate zones lamps
zones = make_zones(angles, lop, wav, spct, zonData )

# Make bins
x = np.loadtxt("Intrants/wav.lst")
n = x.size
ilims = np.genfromtxt("Intrants/integration_limits.dat",skip_header=1)
dl = ilims[1:]-ilims[:-1]
bool_array = (ilims[0]<=wav)*(wav<ilims[-1])
y = np.array(map(np.mean,np.array_split(zones[:,:,bool_array],n,-1),[-1]*n)).transpose(1,2,0)

# Photopic/scotopic spectrum
ratio_ps = float(raw_input("    photopic/scotopic ratio for lamp power ? (0 <= p/(s+p) <= 1) : "))
nspct = ratio_ps*photopic + (1-ratio_ps)*scotopic
nspct = nspct/np.sum(nspct)
nspct = np.array(map(np.mean,np.array_split(nspct[bool_array],n)))

dirname = "Intrants_"+new_name.replace(' ','_')+'/'
if not os.path.exists(dirname):
	os.makedirs(dirname)

# Link unmodified files
fctem = set(glob("Intrants/*fctem*"))
lumlp = set(glob("Intrants/*lumlp*"))
intrs = set(glob("Intrants/*"))

files = intrs-lumlp-fctem
for name in files:
	try:
		os.symlink(os.path.abspath(name), dirname+os.path.basename(name))
	except OSError as e:
		if e[0] != 17:
			raise

# Main treatement
for z in xrange(len(zones)):
	print "Treating zone : %d/%d"%(z+1,len(zones))
	lum_files = sorted(glob("Intrants/*lumlp_%03d*.pgm"%(z+1)))
	spt_files = sorted(glob("Intrants/fctem*%03d*.dat"%(z+1)))

	# Linking files if inventory unchanged (faster)
	if zonData[z]==oldZonData[z]:
		for name in lum_files+spt_files:
			os.symlink(os.path.abspath(name), dirname+os.path.basename(name))
	else:
		head,p,data = load_pgm(lum_files[0])
		data = np.asarray([load_pgm(s)[2].reshape((1,-1)) for s in lum_files])
		lumtot = np.sum(data*(dl*nspct)[:,np.newaxis,np.newaxis],0) # Total luminosity
		spct = np.sum(y[z] * 2*np.pi*np.sin(np.deg2rad(angles))[:,np.newaxis] * (angles[1]-angles[0]), 0)
		K = safe_divide(lumtot,np.sum(spct*nspct*dl))

		ndata = spct[:,np.newaxis,np.newaxis]*K[np.newaxis]

		# Saving new files
		for wl in xrange(n):
			save_pgm(dirname+os.path.basename(lum_files[wl]),head,p,ndata[wl])
			np.savetxt(dirname+os.path.basename(spt_files[wl]),np.concatenate([y[z,:,wl],angles]).reshape((2,-1)).T)

print "Done."

