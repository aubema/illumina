#!/usr/bin/env python

import numpy as np
from pytools import *
from glob import glob
import os

new_name = raw_input("New name : ")
inv_path = raw_input("New inventory path : ")

print "\nLoading data..."

# Angular distribution (normalised to 1)
lop_files = glob("Lights/*.lop")
angles = np.loadtxt(lop_files[0])[:,1]
lop = { os.path.basename(s).split('_',1)[0]:LOP_norm(angles,np.loadtxt(s)[:,0]) for s in lop_files }

# Spectral distribution (normalised with scotopric vision to 1)
wavelenght, viirs = np.loadtxt("Lights/viirs.dat", skiprows=1).T
scotopic = np.loadtxt("Lights/scotopic.dat", skiprows=1)[:,1]
photopic = np.loadtxt("Lights/photopic.dat", skiprows=1)[:,1]

ratio_ps = float(raw_input("    photopic/scotopic ratio for inventory ? (0 <= p/(s+p) <= 1) : "))
norm_spectrum = ratio_ps*photopic + (1-ratio_ps)*scotopic
norm_spectrum = norm_spectrum/np.sum(norm_spectrum)

spct_files = glob("Lights/*.spct")
spct = { os.path.basename(s).split('_',1)[0]:SPD_norm(wavelenght,norm_spectrum,np.loadtxt(s,skiprows=1)[:,1]) for s in spct_files }

zonData = parse_inventory(inv_path)
oldZonData = parse_inventory("inventaire.txt")

print "\nCalculating the generalized lamps..."

# Calculate zones lamps
zones = np.asarray([ zon_norm(angles, wavelenght, sum(l[0]*spct[l[1]]*lop[l[2]][:,np.newaxis] for l in lampData)) for lampData in zonData ])

x = np.loadtxt("Intrants/wav.lst")
n = x.size
ilims = np.genfromtxt("Intrants/integration_limits.dat",skiprows=1)
dl = ilims[1:]-ilims[:-1]
bool_array = (ilims[0]<=wavelenght)*(wavelenght<ilims[-1])
y = np.array(map(np.mean,np.array_split(zones[:,:,bool_array],n,-1),[-1]*n)).transpose(1,2,0)

ratio_ps = float(raw_input("    photopic/scotopic ratio for lamp power ? (0 <= p/(s+p) <= 1) : "))
nspct = ratio_ps*photopic + (1-ratio_ps)*scotopic
nspct = nspct/np.sum(nspct)
nspct = np.array(map(np.mean,np.array_split(nspct[bool_array],n)))

dirname = "Intrants_"+new_name.replace(' ','_')+'/'
if not os.path.exists(dirname):
	os.makedirs(dirname)

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

for z in xrange(len(zones)):
	print "Treating zone : %d/%d"%(z+1,len(zones))
	lum_files = sorted(glob("Intrants/*lumlp_%03d*.pgm"%(z+1)))
	spt_files = sorted(glob("Intrants/fctem*%03d*.dat"%(z+1)))
	if zonData[z]==oldZonData[z]:
		for name in lum_files+spt_files:
			os.symlink(os.path.abspath(name), dirname+os.path.basename(name))
	else:
		head,p,data = load_pgm(lum_files[0])
		data = np.asarray([load_pgm(s)[2].reshape((1,-1)) for s in lum_files])
		lumtot = np.sum(data*(dl*nspct)[:,np.newaxis,np.newaxis],0)
		spct = np.sum(y[z] * 2*np.pi*np.sin(np.deg2rad(angles))[:,np.newaxis] * (angles[1]-angles[0]), 0)
		K = safe_divide(lumtot,np.sum(spct*nspct*dl))

		ndata = spct[:,np.newaxis,np.newaxis]*K[np.newaxis]

		for wl in xrange(n):
			save_pgm(dirname+os.path.basename(lum_files[wl]),head,p,ndata[wl])
			np.savetxt(dirname+os.path.basename(spt_files[wl]),np.concatenate([y[z,:,wl],angles]).reshape((2,-1)).T)

print "Done."

