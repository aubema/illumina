#!/usr/bin/env python
#
# Preprocessing for Illumina
#
# Author : Alexandre Simoneau
#
# January 2019

import numpy as np
import shutil, re, os, yaml
from glob import glob
import pytools as pt
import MultiScaleData as MSD
from scipy.interpolate import interp1d as interp

dir_name = "Inputs/"
shutil.rmtree(dir_name,True)
os.makedirs(dir_name)

with open("inputs_params.in") as f:
	params = yaml.load(f)

# Angular distribution (normalised to 1)
lop_files = glob("Lights/*.lop")
angles = np.arange(181,dtype=float)
lop = { os.path.basename(s).split('_',1)[0]:pt.load_lop(angles,s) for s in lop_files }

# Spectral distribution (normalised with scotopric vision to 1)
wav, viirs = np.loadtxt("Lights/viirs.dat", skiprows=1).T
viirs = pt.spct_norm(wav,viirs)
norm_spectrum = pt.load_spct(wav, np.ones(wav.shape), "Lights/photopic.dat", 1)

spct_files = glob("Lights/*.spct")
spct = { os.path.basename(s).split('_',1)[0]:pt.load_spct(wav,norm_spectrum,s) for s in spct_files }

# lamps distribution
inv_name = "inventory.txt"
zonData = pt.parse_inventory(inv_name,7)

print "Calculating the generalized lamps."

# Calculate zones lamps
zones = pt.make_zones(angles, lop, wav, spct, zonData )

print "Splitting in a few wavelengths."
n = params['nb_bins']
lmin = params['lambda_min']
lmax = params['lambda_max']

bool_array = (lmin<=wav)*(wav<lmax)

limits = np.array(map(np.min,np.array_split(wav[bool_array],n,-1))+[lmax])

lim_file = dir_name + "integration_limits.dat"
with open(lim_file,'w') as f:
	f.write("%d\n"%n)
with open(lim_file,'ab') as f:
	np.savetxt(f,limits[:,np.newaxis])

# Create the desired lamp files
x = np.mean([limits[1:],limits[:-1]],0)
y = np.array(map(np.mean,np.array_split(zones[:,:,bool_array],n,-1),[-1]*n)).transpose(1,2,0)

print "Creating files."

try:
	os.symlink(os.path.abspath("srtm.hdf5"),dir_name+"srtm.hdf5")
except OSError as e:
	if e[0] != 17:
		raise

for l in xrange(n):
	for z in xrange(len(zones)):
		np.savetxt( dir_name+"fctem_wl_%03d_zon_%03d.dat"%(x[l],z+1), np.concatenate([ y[z,:,l],angles ]).reshape((2,-1)).T )

out_name = params['exp_name']

# Removing overlapped light
viirs_dat = MSD.Open("stable_lights.hdf5")

n = viirs_dat._attrs['scale_factor']
n1 = n * (n/2)
n2 = n1 + n
for i in range(1,len(viirs_dat)):
	viirs_dat[i][n1:n2,n1:n2] = 0

viirs_dat.save(dir_name+"stable_lights")

print "Making zone properties files."

circles = viirs_dat.copy() # Same geolocalisation
circles.clear()

zonfile = np.loadtxt(inv_name,usecols=range(7))

# zone number
for i,dat in enumerate(zonfile,1):
	circles.set_circle((dat[0],dat[1]),dat[2]*1000,i)
circles.save(dir_name+"/"+out_name+"_zone")

# obstacle height
for i,dat in enumerate(zonfile,1):
	circles.set_circle((dat[0],dat[1]),dat[2]*1000,dat[3])
circles.save(dir_name+"/"+out_name+"_obsth")

# obstacle distance
for i,dat in enumerate(zonfile,1):
	circles.set_circle((dat[0],dat[1]),dat[2]*1000,dat[4])
circles.save(dir_name+"/"+out_name+"_obstd")

# obstacle opacity
for i,dat in enumerate(zonfile,1):
	circles.set_circle((dat[0],dat[1]),dat[2]*1000,dat[5])
circles.save(dir_name+"/"+out_name+"_obstf")

# lamps height
for i,dat in enumerate(zonfile,1):
	circles.set_circle((dat[0],dat[1]),dat[2]*1000,dat[6])
circles.save(dir_name+"/"+out_name+"_altlp")

print "Linking mie files."

aero_profile = params['aerosol_profile']
RH = params['relative_humidity']
mie_pre = aero_profile+"_RH%02d" % RH

ppath = os.environ['PATH'].split(os.pathsep)
illumpath = filter(lambda s: "illumina" in s and "bin" not in s, ppath)[0]

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

with open(dir_name+"/zon.lst",'w') as zfile:
	zfile.write('\n'.join( map(lambda n:"%03d"%n, xrange(1,len(zones)+1) ))+'\n')
with open(dir_name+"/wav.lst",'w') as zfile:
	zfile.write('\n'.join( map(lambda n:"%03d"%n, x ))+'\n')

print "Interpolating reflectance."

wls = np.unique(np.concatenate([x,[700]]))
circles.clear()

refl_wav = np.loadtxt("modis.dat",usecols=[1])
refl_raw = [ MSD.Open("refl_b%02d.hdf5" % (i+1)) for i in xrange(len(refl_wav)) ]
refl_int = [ circles.copy() for i in xrange(len(wls)) ]

for b in xrange(len(refl_raw[0])):
	dat = np.asarray(map(lambda a: a[b], refl_raw))
	I = interp(modis_wav, dat, axis=0, copy=False,
			   bounds_error=False, fill_value='extrapolate')
	refl = I(wls)
	for i in xrange(len(wls)):
		refl_int[i][b] = refl[i]

for i in xrange(len(refl_int)):
	refl_int[i].save(dir_name+"modis_%03d" % wls[i])

print "Inverting lamp intensity."

#TODO

print "Done."
