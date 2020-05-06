#!/usr/bin/env python3
#
# Preprocessing for Illumina
#
# Author : Alexandre Simoneau
#
# January 2019

import click
import numpy as np
import shutil, re, os, yaml
from glob import glob
import pytools as pt
import MultiScaleData as MSD

from collections import defaultdict as ddict
from make_zones import make_zones
from make_lamps import make_lamps

@click.command()
def inputs():
	"""Prepares the executions inputs.
	"""

	print("Preparing the inputs for the experiment.")

	dir_name = "Inputs/"
	shutil.rmtree(dir_name,True)
	os.makedirs(dir_name)
	shutil.copy("inputs_params.in","Inputs/inputs_params.in")

	with open("inputs_params.in") as f:
		params = yaml.safe_load(f)

	if params['zones_inventory'] is not None and \
		params['lamps_inventory'] is not None:

		print("Validating the inventories.")

		lamps = np.loadtxt(params['lamps_inventory'],usecols=[0,1])
		zones = np.loadtxt(params['zones_inventory'],usecols=[0,1,2])
		zonData = pt.parse_inventory(params['zones_inventory'],7)

		hasLights = [ sum( x[0] for x in z ) != 0 for z in zonData ]

		circles = MSD.from_domain("domain.ini")
		for dat,b in zip(zones,hasLights):
			circles.set_circle((dat[0],dat[1]),dat[2]*1000,b)

		zones_ind = MSD.from_domain("domain.ini")
		for i,dat in enumerate(zones,1):
			zones_ind.set_circle((dat[0],dat[1]),dat[2]*1000,i)

		failed = set()
		for l,coords in enumerate(lamps,1):
			for i in range(len(circles)):
				try:
					col,row = circles._get_col_row(coords,i)
					if circles[i][row,col] and col >= 0 and row >= 0:
						zon_ind = zones_ind[i][row,col]
						failed.add((l,coords[0],coords[1],zon_ind))
				except IndexError:
					continue

		if len(failed):
			for l,lat,lon,zon_ind in sorted(failed):
				print("WARNING: Lamp #%d (%.06g,%.06g) falls within non-null zone #%d" \
					% (l,lat,lon,zon_ind))
			raise SystemExit()

	out_name = params['exp_name']

	print("Loading photometry files.")

	# Angular distribution (normalised to 1)
	lop_files = glob("Lights/*.lop")
	angles = np.arange(181,dtype=float)
	lop = {
		os.path.basename(s).rsplit('.',1)[0].split('_',1)[0] : \
		pt.load_lop(angles,s) for s in lop_files
	}

	# Spectral distribution (normalised with scotopric vision to 1)
	wav, viirs = np.loadtxt( "Lights/viirs.dat", skiprows=1 ).T
	viirs /= np.max(viirs)
	#viirs = pt.spct_norm(wav,viirs)
	# norm_spectrum = pt.load_spct(
	# 	wav,
	# 	np.ones(wav.shape),
	# 	"Lights/photopic.dat",
	# 	1
	# )
	wav, norm_spectrum = np.loadtxt("Lights/photopic.dat", skiprows=1).T
	norm_spectrum /= np.max(norm_spectrum)

	spct_files = glob("Lights/*.spct")
	spct = {
		os.path.basename(s).rsplit('.',1)[0].split('_',1)[0] : \
		pt.load_spct(wav,norm_spectrum,s) \
		for s in spct_files
	}

	print("Splitting in a few wavelengths.")

	n_bins = params['nb_bins']
	lmin = params['lambda_min']
	lmax = params['lambda_max']

	bool_array = (lmin<=wav)*(wav<lmax)

	limits = np.array(list(map(
		np.min, np.array_split(
			wav[bool_array],
			n_bins,
			-1 ) )) + [lmax] )

	lim_file = dir_name + "integration_limits.dat"
	with open(lim_file,'w') as f:
		f.write("%d\n"%n_bins)
	with open(lim_file,'ab') as f:
		np.savetxt(f,limits[:,np.newaxis])

	x = np.mean([limits[1:],limits[:-1]],0)

	print("Interpolating reflectance.")

	asper_files = glob("Lights/*.asper")
	asper = {
		os.path.basename(s).split('.',1)[0] : \
		np.loadtxt(s) \
		for s in asper_files
	}

	for type in asper:
		wl,refl = asper[type].T
		wl *= 1000.
		refl /= 100.
		asper[type] = np.interp(wav, wl, refl)

	sum_coeffs = sum(
		params['reflectance'][type] \
		for type in params['reflectance']
	)
	if sum_coeffs == 0:
		sum_coeffs = 1.

	refl = sum(
		asper[type]*coeff/sum_coeffs \
		for type,coeff \
		in params['reflectance'].items()
	)

	reflect = [ np.mean(a) for a in \
	    np.array_split(
			refl[bool_array],
			n_bins,
			-1
	    )
	]

	with open(dir_name+"/refl.lst",'w') as zfile:
		zfile.write('\n'.join( ["%.06g"%n for n in reflect])+'\n')

	print("Linking mie files.")

	aero_profile = params['aerosol_profile']
	RH = params['relative_humidity']
	mie_pre = aero_profile+"_RH%02d" % RH

	ppath = os.environ['PATH'].split(os.pathsep)
	illumpath = [s for s in ppath if "illumina" in s and "bin" not in s][0]

	mie_path = illumpath + "/Aerosol_optical_prop/"
	mie_files = glob(mie_path+mie_pre+"*.mie.out")
	mie_files = { int(s.split('.')[-3][:3]):s for s in mie_files }
	mie_wl = np.asarray(sorted(mie_files.keys()))
	wl2mie = np.asarray([min(mie_wl, key=lambda i: abs(i-j)) for j in x])

	for i in range(len(wl2mie)):
		name = dir_name+mie_pre.strip('_')+"_0.%03d0um.mie.out"%x[i]
		try:
			shutil.copy2(os.path.abspath(mie_files[wl2mie[i]]),name)
		except OSError as e:
			if e[0] != 17:
				raise

	shutil.copy("srtm.hdf5",dir_name)

	with open(dir_name+"/wav.lst",'w') as zfile:
		zfile.write('\n'.join( ["%03d"%n for n in x])+'\n')

	if params['zones_inventory'] is not None:
		dir_name = ".Inputs_zones/"
		inv_name = params['zones_inventory']
		n_inv = 7
		shutil.rmtree(dir_name,True)
		os.makedirs(dir_name)
		make_zones(dir_name,inv_name,n_inv,n_bins,params,out_name,
			x,lop,angles,wav,spct,viirs,refl,bool_array)

	if params['lamps_inventory'] is not None:
		dir_name = ".Inputs_lamps/"
		shutil.rmtree(dir_name,True)
		os.makedirs(dir_name)
		make_lamps(dir_name,n_bins,params,out_name,
			x,lop,angles,wav,spct,viirs,refl,bool_array)

	print("Unifying inputs.")

	lfiles = { fname.split(os.sep)[-1] for fname in glob(".Inputs_lamps/*") }
	zfiles = { fname.split(os.sep)[-1] for fname in glob(".Inputs_zones/*") }
	for fname in lfiles-zfiles:
		shutil.move(os.path.join(".Inputs_lamps",fname),"Inputs")
	for fname in zfiles-lfiles:
		shutil.move(os.path.join(".Inputs_zones",fname),"Inputs")
	for fname in zfiles&lfiles:
		if "fctem" in fname:
			shutil.move(os.path.join(".Inputs_lamps",fname),"Inputs")
		elif fname.endswith('.lst'):
			with open(os.path.join(".Inputs_lamps",fname)) as f:
				ldat = f.readlines()
			with open(os.path.join(".Inputs_zones",fname)) as f:
				zdat = f.readlines()
			with open(os.path.join("Inputs",fname),'w') as f:
				f.write(''.join(sorted(set(ldat+zdat))))
		elif fname.endswith('.hdf5'):
			ldat = MSD.Open(os.path.join(".Inputs_lamps",fname))
			zdat = MSD.Open(os.path.join(".Inputs_zones",fname))
			for i, dat in enumerate(ldat):
				zdat[i][dat != 0] = dat[dat != 0]
			zdat.save(os.path.join("Inputs",fname))
		else:
			print("WARNING: File %s not merged properly." % fname)
	if "origin.hdf5" not in zfiles:
		origin = MSD.from_domain("domain.ini")
		origin.save("Inputs/origin")
	shutil.rmtree(".Inputs_lamps",True)
	shutil.rmtree(".Inputs_zones",True)

	print("Done.")
