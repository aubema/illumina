#!/usr/bin/env python3

import numpy as np
import pytools as pt
from glob import glob
import os, sys
import shutil
import MultiScaleData as MSD
import yaml, click

from scipy.interpolate import interp1d as interp
from collections import defaultdict as ddict

@click.command()
@click.argument("name")
@click.option("-z","--zones",type=click.Path(exists=True),help="New zones inventory filename.")
@click.option("-l","--lights",type=click.Path(exists=True),help="New discrete lights inventory filename.")
def alternate(name,zones,lights):
    """Generates an alternate scenario at constant lumen.

    This scenatio will be based on the content of the `Inputs` folder and
    will be placed in a folder named `Inputs_NAME`.
    """
    if zones == None and lights == None:
        print("ERROR: At least one of 'zones' and 'lights' must be provided.")
        raise SystemExit

    dirname = "Inputs_%s/" % name

    if os.path.exists(dirname):
        shutil.rmtree(dirname)
    os.makedirs(dirname)

    with open("inputs_params.in") as f:
        params = yaml.safe_load(f)

    if zones is not None and \
    	lights is not None:

    	print("Validating the inventories.")

    	lamps = np.loadtxt(lights,usecols=[0,1])
    	zones = np.loadtxt(params['zones_inventory'],usecols=[0,1,2])
    	zonData = pt.parse_inventory(zones,0)

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


    print("\nLoading data...")

    # Angular distribution (normalised to 1)
    lop_files = glob("Lights/*.lop")
    angles = np.arange(181,dtype=float)
    lop = {
        os.path.basename(s).rsplit('.',1)[0].split('_',1)[0] : \
        pt.load_lop(angles,s) \
    for s in lop_files }

    # Spectral distribution (normalised with scotopric vision to 1)
    wav, viirs = np.loadtxt("Lights/viirs.dat", skiprows=1).T
    viirs = pt.spct_norm(wav,viirs)
    scotopic = pt.load_spct(wav, np.ones(wav.shape), "Lights/scotopic.dat", 1)
    photopic = pt.load_spct(wav, np.ones(wav.shape), "Lights/photopic.dat", 1)

    #ratio_ps = float(raw_input("    photopic/scotopic ratio ? (0 <= p/(s+p) <= 1) : "))
    ratio_ps = 1.
    norm_spectrum = ratio_ps*photopic + (1-ratio_ps)*scotopic

    spct_files = glob("Lights/*.spct")
    spct = {
        os.path.basename(s).rsplit('.',1)[0].split('_',1)[0] : \
        pt.load_spct(wav,norm_spectrum,s) \
    for s in spct_files }

    # Make bins
    x = np.loadtxt("Inputs/wav.lst",ndmin=1).tolist()
    n_bins = len(x)
    shutil.copy("Inputs/integration_limits.dat",dirname)
    ilims = np.genfromtxt("Inputs/integration_limits.dat",skip_header=1)
    dl = ilims[1:]-ilims[:-1]
    bool_array = (ilims[0]<=wav)*(wav<ilims[-1])

    ppath = os.environ['PATH'].split(os.pathsep)
    illumpath = [s for s in ppath if "illumina" in s and "bin" not in s][0]

    out_name = params['exp_name']

    asper_files = glob("Lights/*.aster")
    asper = {
    	os.path.basename(s).split('.',1)[0] : \
    	np.loadtxt(s) \
    	for s in asper_files
    }

    for type in asper:
    	wl,refl = asper[type].T
    	wl *= 1000.
    	refl /= 100.
    	asper[type] = interp(
    		wl, refl,
    		bounds_error=False,
    		fill_value=0.
    	)(wav)

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

    with open(dirname+"/refl.lst",'w') as zfile:
    	zfile.write('\n'.join( ["%.06g"%n for n in reflect])+'\n')

    # Photopic/scotopic spectrum
    #   ratio_ps = float(raw_input("    photopic/scotopic ratio for lamp power ? (0 <= p/(s+p) <= 1) : "))
    nspct = ratio_ps*photopic + (1-ratio_ps)*scotopic
    nspct = nspct/np.sum(nspct)
    nspct = np.array(list(map(np.mean,np.array_split(nspct[bool_array],n_bins))))

    for aero_file in glob("Inputs/*.txt"):
        shutil.copy(aero_file,aero_file.replace("Inputs",dirname))

    # print("Linking mie files.")
    #
    # aero_profile = params['aerosol_profile']
    # RH = params['relative_humidity']
    # mie_pre = aero_profile+"_RH%02d" % RH
    #
    # ppath = os.environ['PATH'].split(os.pathsep)
    # illumpath = [s for s in ppath if "illumina" in s and "bin" not in s][0]
    #
    # mie_path = illumpath + "/Aerosol_optics/"
    # mie_files = glob(mie_path+mie_pre+"*.mie.out")
    # mie_files = { int(s.split('.')[-3][:3]):s for s in mie_files }
    # mie_wl = np.asarray(sorted(mie_files.keys()))
    # wl2mie = np.asarray([min(mie_wl, key=lambda i: abs(i-j)) for j in x])
    #
    # for i in range(len(wl2mie)):
    # 	name = dirname+mie_pre.strip('_')+"_0.%03d0um.mie.out"%x[i]
    # 	try:
    # 		shutil.copy2(os.path.abspath(mie_files[wl2mie[i]]),name)
    # 	except OSError as e:
    # 		if e[0] != 17:
    # 			raise

    shutil.copy("srtm.hdf5",dirname)

    with open(dirname+"/wav.lst",'w') as zfile:
    	zfile.write('\n'.join( ["%03d"%n for n in x])+'\n')

    if zones is not None:
        dir_name = ".Inputs_zones/"
        inv_name = zones
        n_inv = 0
        shutil.rmtree(dir_name,True)
        os.makedirs(dir_name)
        exec(compile(open(os.path.join(illumpath,"make_zones.py"), "rb").read(), os.path.join(illumpath,"make_zones.py"), 'exec'))

        oldlumlp = MSD.from_domain("domain.ini")
        for fname in glob("Inputs/*lumlp*"):
            ds = MSD.Open(fname)
            wl = int(fname.split('_')[1])
            for i, dat in enumerate(ds):
                oldlumlp[i] += dat * nspct[x.index(wl)] * dl[x.index(wl)]

        newlumlp = MSD.from_domain("domain.ini")
        for fname in glob(os.path.join(dir_name,"*lumlp*")):
            ds = MSD.Open(fname)
            wl = int(fname.split('_')[2])
            for i, dat in enumerate(ds):
                newlumlp[i] += dat * nspct[x.index(wl)] * dl[x.index(wl)]

        ratio = MSD.from_domain("domain.ini")
        for i in range(len(ratio)):
            ratio[i] = pt.safe_divide(oldlumlp[i],newlumlp[i])

        for fname in glob(os.path.join(dir_name,"*lumlp*")):
            ds = MSD.Open(fname)
            for i, dat in enumerate(ratio):
                ds[i] *= dat
            ds.save(fname)

    if lights is not None:
        params['lamps_inventory'] = lights
        dir_name = ".Inputs_lamps/"
        shutil.rmtree(dir_name,True)
        os.makedirs(dir_name)
        exec(compile(open(os.path.join(illumpath,"make_lamps.py"), "rb").read(), os.path.join(illumpath,"make_lamps.py"), 'exec'))

    print("Unifying inputs.")

    lfiles = { fname.split(os.sep)[-1] for fname in glob(".Inputs_lamps/*") }
    zfiles = { fname.split(os.sep)[-1] for fname in glob(".Inputs_zones/*") }
    for fname in lfiles-zfiles:
    	shutil.move(os.path.join(".Inputs_lamps",fname),dirname)
    for fname in zfiles-lfiles:
    	shutil.move(os.path.join(".Inputs_zones",fname),dirname)
    for fname in zfiles&lfiles:
    	if "fctem" in fname:
    		shutil.move(os.path.join(".Inputs_lamps",fname),dirname)
    	elif fname.endswith('.lst'):
    		with open(os.path.join(".Inputs_lamps",fname)) as f:
    			ldat = f.readlines()
    		with open(os.path.join(".Inputs_zones",fname)) as f:
    			zdat = f.readlines()
    		with open(os.path.join(dirname,fname),'w') as f:
    			f.write(''.join(sorted(set(ldat+zdat))))
    	elif fname.endswith('.hdf5'):
    		ldat = MSD.Open(os.path.join(".Inputs_lamps",fname))
    		zdat = MSD.Open(os.path.join(".Inputs_zones",fname))
    		for i, dat in enumerate(ldat):
    			zdat[i][dat != 0] = dat[dat != 0]
    		zdat.save(os.path.join(dirname,fname))
    	else:
    		print("WARNING: File %s not merged properly." % fname)
    shutil.rmtree(".Inputs_lamps",True)
    shutil.rmtree(".Inputs_zones",True)

    print("Done.")
