#!/usr/bin/env python

import numpy as np, matplotlib.pyplot as plt, re, os
import subprocess as sub, shutil
from collections import OrderedDict as odict, defaultdict as ddict
from glob import glob
from pytools import *
	
# Load data
# IMPORTANT : x axis of all similar data must be the same

shutil.rmtree("Intrants",True)
shutil.rmtree("Lamps",True)
print "Loading data..."

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

# lamps distribution
inv_name = raw_input("    LOP inventory filename : ")
zonData = parse_inventory(inv_name,7)

# make zon file
with open(inv_name) as f:
	zonfile = strip_comments(f.readlines())
zonfile = map(lambda s: s.split()[:7], zonfile)

print "Calculating the generalized lamps..."

# Calculate zones lamps
zones = make_zones(angles, lop, wav, spct, zonData )

print "Saving data..."

dirname = "Lamps"
if not os.path.exists(dirname):
	os.makedirs(dirname)

for i in xrange(len(zones)):
	bin_zon = np.zeros((len(angles)+1,len(wav)+1))
	bin_zon[0,0] = len(wav+1)
	bin_zon[1:,0] = angles
	bin_zon[0,1:] = wav
	bin_zon[1:,1:] = zones[i]

	np.savetxt("Lamps/zone%i_lamp.dat"%(i+1), bin_zon[1:])
	bin_zon.astype(np.float32).tofile("Lamps/zone%i_lamp.bin"%(i+1))

print "Plotting..."

sub.call(["zones_plot.sh","%d"%len(zones)])

print "Splitting in a few wavelengths..."
n = int(raw_input("    Number of wavelengths to use : "))
lmin = float(raw_input("    lambda min : "))
lmax = float(raw_input("    lambda max : "))

bool_array = (lmin<=wav)*(wav<lmax)

limits = np.array(map(np.min,np.array_split(wav[bool_array],n,-1))+[lmax])

filename = "integration_limits.dat"
with open(filename,'w') as f:
	f.write("%d\n"%n)
with open(filename,'ab') as f:
	np.savetxt(f,limits[:,np.newaxis])

# Create the desired lamp files
x = np.mean([limits[1:],limits[:-1]],0)
y = np.array(map(np.mean,np.array_split(zones[:,:,bool_array],n,-1),[-1]*n)).transpose(1,2,0)

print "Creating files..."

for l in xrange(n):
	dirname = "Intrants/"
	if not os.path.exists(dirname):
		os.makedirs(dirname)
	for z in xrange(len(zones)):
		np.savetxt( dirname+"fctem_wl_%03d_zon_%03d.dat"%(x[l],z+1), np.concatenate([ y[z,:,l],angles ]).reshape((2,-1)).T )

ans = raw_input("    Preparing files for viirs2lum ? ([y]/n) ")

stop = False
try:
	if ans[0] in ['N','n']:
		stop = True
except IndexError:
	pass

if not stop:
	# Questions for viirs2lum
	out_name = raw_input("Output root name of the experiment [this name will be used for all the subsequent files]?\n    ")
	pgm_name = raw_input("viirs-dnb file name? [e.g. stable_lights.pgm]\n    ")
	modis_name = raw_input("modis reflectance file list file name? [e.g. modis.dat]\n    ")
	modis_dir = raw_input("modis directory? [e.g. pgms]\n    ")
	zon_name = out_name+".zon"
	srtm_name = raw_input("elevation file name? [e.g. srtm.pgm]\n    ")
	
	dir_name = "./Intrants/"
	tmp_names = {'sat':"stable_lights.pgm", 'modis':"modis.dat", 'viirs':"viirs.dat", 'zon':"zone.zon"}
	tmp_list = set(tmp_names.values()) # List of files to be removed after execution
	tmp_names['integ'] = "integration_limits.dat"
	tmp_names['srtm'] = "srtm.pgm"

        ans = raw_input("Cutoff low values in the viirs-dnb data? ([y]/n) ")
        viirs_lowcut = True
        try:
            if ans[0] in ['N','n']:
                viirs_lowcut = False
        except IndexError:
            pass

        viirs_head,viirs_p,viirs_dat = load_pgm(pgm_name)

        if viirs_lowcut:
            try:
                val = float(raw_input("    Cutoff value [default peak value/4] : "))
            except ValueError:
                v,c = np.unique(viirs_dat,return_counts=True)
                val = v[c>=c.max()/4][-1]

            viirs_reject = viirs_dat.copy()
            viirs_reject[viirs_reject >= val] = 0
            viirs_dat[viirs_dat < val] = 0
            
            save_pgm(dir_name+"NaturalAndScatteredLight.pgm", viirs_head, viirs_p, viirs_reject)

        save_pgm(dir_name+tmp_names['sat'], viirs_head, viirs_p, viirs_dat)

    	# Creating symbolic links to necessary reflectance and photometry files
	modis_files = np.genfromtxt(modis_name,skip_header=1,usecols=1,dtype=str)
	modis_files = map(lambda s: modis_dir+"/"+s,modis_files)
	zon_files = [ "Lamps/zone%d_lamp.dat" % (i+1) for i in xrange(len(zones)) ]
	for filename in np.concatenate([modis_files,zon_files]):
		name = dir_name+os.path.basename(filename)
		try:
			os.symlink(os.path.abspath(filename),name)
		except OSError as e:
			if e[0] != 17:
				raise
		tmp_list.add(name)

    	# Creating the zon file based on the inventory
	with open(zon_name,'w') as f:
		f.write("%d\n" % len(zonfile))
		for i in xrange(len(zonfile)):
			zonfile[i].insert(3,os.path.basename(zon_files[i]))
			f.write((("%s\t"*len(zonfile[i]))[:-1]+"\n") % tuple(zonfile[i]))

    	# Linking useful files
#	os.symlink(os.path.abspath(pgm_name),dir_name+tmp_names['sat'])
	os.symlink(os.path.abspath(modis_name),dir_name+tmp_names['modis'])
	os.symlink(os.path.abspath("Lights/viirs.dat"),dir_name+tmp_names['viirs'])
	os.symlink(os.path.abspath(zon_name),dir_name+tmp_names['zon'])
	os.symlink(os.path.abspath("integration_limits.dat"),dir_name+tmp_names['integ'])
	os.symlink(os.path.abspath(srtm_name),dir_name+tmp_names['srtm'])

	print "Linking mie files..."
	
	mie_pre = raw_input("    Mie file prefix : ")
	illum_dir = sorted([s for s in os.environ['PATH'].split(':') if 'illumina' in s.lower() ], key=lambda s:len(s))[0]
	mie_path = illum_dir + "/Aerosol_optical_prop/"
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
#		tmp_list.add(name)
	
ans = raw_input("    Executing viirs2lum ? ([y]/n) ")

stop = False
try:
	if ans[0] in ['N','n']:
		stop = True
except IndexError:
	pass
	
if not stop:
	print "Launching Fortran..."

	os.chdir(dir_name)
	p = sub.Popen("viirs2lum", stdin=sub.PIPE)
	param = out_name+"\n"+os.path.basename(tmp_names['sat'])+"\n"+os.path.basename(tmp_names['zon'])+"\n"
	p.communicate(param)
	
	print "Fortran done."
	
	for filename in tmp_list:
		os.remove(os.path.basename(filename))

	#mie_files = sorted(glob("*.mie.out"))
	# Writing both .lst files
	with open("zon.lst",'w') as zfile:
		zfile.write('\n'.join( map(lambda n:"%03d"%n, xrange(1,len(zones)+1) ))+'\n')
	with open("wav.lst",'w') as zfile:
		zfile.write('\n'.join( map(lambda n:"%03d"%n, x ))+'\n')
  
	os.chdir("..")

print "Done."


