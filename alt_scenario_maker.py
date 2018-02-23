#!/usr/bin/env python

import numpy as np
from pytools import *
from glob import glob
import os, sys

print "What type of scenario do you want to make?"
print "\t1. Different photometry at constant lumen"
print "\t2. Split lamp power based on spectral composition"

try:
    alt_type = int(raw_input())
except ValueError:
    sys.exit("\nUnrecognized input. Aborted.")

if alt_type not in range(1,3):
    sys.exit("\nInvalid type. Aborted.")

print ""

if alt_type == 1:
    new_name = raw_input("New name: ")
    old_invp = raw_input("Old inventory path: ")
    inv_path = raw_input("New inventory path: ")

if alt_type == 2:
    inv_path = raw_input("Inventory path: ")

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

#ratio_ps = float(raw_input("    photopic/scotopic ratio ? (0 <= p/(s+p) <= 1) : "))
ratio_ps = 1.
norm_spectrum = ratio_ps*photopic + (1-ratio_ps)*scotopic

spct_files = glob("Lights/*.spct")
spct = { os.path.basename(s).split('_',1)[0]:load_spct(wav,norm_spectrum,s) for s in spct_files }

if alt_type == 1:
    zonData = parse_inventory(inv_path)
    oldZonData = parse_inventory(old_invp,7)
if alt_type == 2:
    zonData = parse_inventory(inv_path,7)

print "\nCalculating the generalized lamps..."

# Calculate zones lamps
zones = make_zones( angles, lop, wav, spct, zonData )

# Make bins
x = np.loadtxt("Intrants/wav.lst")
n = x.size
ilims = np.genfromtxt("Intrants/integration_limits.dat",skip_header=1)
dl = ilims[1:]-ilims[:-1]
bool_array = (ilims[0]<=wav)*(wav<ilims[-1])
y = np.array(map(np.mean,np.array_split(zones[:,:,bool_array],n,-1),[-1]*n)).transpose(1,2,0)

# Photopic/scotopic spectrum
if alt_type == 1:
#   ratio_ps = float(raw_input("    photopic/scotopic ratio for lamp power ? (0 <= p/(s+p) <= 1) : "))
    nspct = ratio_ps*photopic + (1-ratio_ps)*scotopic
    nspct = nspct/np.sum(nspct)
    nspct = np.array(map(np.mean,np.array_split(nspct[bool_array],n)))

    dirnames = ["Intrants_"+new_name.replace(' ','_')+'/']
if alt_type == 2:
    lamp_types = set(map(lambda x: x[1], sum(zonData,[])))
    dirnames = [ "Intrants_"+lamp_t+'/' for lamp_t in lamp_types]

for dirname in dirnames:
    if not os.path.exists(dirname):
        os.makedirs(dirname)

# Link unmodified files
fctem = set(glob("Intrants/*fctem*"))
lumlp = set(glob("Intrants/*lumlp*"))
intrs = set(glob("Intrants/*"))

files = intrs-lumlp-fctem
for dirname in dirnames:
    for name in files:
        try:
            os.symlink(os.path.abspath(name), dirname+os.path.basename(name))
        except OSError as e:
            if e[0] != 17:
                raise

# Main treatement
if alt_type == 1:
    for z in xrange(len(zones)):
        print "Treating zone : %d/%d"%(z+1,len(zones))
        lum_files = sorted(glob("Intrants/*lumlp_%03d*.pgm"%(z+1)))
        spt_files = sorted(glob("Intrants/fctem*%03d*.dat"%(z+1)))

        # Linking files if inventory unchanged (faster)
        if zonData[z]==oldZonData[z]:
            for name in lum_files+spt_files:
                os.symlink(os.path.abspath(name), dirnames[0]+os.path.basename(name))
        else:
            head,p,data = load_pgm(lum_files[0])
            data = np.asarray([load_pgm(s)[2].reshape((1,-1)) for s in lum_files])
            lumtot = np.sum(data*(dl*nspct)[:,None,None],0) # Total luminosity
            spct = np.sum( y[z] * 2*np.pi*np.sin(np.deg2rad(angles))[:,None]
                         * (angles[1]-angles[0]), 0 )
            K = safe_divide(lumtot,np.sum(spct*nspct*dl))

            ndata = spct[:,None,None]*K[None]

            # Saving new files
            for wl in xrange(n):
                save_pgm( dirnames[0] + os.path.basename(lum_files[wl]), head, p, ndata[wl] )
                np.savetxt( dirnames[0] + os.path.basename(spt_files[wl]),
                            np.concatenate([y[z,:,wl],angles]).reshape((2,-1)).T )
if alt_type == 2:
    for lamp_t in lamp_types:
        print "Treating lamp : " + lamp_t
        dirname = "Intrants_"+lamp_t+'/'

        newZonData = [ filter( lambda z: z[1]==lamp_t, zone ) for zone in zonData ]
        newZonData[newZonData==[]] = [[0,lamp_t,'0']]
        weights = np.asarray(map(lambda z:sum(map(lambda x:x[0],z)),newZonData))

        newZones = make_zones( angles, lop, wav, spct, newZonData ) * weights[:,None,None]
        ny = np.array(map(np.mean,np.array_split(newZones[:,:,bool_array],n,-1),[-1]*n)).transpose(1,2,0)

        spct_old = np.sum( y  * 2*np.pi*np.sin(np.deg2rad(angles))[:,None]
                         * (angles[1]-angles[0]), 1 )
        spct_new = np.sum( ny * 2*np.pi*np.sin(np.deg2rad(angles))[:,None]
                         * (angles[1]-angles[0]), 1 )

        for z in xrange(len(zones)):
            print "  Treating zone : %d/%d"%(z+1,len(zones))
            lum_files = sorted(glob("Intrants/*lumlp_%03d*.pgm"%(z+1)))
            spt_files = sorted(glob("Intrants/fctem*%03d*.dat"%(z+1)))

            # Linking files if weigth = 1 (faster)
            if weights[z] == 1:
                for name in lum_files+spt_files:
                    os.symlink(os.path.abspath(name), dirname+os.path.basename(name))
            else:
                #head,p,data = load_pgm(lum_files[0])
                #data = np.asarray([load_pgm(s)[2] for s in lum_files])

                # Saving new files
                for wl in xrange(n):
                    head,p,data = load_pgm(lum_files[wl])
                    ndata = data * spct_new[z,wl]/spct_old[z,wl]

                    save_pgm( dirname + os.path.basename(lum_files[wl]), head, p, ndata )
                    np.savetxt( dirname + os.path.basename(spt_files[wl]),
                                np.concatenate([ny[z,:,wl],angles]).reshape((2,-1)).T )

print "Done."
