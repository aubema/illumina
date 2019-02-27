#!/usr/bin/env python2

import numpy as np
import pytools as pt
from glob import glob
import os, sys
from shutil import rmtree
import MultiScaleData as MSD
import argparse

parser = argparse.ArgumentParser(
    description="Generates alternatives scenarios based on the current one."
)
parser.add_argument("mode",
    help="Execution mode."\
    "`conversion` replace the current light inventory by a new one "\
    "while keeping the total emited lumen constant. "\
    "`split` separates the inventory by lighting technology. ",
    choices=['conversion','split']
)
parser.add_argument("name", nargs='?',
    help="Name of the new scenario if the mode is set to `conversion`. "\
    "The new inventory data needs to be in `inventory_NAME.txt`. "\
    "If the mode is set to `split`, only the specified technology will be extracted. "\
    "If that technology is not present in the inventory, the output will be all zeros."
)

p = parser.parse_args()

if p.mode == "conversion":
    if p.name is None:
        parser.error("conversion requires name")
    new_name = p.name
    old_invp = "inventory.txt"
    inv_path = "inventory_%s.txt" % new_name
elif p.mode == "split":
    inv_path = "inventory.txt"

print "\nLoading data..."

# Angular distribution (normalised to 1)
lop_files = glob("Lights/*.lop")
angles = np.arange(181,dtype=float)
lop = {
    os.path.basename(s).split('_',1)[0] : \
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
    os.path.basename(s).split('_',1)[0] : \
    pt.load_spct(wav,norm_spectrum,s) \
for s in spct_files }

if p.mode == "conversion":
    zonData = pt.parse_inventory(inv_path)
    oldZonData = pt.parse_inventory(old_invp,7)
elif p.mode == "split":
    zonData = pt.parse_inventory(inv_path,7)

print "\nCalculating the generalized lamps..."

# Calculate zones lamps
zones = pt.make_zones( angles, lop, wav, spct, zonData )

# Make bins
x = np.loadtxt("Inputs/wav.lst")
n = x.size
ilims = np.genfromtxt("Inputs/integration_limits.dat",skip_header=1)
dl = ilims[1:]-ilims[:-1]
bool_array = (ilims[0]<=wav)*(wav<ilims[-1])
y = np.array(map(np.mean,
        np.array_split(zones[:,:,bool_array],n,-1),[-1]*n
    )).transpose(1,2,0)

# Photopic/scotopic spectrum
if p.mode == "conversion":
#   ratio_ps = float(raw_input("    photopic/scotopic ratio for lamp power ? (0 <= p/(s+p) <= 1) : "))
    nspct = ratio_ps*photopic + (1-ratio_ps)*scotopic
    nspct = nspct/np.sum(nspct)
    nspct = np.array(map(np.mean,np.array_split(nspct[bool_array],n)))

    dirnames = ["Inputs_"+new_name.replace(' ','_')+'/']
elif p.mode == "split":
    lamp_types = set(map(lambda x: x[1], sum(zonData,[])))
    if p.name is not None:
        if p.name not in lamp_types:
            print "WARNING: Technology %s not present in inventory." % p.name
        lamp_types = [ p.name ]
    dirnames = [ "Inputs_"+lamp_t+'/' for lamp_t in lamp_types]

for dirname in dirnames:
    if os.path.exists(dirname):
        rmtree(dirname)
    os.makedirs(dirname)

# Link unmodified files
fctem = set(glob("Inputs/*fctem*"))
lumlp = set(glob("Inputs/*lumlp*"))
intrs = set(glob("Inputs/*"))

files = intrs-lumlp-fctem
for dirname in dirnames:
    for name in files:
        os.symlink(
            os.path.abspath(name),
            dirname+os.path.basename(name)
        )

# Main treatement
if p.mode == "conversion":
    for z in xrange(len(zones)):
        print "Treating zone : %d/%d"%(z+1,len(zones))
        lum_files = sorted(glob("Inputs/*lumlp_%03d*.hdf5"%(z+1)))
        spt_files = sorted(glob("Inputs/fctem*%03d*.dat"%(z+1)))

        # Linking files if inventory unchanged (faster)
        if zonData[z]==oldZonData[z]:
            for name in lum_files+spt_files:
                os.symlink(
                    os.path.abspath(name),
                    dirnames[0]+os.path.basename(name)
                )
        else:
            attrs = MSD.Open(lum_files[0])._attrs
            data = np.asarray([MSD.Open(s) for s in lum_files])
            lumtot = np.sum(data*(dl*nspct)[:,None],0) # Total luminosity
            spct = np.sum( y[z] * 2*np.pi*np.sin(np.deg2rad(angles))[:,None]
                * (angles[1]-angles[0]), 0 )
            K = pt.safe_divide(lumtot,np.sum(spct*nspct*dl))

            ndata = spct[:,None]*K[None]

            # Saving new files
            for wl,dat in enumerate(ndata):
                dat = MSD.MultiScaleData(attrs,dat)
                dat.save( dirnames[0] + os.path.basename(lum_files[wl]) )
                np.savetxt( dirnames[0] + os.path.basename(spt_files[wl]),
                    np.concatenate([y[z,:,wl],angles]).reshape((2,-1)).T )
elif p.mode == "split":
    for lamp_t in lamp_types:
        print "Treating lamp : " + lamp_t
        dirname = "Inputs_"+lamp_t+'/'

        newZonData = [ [l for l in zone if l[1]==lamp_t] for zone in zonData ]
        for i in range(len(newZonData)):
            if newZonData[i]==[]:
                newZonData[i] = [[0,lamp_t,'0']]
        weights = np.asarray(map(lambda z:sum(map(lambda x:x[0],z)),newZonData))

        newZones = pt.make_zones( angles, lop, wav, spct, newZonData ) \
            * weights[:,None,None]
        ny = np.array(map(np.mean,
                np.array_split(newZones[:,:,bool_array],n,-1),[-1]*n
            )).transpose(1,2,0)

        spct_old = np.sum( y  * 2*np.pi*np.sin(np.deg2rad(angles))[:,None]
            * (angles[1]-angles[0]), 1 )
        spct_new = np.sum( ny * 2*np.pi*np.sin(np.deg2rad(angles))[:,None]
            * (angles[1]-angles[0]), 1 )

        for z in xrange(len(zones)):
            print "  Treating zone : %d/%d"%(z+1,len(zones))
            lum_files = sorted(glob("Inputs/*lumlp_%03d*.hdf5"%(z+1)))
            spt_files = sorted(glob("Inputs/fctem*%03d*.dat"%(z+1)))

            # Linking files if weigth = 1 (faster)
            if weights[z] == 1:
                for name in lum_files+spt_files:
                    os.symlink(
                        os.path.abspath(name),
                        dirname+os.path.basename(name)
                    )
            else:
                # Saving new files
                for wl in xrange(n):
                    data = MSD.Open(lum_files[wl])
                    ndata = data * spct_new[z,wl]/spct_old[z,wl]

                    ndata.save( dirname + os.path.basename(lum_files[wl]) )
                    np.savetxt( dirname + os.path.basename(spt_files[wl]),
                        np.concatenate([ny[z,:,wl],angles]).reshape((2,-1)).T )

print "Done."
