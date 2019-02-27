#!/usr/bin/python2

# Script used to make polar plots of sky map
#
# Author: Alexandre Simoneau
# Date: 2017-02-08
#
# Usage: Call in a directory containing data extracted using
# 'extract-output-data.bash' with the '-s' option

import numpy as np
from glob import glob
from collections import defaultdict as ddict
import matplotlib.pyplot as plt

# Read data from files

data = ddict(lambda:ddict(lambda:ddict(list)))

for fname in glob("*.txt"):
    up = int(fname.split('.')[0].split('_')[-1])
    with open(fname) as f:
        d = f.readlines()
    for l in d:
        p, val = l.split()
        pos,scen,ta,wl,el,az,rd,sd = p.split('-')
        wl = int(wl[2:])
        el = float(el[2:])
        az = float(az[2:])
        val = float(val)
        
        data[pos][up][wl].append([el,az,val])

# Plot the data

for p in data:
    for u in data[p]:
        for w in data[p][u]:
            arr = np.asarray(sorted(data[p][u][w]))
            els = list(np.unique(arr[:,0]))
            azs = list(np.unique(arr[:,1]))
            d = np.zeros((len(els),len(azs)))

            for e,a,v in arr:
                d[els.index(e),azs.index(a)] = v
            
            # To have a full circle
            x = np.zeros((len(els),len(azs)+1))
            x[:,:-1] = d        # Copy data
            x[:,-1] = d[:,0]    # Full circle
            if els[-1] == 90:   # Uniform zenith
                x[-1,:] = d[-1][d[-1]!=0][0]
            x *= 1e6            # Normalisation

            theta = np.radians(azs+[azs[0]+360])
            r = 90 - np.asarray(els)

            # Control the interpolation resolution
            contour_levels = np.linspace(0,np.max(x),501)

#            data[p][u][w] = dict(el=els,az=azs,data=d)
            
            plt.figure()
            ax = plt.subplot(111,polar=True)
            ax.set_theta_zero_location('N')
            ax.xaxis.set_ticklabels(['N','NE','E','SE','S','SW','W','NW'])
            plt.contourf(theta,r,x,contour_levels)
#            plt.imshow(d,origin="lower",interpolation="None")
#            plt.xlabel("Azimutal angle")
#            plt.ylabel("Elevation angle")
            plt.title("Sky at %s for %d%% uplight at %dnm" % (p,u,w) )
            plt.colorbar()
            plt.savefig("%s_up%d_wl%d.png" % (p,u,w))
            plt.close()
