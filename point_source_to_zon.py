#!/usr/bin/env python3

import numpy as np
import pyproj, yaml
from collections import defaultdict as ddict

# Read domain.ini
dom_name = input("    Domain parameters filename : ")
with open(dom_name) as f:
    domain = yaml.safe_load(f)

# Read light inventory
# format:
#    lat lon pow hobs dobs fobs hlamp S_type P_type
inv_name = input("    LOP inventory filename : ")

inv = np.loadtxt(inv_name,dtype=object)
inv[:,:7] = inv[:,:7].astype(float)

# Classify point by pixel
data = ddict(list)

for light in inv:
    lat,lon = light[:2]

    p1 = pyproj.Proj("epsg:4326") # WGS84
    p2 = pyproj.Proj(domain['srs'])
    x, y = pyproj.transform(p1, p2, lon, lat, always_xy=True)

    for i in range(domain['nb_layers']):
        extent = domain['extents'][i]

        if extent['xmin'] < x < extent['xmax'] and \
           extent['ymin'] < y < extent['ymax']:
            level = extent['layer']
            col = int((x-extent['xmin'])/extent['pixel_size'] - 0.5)
            row = int((y-extent['ymin'])/extent['pixel_size'] - 0.5)
            break

    data[level,col,row].append(light[2:])

# Generate zone inventory
out_name = input("    Output filename : ")
with open(out_name,'w') as inv_file:
    inv_file.write("# X	Y	R	hobs	dobs	fobst   hlamp	Zone inventory\n")
    for level,col,row in data:
        lights = np.asarray(data[level,col,row])
        pow_tot = np.sum(lights[:,0])

        ho = np.average(lights[:,1],weights=lights[:,0])
        do = np.average(lights[:,2],weights=lights[:,0])
        fo = np.average(lights[:,3],weights=lights[:,0])
        hl = np.average(lights[:,4],weights=lights[:,0])

        frac = ddict(float)
        for i in range(len(lights)):
            frac[tuple(lights[i][-2:])] += 100*lights[i][0]/pow_tot

        extent = domain['extents'][domain['nb_layers']-level-1]
        x = extent['pixel_size'] * (col - 0.5) + extent['xmin']
        y = extent['pixel_size'] * (row - 0.5) + extent['ymin']

        lon,lat = pyproj.transform(p2,p1,x,y,always_xy=True)

        data_line = ("%.06f\t"*2+"%g\t"*5+"%s\n") % \
                    ( lat, lon,
                      extent['pixel_size']/2000., # diameter in km to radii in m
                      ho, do, fo, hl,
                      ' '.join( ["%g_%s_%s" %
                                (frac[i],i[0],i[1]) for i in frac] ) )

        inv_file.write(data_line)
