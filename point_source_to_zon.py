#!/usr/bin/env python

import numpy as np
import pyproj
from collections import defaultdict as ddict

# Read domain.ini
dom_name = raw_input("    Domain parameters filename : ")
def params_parser(fname):
    with open(fname) as f:
        lines = filter(lambda line: ':' in line, f.readlines())

    params = map(lambda s: map(str.strip, s.split(':',1)), lines)
    return { p[0]:p[1] for p in params }

domain = params_parser(dom_name)

domain['pixsize'] = float(domain['pixsize'])

domain['xmin'],domain['ymin'],domain['xmax'],domain['ymax'] = map(float,domain['bbox'].split())

domain['xsize'] = int((domain['xmax']-domain['xmin'])/domain['pixsize'])
domain['ysize'] = int((domain['ymax']-domain['ymin'])/domain['pixsize'])

# Read light inventory
# format:
#    lat lon pow hobs dobs fobs hlamp S_type P_type
inv_name = raw_input("    LOP inventory filename : ")

inv = np.loadtxt(inv_name,dtype=object)
inv[:,:7] = inv[:,:7].astype(float)

# Classify point by pixel
data = ddict(list)

for light in inv:
    lat,lon = light[:2]

    p1 = pyproj.Proj(init="epsg:4326") # WGS84
    p2 = pyproj.Proj(init=domain['srs'])
    x, y = pyproj.transform(p1, p2, lon, lat)

    col = int((x-domain['xmin'])/domain['pixsize']) + 1
    row = int((y-domain['ymin'])/domain['pixsize']) + 1

    data[col,row].append(light[2:])

# Generate zone inventory
out_name = raw_input("    Output filename : ")
with open(out_name,'w') as inv_file:
    inv_file.write("# X	Y	R	hobs	dobs	fobst   hlamp	Zone inventory		# Comment\n")
    for col,row in data:
        lights = np.asarray(data[col,row])
        pow_tot = np.sum(lights[:,0])

        ho = np.average(lights[:,1],weights=lights[:,0])
        do = np.average(lights[:,2],weights=lights[:,0])
        fo = np.average(lights[:,3],weights=lights[:,0])
        hl = np.average(lights[:,4],weights=lights[:,0])

        frac = ddict(float)
        for i in range(len(data[col,row])):
            frac[tuple(data[col,row][i][-2:])] += 100*data[col,row][i][0]/pow_tot

        x = domain['pixsize'] * (col - 0.5) + domain['xmin']
        y = domain['pixsize'] * (row - 0.5) + domain['ymin']

        lon,lat = pyproj.transform(p2,p1,x,y)

        data_line = ("%g\t"*7+"%s\n") % \
                    ( lat, lon,
                      domain['pixsize']/2/1000,
                      ho, do, fo, hl,
                      ' '.join( map(lambda i: "%g_%s_%s" %
                                (frac[i],i[0],i[1]), frac) ) )

        inv_file.write(data_line)
