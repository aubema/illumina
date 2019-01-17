#!/usr/bin/env python

import pyproj
import yaml

print "Srcipt to automatically generate 'domain.ini'"
print "Defalt values are written within [brackets]"
print ""

scale_factor = 5
nb_layers = 7
scale_min = 1

domain = dict()

# Input observatory coordinates
obs_coords = raw_input("Enter observatory coordinates (lat,lon): ")
lat,lon = map(float,obs_coords.strip("()").split(","))

domain["observer"] = dict()
domain["observer"]["latitude"] = lat
domain["observer"]["longitude"] = lon

# Define projection
default_srs = "32" + ("6" if lat >= 0 else "7") + "%d" % (lon/6+31) # WGS84/UTM
srs = raw_input("\nEPSG projection to use [%s]: " % default_srs)
try:
    srs = str(int(srs))
except ValueError:
    print "  Using %s" % default_srs
    srs = default_srs

wgs84 = pyproj.Proj(init="epsg:4326")
proj  = pyproj.Proj(init="epsg:%s" % srs)

domain["projection"] = "epsg:%s" % srs

npix = scale_factor**2
x0,y0 = pyproj.transform(wgs84,proj,lon,lat)

domain["nb_pixels"] = npix
domain["observer"]["x"] = x0
domain["observer"]["y"] = y0

domain["extents"] = list()

for i in range(nb_layers):
    size = scale_min * scale_factor**i

    xmin = x0 - size*npix/2.
    xmax = x0 + size*npix/2.
    ymin = y0 - size*npix/2.
    ymax = y0 + size*npix/2.

    domain["extents"].append(dict())
    domain["extents"][-1]["pixel_size"] = size
    domain["extents"][-1]["bounding_box"] = dict()
    domain["extents"][-1]["bounding_box"]["xmin"] = xmin
    domain["extents"][-1]["bounding_box"]["xmax"] = xmax
    domain["extents"][-1]["bounding_box"]["ymin"] = ymin
    domain["extents"][-1]["bounding_box"]["ymax"] = ymax

with open("domain.ini",'w') as f:
    yaml.dump(domain,f,default_flow_style=False)

# print lon/lat bbox formatted for earthdata
SE = pyproj.transform(proj,wgs84,xmin,ymin)
SW = pyproj.transform(proj,wgs84,xmin,ymax)
NE = pyproj.transform(proj,wgs84,xmax,ymin)
NW = pyproj.transform(proj,wgs84,xmax,ymax)

N = max(NE[1],NW[1])
S = min(SE[1],SW[1])
E = max(NE[0],SE[0])
W = min(NW[0],SW[0])

print "\nBounding box:"
print "  SW: %f,%f" % (S,W)
print "  NE: %f,%f" % (N,E)
