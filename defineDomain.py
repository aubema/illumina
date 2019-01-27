#!/usr/bin/env python

import pyproj
import yaml
import math

def eng_format(x, unit=''):
    # Credit: 200_success on StackOverflow
    # https://codereview.stackexchange.com/a/50971
    #
    # U+03BC is Greek lowercase mu
    UNITS = [' ', ' k', ' M', ' G'] + \
            ([None] * 10) + \
            [' f', ' p', ' n', u' \u03bc', ' m']

    power_of_1000 = int(math.floor(math.log10(x) // 3))
    exponent = 3 * power_of_1000
    prefix = UNITS[power_of_1000]
    if prefix is None:
        prefix = '*10^%d ' % exponent

    significand = x * 10**(-exponent)
    return '%.2f%s%s' % (significand, prefix, unit)

with open("domain_params.in") as f:
    domain = yaml.load(f)

lat = domain.pop("latitude")
lon = domain.pop("longitude")

domain["observer"] = dict(latitude=lat, longitude=lon)

# Define projection
if domain["srs"] == "auto":
    default_srs = \
        "epsg:32" + \
        ("6" if lat >= 0 else "7") + \
        "%d" % (lon/6+31) # WGS84/UTM
    domain["srs"] = default_srs

wgs84 = pyproj.Proj(init="epsg:4326")
proj  = pyproj.Proj(init=domain['srs'])

npix = domain["scale_factor"]**2
x0,y0 = pyproj.transform(wgs84,proj,lon,lat)

domain["nb_pixels"] = npix

domain["extents"] = list()

for i in range(domain["nb_layers"]):
    size = domain["scale_min"] * domain["scale_factor"]**i

    print "Layer",i
    print "Pixel size:", eng_format(size,'m')
    print "Domain size:", eng_format(size*npix,'m')
    print ""

    xmin = x0 - size*npix/2.
    xmax = x0 + size*npix/2.
    ymin = y0 - size*npix/2.
    ymax = y0 + size*npix/2.

    extent = dict()
    extent["layer"] = i
    extent["pixel_size"] = size
    bbox = dict(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
    extent.update(bbox)
    domain["extents"].append(extent)

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

print "Bounding box:"
print "  SW: %f,%f" % (S,W)
print "  NE: %f,%f" % (N,E)
