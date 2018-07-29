#!/usr/bin/env python

import pyproj

print "Srcipt to automatically generate 'domain.ini'"
print "Defalt values are written within [brackets]"
print ""

# Input observatory coordinates
obs_coords = raw_input("Enter observatory coordinates (lat,lon): ")
lat,lon = map(float,obs_coords.strip("()").split(","))

# Define projection
default_srs = "32" + ("6" if lat >= 0 else "7") + "%d" % (lon/6+31) # WGS84/UTM
srs = raw_input("\nEPSG projection to use [%s]: " % default_srs)
try:
    srs = str(int(srs))
except ValueError:
    print "  Using %s" % default_srs
    srs = default_srs

# Define resolution
default_pixsize = 1000
pixsize = raw_input("Pixel resolution (in meters) [%d]: " % default_pixsize)
try:
    pixsize = int(pixsize)
except ValueError:
    print "  Using %d" % default_pixsize
    pixsize = default_pixsize

# Define domain size
default_xsize = 200
xsize = raw_input("Domain x radius (km) [%d]: " % default_xsize)
try:
    xsize = int(xsize)
except ValueError:
    print "  Using %d" % default_xsize
    xsize = default_xsize

default_ysize = 200
ysize = raw_input("Domain y radius (km) [%d]: " % default_ysize)
try:
    ysize = int(ysize)
except ValueError:
    print "  Using %d" % default_ysize
    ysize = default_ysize

# Create domain.ini file
wgs84 = pyproj.Proj(init="epsg:4326")
proj  = pyproj.Proj(init="epsg:%s" % srs)

x0,y0 = pyproj.transform(wgs84,proj,lon,lat)
xmin = x0 - xsize*1000 - pixsize/2
xmax = x0 + xsize*1000 + pixsize/2
ymin = y0 - ysize*1000 - pixsize/2
ymax = y0 + ysize*1000 + pixsize/2

with open("domain.ini",'w') as f:
    f.write("bbox: %d %d %d %d\n" % (xmin, ymin, xmax, ymax))
    f.write("pixsize: %d\n" % pixsize)
    f.write("srs: epsg:%s\n" % srs)

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
