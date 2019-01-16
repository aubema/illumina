#!/usr/bin/env python

import pyproj

print "Srcipt to automatically generate 'domain.ini'"
print "Defalt values are written within [brackets]"
print ""

with open("domain.ini",'w') as f:
    f.write("---\n") # Begin YAML file

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

wgs84 = pyproj.Proj(init="epsg:4326")
proj  = pyproj.Proj(init="epsg:%s" % srs)

with open("domain.ini",'a') as f:
    f.write("projection: epsg:%s\n" % srs)
    f.write("extents:\n")

npix = 25
x0,y0 = pyproj.transform(wgs84,proj,lon,lat)
for i in range(7):
    size = 5**i

    xmin = x0 - size*npix/2.
    xmax = x0 + size*npix/2.
    ymin = y0 - size*npix/2.
    ymax = y0 + size*npix/2.

    with open("domain.ini",'a') as f:
        f.write("  - pixsize: %d\n" % size)
        f.write("    bounding_box:\n")
        f.write("      xmin: %d\n" % xmin)
        f.write("      xmax: %d\n" % xmax)
        f.write("      ymin: %d\n" % ymin)
        f.write("      ymax: %d\n" % ymax)


with open("domain.ini",'a') as f:
    f.write("...\n") # End YAML file

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
