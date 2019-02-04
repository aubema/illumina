#!/usr/bin/env python

from glob import glob
import gdal, yaml, h5py
import hdftools

def warp(srcfiles, projection=None, extent=None):
    bounding_box = [
        extent["xmin"],
        extent["ymin"],
        extent["xmax"],
        extent["ymax"] ]

    vrt = gdal.BuildVRT('',srcfiles)
    ds = gdal.Warp( '', vrt,
        format="VRT",
        dstSRS=projection,
        dstNodata=0.,
        outputBounds=bounding_box,
        xRes=extent["pixel_size"],
        yRes=extent["pixel_size"],
        resampleAlg="cubicspline" )

    return ds.GetRasterBand(1).ReadAsArray()

def MYD09A1_band_name(fname, band_n):
    sub_ds = gdal.Open(fname).GetSubDatasets()
    return next( s[0] for s in sub_ds if ("b%02d" % band_n) in s[0] )

def save(params, data, dstname, scale_factor=1.):
    scaled_data = [ d*scale_factor for d in data ]
    ds = hdftools.from_domain(params,data)
    ds.save(dstname)

with open(glob("*.ini")[0]) as f:
    params = yaml.load(f)

files = glob("SRTM/*.hgt")
print "    ".join(map(str,files))
data = [ warp(files, params['srs'], extent) \
    for extent in params['extents'] ]
save(params, data, "srtm")

for band_n in xrange(1,8):
    files = glob("MODIS/*.hdf")
    fname = "refl_b%02d" % band_n
    band_names = map( lambda f: MYD09A1_band_name(f, band_n), files )
    print "    ".join(map(str,band_names))
    data = [ warp(band_names, params['srs'], extent) \
        for extent in params['extents'] ]
    save(params, data, fname, scale_factor=0.0001)

files = glob("VIIRS-DNB/*.tif")
print "    ".join(map(str,files))
data = [ warp(files, params['srs'], extent) \
    for extent in params['extents'] ]
save(params, data, "stable_lights")
