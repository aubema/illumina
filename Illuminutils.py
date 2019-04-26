#!/usr/bin/env python2

from glob import glob
import ogr, osr, gdal
import yaml, h5py
import hdftools, math
import os, tempfile, zipfile

def warp(srcfiles, projection, extent):
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

def rasterize(srcshape, projection, extent):
    width  = int(math.ceil(abs(extent["xmax"] - extent["xmin"]) / extent["pixel_size"]))
    height = int(math.ceil(abs(extent["ymax"] - extent["ymin"]) / extent["pixel_size"]))
    ds = gdal.GetDriverByName('GTiff').Create(os.path.join("WM.tiff"), width, height, 1, gdal.GDT_Byte)
    ds.SetGeoTransform([extent["xmin"], extent["pixel_size"], 0.0, extent["ymax"], 0.0, -extent["pixel_size"]])
    srs = osr.SpatialReference()
    srs.ImportFromProj4(projection)
    ds.SetProjection(srs.ExportToWkt())
    ds.GetRasterBand(1).Fill(0)
    gdal.RasterizeLayer(ds, [1], srcshape, burn_values=[1])
    return ds.GetRasterBand(1).ReadAsArray()

def save(params, data, dstname, scale_factor=1.):
    scaled_data = [ d*scale_factor for d in data ]
    ds = hdftools.from_domain(params,data)
    ds.save(dstname)

with open(glob("*.ini")[0]) as f:
    params = yaml.safe_load(f)

files = glob("SRTM/*.hgt")
if not len(files):
    print "ERROR: Could not find SRTM file(s), aborting."
    raise SystemExit
print "    ".join(map(str,files))
data = [ warp(files, params['srs'], extent) \
    for extent in params['extents'] ]
save(params, data, "srtm")

files = glob("VIIRS-DNB/*.tif")
if not len(files):
    print "WARNING: Did not find VIIRS file(s)."
    print "If you don't intend to use zones inventory, you cans safely ignore this."
else:
    print "    ".join(map(str,files))
    data = [ warp(files, params['srs'], extent) \
        for extent in params['extents'] ]
    save(params, data, "stable_lights")


    zip = zipfile.ZipFile('hydropolys.zip')
    tempdir = tempfile.mkdtemp(prefix='.',dir='.')

    print "Unzipping..."
    zip.extractall(tempdir)
    shape = ogr.Open(os.path.join(tempdir, 'hydropolys.shp')).GetLayer()

    data = [ rasterize(shape, params['srs'], extent) \
        for extent in params['extents'] ]
    save(params, data, "water_mask")

    # Clean up after ourselves.
    for file in os.listdir(tempdir):
      os.unlink(os.path.join(tempdir, file))
    os.rmdir(tempdir)
