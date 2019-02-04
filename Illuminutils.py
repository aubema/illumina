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

class Illuminutils:
    def __init__(self, ini):
        with open(ini) as f:
            self.params = yaml.load(f)

    def save(self, data, dstname, scale_factor=1.):
        scaled_data = [ d*scale_factor for d in data ]
        ds = hdftools.from_domain(self.params,data)
        ds.save(dstname)

    def srtm(self, folder):
        files = glob(folder+"/*.hgt")
        print "    ".join(map(str,files))
        data = [ warp(files, self.params['srs'], extent) \
            for extent in self.params['extents'] ]
        self.save(data,"srtm")

    def modis(self, folder, band_n):
        files = glob(folder+"/*.hdf")
        fname = "refl_b%02d" % band_n
        band_names = map( lambda f: MYD09A1_band_name(f, band_n), files )
        print "    ".join(map(str,band_names))
        data = [ warp(band_names, self.params['srs'], extent) \
            for extent in self.params['extents'] ]
        self.save(data, fname, scale_factor=0.0001)

    def viirs(self, folder):
        files = glob(folder+"/*.tif")
        print "    ".join(map(str,files))
        data = [ warp(files, self.params['srs'], extent) \
            for extent in self.params['extents'] ]
        self.save(data,"stable_lights")

if __name__ == "__main__":
    I = Illuminutils(glob("*.ini")[0])

    I.srtm("SRTM")
    I.viirs("VIIRS-DNB")
    for i in range(1,8):
        I.modis("MODIS",i)
