#!/usr/bin/env python

import os, subprocess as sub
from glob import glob
import gdal, pyproj
from pytools import save_pgm as _save_pgm

def params_parser(fname):
    with open(fname) as f:
        lines = filter(lambda line: ':' in line, f.readlines())

    params = map(lambda s: map(str.strip, s.split(':',1)), lines)
    return { p[0]:p[1] for p in params }

def warp(srcfiles, dstfile, srs=None, extent=None, resolution=None, resampling_method='cubicspline'):
    """gdalwarp interface

    Parameters equivalence:
        srcfiles: srcfile
        dstfile: dstfile
        srs: t_srs
        extent: te
        resolution: tr
        resampling_method: r

    List parameters (srcfiles, extent) must be given 
    as a list or a space separated string.
    """

    if srcfiles:
        if type(srcfiles) in [ list, tuple ]:
            srcfiles = ' '.join(srcfiles)

    if extent:
        if type(extent) in [ list, tuple ]:
            extent = ' '.join(map(str,extent))

    if resolution:
        resolution = ("%s " % resolution) *2

    cmd = ["gdalwarp"]

    if srs:
        cmd += ['-t_srs "%s"' % srs]
    if extent:
        cmd += ["-te " + extent]
    if resolution:
        cmd += ["-tr " + resolution]
    cmd += ["-r " + resampling_method]
    cmd += [srcfiles, dstfile]

    cmd_string = ' '.join(cmd)

    print "Executing command :"
    print cmd_string

    return sub.call(cmd_string,shell=True)

def get_MYD09A1_band_name(fname, band_n):
    out = sub.Popen(["gdalinfo",fname], stdin=sub.PIPE, stdout=sub.PIPE).communicate()[0]
    line = filter(lambda s: "SUBDATASET_%d_NAME" % band_n in s, out.split('\n'))[0]
    return line.split('=',1)[1]

class Illuminutils:
    cell_height = [ 0.5,0.6,0.72,0.86,1.04,1.26,1.52,1.84,2.22,2.68,3.24,3.92,
                    4.74,5.72,6.9,8.34,10.08,12.18,14.72,17.78,21.48,25.94,
                    31.34,37.86,45.74,55.26,66.76,80.64,97.42,117.68,142.16,
                    171.72,207.44,250.58,302.7,365.66,441.72,533.6,644.58,
                    778.66,940.62,1136.26,1372.6,1658.1,2002.98,2419.6,2922.88,
                    3530.84,4265.26,5152.44]

    def __init__(self, ini):
        P = params_parser(ini)

        P['pixsize'] = float(P['pixsize'])

        P['xmin'],P['ymin'],P['xmax'],P['ymax'] = map(float,P['bbox'].split())

        P['xsize'] = int((P['xmax']-P['xmin'])/P['pixsize'])
        P['ysize'] = int((P['ymax']-P['ymin'])/P['pixsize'])

        self.params = P

    def pgm_convert(self, srcfile, dstfile, maxint=65535, scale_factor=1.):
        gdal_file =gdal.Open(srcfile)
        band = gdal_file.GetRasterBand(1)
        data = band.ReadAsArray() * scale_factor
        
        pixsize = self.params['pixsize']
        xmin = self.params['xmin'] + 0.5 * pixsize 
        ymin = self.params['ymin'] + 0.5 * pixsize 
        lon0, lat0 = pyproj.Proj("+init="+self.params['srs'])(xmin, ymin, inverse=True)

        head = { 'x0':xmin, 'y0':ymin, 'lon0':lon0, 'lat0':lat0,
                 'pixsize':pixsize, 'srs':self.params['srs'] }
        p = data.shape[::-1] + (maxint,)

        _save_pgm(dstfile, head, p, data)

    def srtm(self, folder, cleanup=True):
        files = glob(folder+"/*.hgt")
        warp(files, "srtm.tif", self.params['srs'], self.params['bbox'], self.params['pixsize'])
        self.pgm_convert("srtm.tif","srtm.pgm")
        if cleanup:
            os.remove("srtm.tif")

    def modis(self, folder, band_n, cleanup=True):
        files = glob(folder+"/*.hdf")
        fname = "refl_b%02d" % band_n
        band_names = map(lambda f: get_MYD09A1_band_name(f, band_n), files)
        warp(band_names, fname+".tif", self.params['srs'], self.params['bbox'], self.params['pixsize'])
        self.pgm_convert(fname+".tif", fname+".pgm", scale_factor=0.0001)
        if cleanup:
            os.remove(fname+".tif")
        

    def viirs(self, folder, cleanup=True):
        files = glob(folder+"/*.tif")
        warp(files, "stable_lights.tif", self.params['srs'], self.params['bbox'], self.params['pixsize'])
        self.pgm_convert("stable_lights.tif","stable_lights.pgm")
        if cleanup:
            os.remove("stable_lights.tif")

    def lookup(self, fname='row_column_to_x_y_lon_lat.csv'):
        with open(fname, 'w') as lookup_file:
            pixsize = self.params['pixsize']
            lookup_file.write('"col","row","x","y","lon","lat"\n')
            for col in xrange(self.params['xsize']):
                for row in xrange(self.params['ysize']):
                    x = self.params['xmin'] + (col + 0.5) * pixsize
                    y = self.params['ymax'] - (row + 0.5) * pixsize
                    p1 = pyproj.Proj(init=self.params['srs'])
                    p2 = pyproj.Proj(init="epsg:4326") # WGS84
                    lon, lat = pyproj.transform(p1, p2, x, y)
                    line = "%i,%i,%f,%f,%f,%f\n" % (col, row, x, y, lon, lat)
                    lookup_file.write(line)


if __name__ == "__main__":
    I = Illuminutils(glob("*.ini")[0])

    I.lookup()
    I.srtm("SRTM")
    I.viirs("VIIRS-DNB")
    for i in [1,2,3,4]:
        I.modis("MODIS",i)


