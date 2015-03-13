# -*- coding: utf-8 -*-
""""Transform raster to illumina pgm format

:Author: Jean-Denis Giguère
:Licence: GPL v3 or later
:Date: 2013-02-16
"""

import logging
import subprocess
import sys
from osgeo import gdal
from osgeo.gdalconst import *
import numpy

class Pgm(object):
    "A raster to transform to special Illumina PGM format"
    def __init__(self):
        self._srcfile = None
        self._dstfile = None
        self.extra_scale_factor = 1
        self.extra_offset = 0
        #TODO: Build spatial metadata from ini
        self._spatial_metadata = SpatialMetadata()
        #TODO: Check values below
        self._spatial_metadata.lat0 = 39.0273
        self._spatial_metadata.lon0 = -6.00354
        self._spatial_metadata.pixsize = 1000
        self._spatial_metadata.srs_string = '+init=epsg:3488'
        self._spatial_metadata.width = 450
        self._spatial_metadata.height = 300
        
    
        
    @property
    def srcfile(self):
        "Raster to reproject"
        return self._srcfile
        
    @srcfile.setter
    def srcfile(self,  fname):
        self._srcfile = fname
        
    @property
    def dstfile(self):
        "Name of reprojected raster"
        return self._dstfile
        
    @dstfile.setter
    def dstfile(self,  fname):
        self._dstfile = fname
        
    def write(self):
        self._gdal_array = array_from_gdal(self.srcfile)
        self._pgm_array = PgmReadyArray(self._gdal_array,  self.extra_scale_factor,  self.extra_offset)
        logging.info("Raster offset is %f \n" % (self._pgm_array.offset,))
        logging.info("Raster scale factor is %f \n" % (self._pgm_array.scale_factor,))
        logging.debug("Raster pgm representation...")
        logging.debug(self._pgm_array.pgm_array)
        logging.debug("Saving the pgm")
        pgm_from_pgm_ready_array(self._pgm_array, self._spatial_metadata, self.dstfile)



class PgmReadyArray:
    def __init__(self, array, scale_factor, offset):
        array = array*scale_factor + offset
        self.maxint = 65535
        self.min = array.min()
        logging.debug('Le min est %f' % (self.min, ))
        self.max = array.max()
        logging.debug('Le max est %f' % (self.max, ))
        self.offset = self.min
        self.range = self.max - self.min
        self.scale_factor = self.range / float(self.maxint) 
        array_float = numpy.array(array, dtype=float)
        array_wo_offset = array_float - float(self.offset)
        array_scaled_wo_offset = array_wo_offset/self.scale_factor
        self.pgm_array = numpy.array(array_scaled_wo_offset, dtype=numpy.uint16)

class SpatialMetadata:
    'Stub class replacing SpatialParams without calculation'
    def __init___(self):
        self.lat0 = None
        self.lon0 = None
        self.pixsize = None
        self.srs_string = None
        self.width = None
        self.height = None


def array_from_gdal(filename):
    'Return array from gdal readable file'
    dataset = gdal.Open(filename, GA_ReadOnly)
    band = dataset.GetRasterBand(1)
    input_array = band.ReadAsArray()

    return input_array

def pgm_from_pgm_ready_array(pgm_ready_array, spatial_metadata, filename):
    file_obj = open(filename, 'w')
    file_obj.write('P2\n')
    file_obj.write("# lat0 %f\n" % (spatial_metadata.lat0,) )
    file_obj.write("# lon0 %f\n" % (spatial_metadata.lon0,) )
    file_obj.write("# pixsize %f\n" % (spatial_metadata.pixsize,) )
    file_obj.write("# srs %s\n" % (spatial_metadata.srs_string,) )
    file_obj.write("# gain %f\n" % (pgm_ready_array.scale_factor,) )
    file_obj.write("# offset %f\n" % (pgm_ready_array.offset,) )
    file_obj.write("%i %i 65535\n" % (spatial_metadata.width, 
                spatial_metadata.height) )

    #Waring! Numpy array indexing differ from fortran indexing!
    #You should consider each other as transposed array
    for height_iter in range(pgm_ready_array.pgm_array.shape[0]):
        for width_iter in range(pgm_ready_array.pgm_array.shape[1]):
            current_pixel = pgm_ready_array.pgm_array[height_iter, width_iter]
            file_obj.write("%i\n" % (current_pixel,) )

    file_obj.close()

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)

 

