#!/usr/bin/env python2
"""gdal2pgm

:Author: Jean-Denis Giguere

Please use this command to get information about usage
$ gdal2pgm -h

:Author: Jean-Denis Giguere

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""


from optparse import OptionParser
import sys
import os
import math

import pyproj
try:
    import osgeo.gdal as gdal
    from osgeo.gdalconst import *
except ImportError:
    import gdal
    from gdalconst import *

class Gdal2pgmException(Exception):
    def __init__(self, val=''):
        Exception.__init__(self, val)

class ArgumentError(SystemExit):
    """Exception to raise when bad or inclompete parameter are given
    to gdal2pgm
    """
    def __init__(self, value='DEFAULT'):
        self.code = {'DEFAULT': 100,
                     'NOARGS': 101,
                     'PIXSIZ': 102,
                     'BBOX': 103,
                     'SPATIAL': 106,
                  }
        self.message = {100: 'Generic Error',
                        101: 'No argument provided',
                        102: 'Pixel size is not valid',
                        103: 'Invalid bounding box',
                        106: 'Impossible to determine spatial domain',
                       }
        sys.stderr.write('Argument error: ' + self.message[self.code[value]] +
                         '\n')
        SystemExit.__init__(self,self.code[value])
    def __str__(self):
        return(self.message[self.code[value]])

class SpatialParams(object):
    def __init__(self, srs=None, lat0=None, lon0=None, pixsiz=None, 
                 width=None, height=None, bbox=None):
        self.define_srs(srs)
        self.lat0=lat0
        self.lon0=lon0
        self.pixsiz=pixsiz
        self.width=width
        self.height=height
        self.bbox=bbox

    def __str__(self):
        return "lat0:%.5f, lon0:%.5f, pixsiz:%f, width:%i, height:%i, \
                bbox:" % (self.lat0, self.lon0, self.pixsiz, self.width,
                          self.height) + self.bbox.__str__()

    def define_srs(self,srs):
        """Use proj string to define Proj object
        
        See pyproj documentation for more information
        """
        if srs is not None:
            self.srs = pyproj.Proj(init=srs)
        else:
            self.srs = None

    def compute(self):
        if self.lat0 is None:
            if self.srs is None:
                raise Gdal2pgmException('No srs avalaible')
            southwest_coords = self.bbox[0:2]
            self.lat0 = self.srs(southwest_coords[0], southwest_coords[1],
                                 inverse=True)[1]
        if self.lon0 is None:
            southwest_coords = self.bbox[0:2]
            self.lon0 = self.srs(southwest_coords[0], southwest_coords[1],
                                 inverse=True)[0]
        if self.bbox is None:
            if self.lat0 is None:
                raise Gdal2pgmException('lat0 is not available')
            if self.lon0 is None:
                raise Gdal2pgmException('lon0 is not available')
            if self.pixsiz is None:
                raise Gdal2pgmException('pixsiz is not available')
            if self.width is None:
                raise Gdal2pgmException('width is not available')
            if self.height is None:
                raise Gdal2pgmException('height is not available')
            southwest_coords = self.srs(self.lon0, self.lat0)
            width_in_meter = self.pixsiz * self.width
            height_in_meter = self.pixsiz * self.height
            self.bbox = (southwest_coords[0], southwest_coords[1],
                         southwest_coords[0] + width_in_meter,
                         southwest_coords[1] + height_in_meter)

        if self.width is None:
            # We suppose that srs units is meter!
            self.width = int(math.ceil( 
                float(self.bbox[2] - self.bbox[0])/self.pixsiz)
            )

        if self.height is None:
            # We suppose that srs units is meter!
            self.height = int(math.ceil( 
                float(self.bbox[3] - self.bbox[1])/self.pixsiz)
            )
        return self

    def corrected_bbox(self):
        """Return a tuple containing a bounding box calculed using
        south-west corner, pixel size, width and height"""

        minx = self.bbox[0]
        miny = self.bbox[1]
        maxx = minx + (self.width * self.pixsiz)
        maxy = miny + (self.height * self.pixsiz)

        corr_bbox = (minx, miny, maxx, maxy)

        return corr_bbox



def no_args_given(argv):
    """Return true if there is no argument"""
    
    if len(argv) == 0:
        return(True)
    else:
        return(False)

def check_pixsiz(option, opt_str, value, parser):
    """Check if pixel size is adequate
    
    This function is use as callback by OptionParser
    """
    is_too_small = (value <= 0)
    doesnt_exist = (value is None)
    if is_too_small or doesnt_exist:
        raise ArgumentError('PIXSIZ')
    else:
        setattr(parser.values, option.dest, value)

def parse_bbox(option, opt_str, value, parser):
    """Parse bounding box argument

    This function is use as callback by OptionParser
    """
    bbox_as_string = value
    bbox = tuple([float(coord) for coord in bbox_as_string.split(',')])
    if (bbox[0] >= bbox[2]) or (bbox[1] >= bbox[3]):
        raise ArgumentError('BBOX')
    else:
        setattr(parser.values, option.dest, bbox)


def parseOptions():
    parser = OptionParser()
    parser.add_option('-r', '--raster', type='string', dest='raster_name',
                     help='Name of raster containing data')
    parser.add_option('--lat0', type='float', dest='lat0', 
                      help='N-W pixel center latitude')
    parser.add_option('--lon0', type='float', dest='lon0',
                      help='N-W pixel center longitude')
    parser.add_option('-p', '--pixsiz', dest='pixsiz', action='callback',
                      callback=check_pixsiz, type='float', 
                      help='Pixel size (in meter)')
    parser.add_option('-b', '--bbox', dest='bbox', action='callback',
                      callback=parse_bbox, type='string',
                      help='Bounding box. Value are separated by commas')
    parser.add_option('--srs', type='string', dest='srs',
                      help='Spatial reference system string')
    (options, args) = parser.parse_args()
    if no_args_given(sys.argv[1:]): 
        print parser.get_usage()
        raise ArgumentError('NOARGS')
    if options.pixsiz is None:
        raise ArgumentError('PIXSIZ')
    else:
        pass

    return options

def spatialParamsFromOptions(options):
    using_bbox = ( options.bbox is not None and 
                  options.pixsiz is not None )
    using_spatial_params = (options.lat0 is not None and
                            options.lon0 is not None and
                            options.pixsiz is not None and
                            options.width is not None and
                            options.height is not None )
    if using_bbox:
        unparse_spatial_params = SpatialParams(
            bbox=options.bbox, pixsiz=options.pixsiz, srs=options.srs)
    elif using_spatial_params:
        unparse_spatial_params = SpatialParams(
            lat0=options.lat0, lon0=options.lon0, pixsiz=options.pixsiz,
            width=options.width, height=options.height, srs=options.srs)
    else:
        raise ArgumentError('SPATIAL')

    spatial_param = unparse_spatial_params.compute()
    return spatial_param

def checkRasterDomain(geotransform, xsize, ysize, spatial_params):
    """Check if raster cover entire spatial domain"""
    northwest_coords = (geotransform[0], geotransform[3])
    southeast_coords = (geotransform[0] + xsize*geotransform[1],
                        geotransform[3] + ysize*geotransform[5])

    north_is_ok = northwest_coords[1] > spatial_params.bbox[3]
    east_is_ok = southeast_coords[0] > spatial_params.bbox[2]
    south_is_ok = southeast_coords[1] < spatial_params.bbox[1]
    west_is_ok = northwest_coords[0] < spatial_params.bbox[0]
    domain_is_ok = north_is_ok and east_is_ok and south_is_ok and west_is_ok
    if not domain_is_ok:
        raise Gdal2pgmException("Raster don't cover spatial domain")


def checkRaster(raster_name, spatial_params):
    """Check if raster sastisfy all conditions to create new pgm"""
    if raster_name is None:
        raise Gdal2pgmException('Raster name must be provided')
    if not os.access(raster_name, os.R_OK):
        raise ArgumentError()
    dataset = gdal.Open(raster_name, GA_ReadOnly)
    if dataset is None:
        raise Gdal2pgmException('No dataset found in file')
    checkRasterDomain(dataset.GetGeoTransform(), dataset.RasterXSize,
                      dataset.RasterYSize, spatial_params)


def getValuesFromRaster(raster_name, spatial_params):
    """We must compute position of bounding box relative
    to image. We need:
        * xoff - x offset (in pixels relative to north-west corner)
        * yoff - y offset (in pixels relative to north-west corner)
        * win_xsize - number of pixels of original raster to copy (x)
        * win_ysize - number of pixels of original raster to copy (y)
        * buf_xsize - number of pixels of destination raster (x)
        * buf_ysize - number of pixels of destination raster (y)

    This function return a numarray with raster values.
    """
    pass



def createPgm(raster_name, spatial_params):
    raster_values = getValuesFromRaster(raster_name, spatial_params)
                                               

def main():
    options = parseOptions()
    spatial_params = spatialParamsFromOptions(options)
    checkRaster(options.raster_name, spatial_params)
    createPgm(options.raster_name, spatial_params)

if __name__ == '__main__':
    main()
