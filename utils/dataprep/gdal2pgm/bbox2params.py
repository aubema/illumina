#!/usr/bin/env python
"""bbox2params.py - Compute spatial parameters required by illumina

When using Illumina, you need to know the following information:
    * Coordinates of the south-west corner of the spatial domain
    * Pixels size (pixels must be square)
    * Width of spatial domain, in pixels
    * Height of spatial domain, in pixels

Domain definition is given by:
    * Bounding box of spatial domain
    - or -
    * Domain center in latitude, longitude
    * Domain size

bbox2params compute required paramaters using:
    * Spatial reference system of the given bounding box
    * Pixels size

For help on usage, please use
bbox2params.py -h

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

import sys
import gdal2pgm
import pyproj

from optparse import OptionParser, OptionGroup

def parseOptions():
    """Parse command line options

    Please refer to the standard python module optparse
    """
    parser = OptionParser()
    parser.add_option('-p', '--pixsiz', dest='pixsiz', action='callback',
                      callback=gdal2pgm.check_pixsiz, type='float',
                      help='Pixel size (in meter)')
    parser.add_option('-b', '--bbox', dest='bbox', action='callback',
                      callback=gdal2pgm.parse_bbox, type='string',
                      help='Bounding box. Value are separated by commas')
    parser.add_option('--srs', type='string', dest='srs',
                      help='Spatial reference system string')
    parser.add_option('-l', dest='print_ll', action='store_true',
                       help='Print bbox coordinates in lat/lon')
    domain_from_center_grp = OptionGroup(parser, "Domain definiton by center",
                    "Instead of giving the bounding box, you can compute it "
                    "from a center point.")
    domain_from_center_grp.add_option('--center-latitude',
                     dest='center_latitude', type='float',
    	    	     help='Latitude of the center of the domain')
    domain_from_center_grp.add_option('--center-longitude',
                     dest='center_longitude', type='float',
                     help='Latitude of the center of the domain')
    domain_from_center_grp.add_option('--width', type='float',
                     help='Domain width [in meters]')
    domain_from_center_grp.add_option('--height', type='float',
                     help='Domain height [in meters]')
    parser.add_option_group(domain_from_center_grp)

    (options, args) = parser.parse_args()

    return options

def addExtraOptions(options):
    """gdal2pgm use extra options, we must add some of them
    for compatibility
    """
    options.lat0 = None
    options.lon0 = None

    return options

def bbox_corners_string(bbox, srs):
    s_w_corner = srs(bbox[0],bbox[1], inverse=True)
    s_e_corner = srs(bbox[2],bbox[1], inverse=True)
    n_w_corner = srs(bbox[0],bbox[3], inverse=True)
    n_e_corner = srs(bbox[2],bbox[3], inverse=True)
    bbox_string = 'Bounding box corners coordinates: \n' + \
            "South-west corner: (%.5f, %.5f) \n" % s_w_corner  + \
            "South-east corner: (%.5f, %.5f) \n" % s_e_corner + \
            "North-west corner: (%.5f, %.5f) \n" % n_w_corner + \
            "North-east corner: (%.5f, %.5f) \n" % n_e_corner

    return bbox_string


def center_ll_2_bbox(options):
    """Compute the bbox parameter given a center latitude, longitude
    and a spatial reference system"""
    if (options.center_latitude is None or options.center_longitude is None):
        print "Only one parameter of center latitude and center longitude is given."
        print "Please provide both or use the bbox paramater"
        sys.exit(1)
    if (options.width is None or options.height is None):
        print "Width and height are required to compute bounding box using", \
                "center coordinates"
    srs_proj = pyproj.Proj(init=options.srs)
    srs_center = srs_proj(options.center_longitude, options.center_latitude)
    bbox = [srs_center[0] - float(options.width)/2,
            srs_center[1] - float(options.height)/2,
            srs_center[0] + float(options.width)/2,
            srs_center[1] + float(options.height)/2]
    bbox_str = "%f,%f,%f,%f " % (bbox[0], bbox[1], bbox[2], bbox[3])

    return bbox



def main():
    options = parseOptions()
    options = addExtraOptions(options)

    if (options.srs is None):
        print 'Spatial reference system is required'
        sys.exit(2)

    if (options.center_latitude is not None or
            options.center_longitude is not None):
    	    print 'Latitude is given'
    	    options.bbox = center_ll_2_bbox(options)

    spatial_params = gdal2pgm.spatialParamsFromOptions(options)
    print """South-west corner latitude (lat0): %.5f
    South-west corner (lon0): %.5f
    Pixel size (pixsiz): %f
    Spatial domain width: %i
    Spatial domain height: %i
    """ % (spatial_params.lat0, spatial_params.lon0,
           spatial_params.pixsiz, spatial_params.width,
           spatial_params.height)
    print 'Corrected bounding box is ' + \
            spatial_params.corrected_bbox().__str__() + '\n'
    if options.print_ll:
        print bbox_corners_string(spatial_params.bbox, spatial_params.srs)

if __name__ == '__main__':
    main()

