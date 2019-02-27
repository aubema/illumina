#!/usr/bin/env python
# -*- coding: utf-8 -*-
""""Coordinate lookup utility fonctions

Input:
    Illumina Parameters ini file (hardcoded)
    Projection definition (hardcoded)

Output:
    Lookup table between row and column number and coordinates

:Author: Jean-Denis Giguère
:Licence: GPL v3 or later
:Date: 2013-02-10
"""


import math
import pyproj
from illuminutils.domain import cell_height, cell_center_height
from illuminutils.experiments.illuminini import IlluminaParameters

params = IlluminaParameters('palomar.ini')

#TODO: Read from parameter
projection = pyproj.Proj('+init=epsg:3488')

def get_lat_lon_height_from_rcz(row,column,zlevel):
    pass

def get_row_column_zlevel_from_llz(latitude,  longitude,  height):
    pass

def get_row_column_zlevel_from_xyz(x, y,  height,  style='fortran'):
    #TODO: Make robust
    minx = params.bbox[0]
    miny = params.bbox[1]
    pixsize = params.pixsize
    column = int(math.floor((x - minx)/pixsize))
    row = int(math.floor((y - miny)/pixsize))
    if style == 'fortran':
        column = column + 1
        row = row + 1
    zlevel = None
    return (row, column, zlevel)

if __name__ == '__main__':
    lookup_fname = 'row_column_to_x_y_lon_lat.csv'
    lookup_file = open(lookup_fname, 'w')
    minx = params.bbox[0]
    maxy = params.bbox[3]
    pixsize = params.pixsize
    #TODO: From read column et row number in ini
    lookup_file.write('"col","row","x","y","lon","lat"\n')
    for col in range(400):
        for row in range(300):
            x = minx + (col + 0.5) * pixsize
            y = maxy - (row + 0.5) * pixsize
            (lon, lat) = projection(x, y, inverse=True)
            line = "%i,%i,%f,%f,%f,%f\n" % (col, row, x, y, lon, lat)
            lookup_file.write(line)
    lookup_file.close()

