#!/usr/bin/env python2

import subprocess

from srtm_domain import tiles

if __name__ == '__main__':
    srtm_filenames = [tile + '.hgt' for tile in tiles] # List comprehension
    nodata = -32678 # From srtm documentation
    mosaic_name = 'srtm1.tif'
    merge_command = ['gdal_merge.py', '-o', mosaic_name,  '-of', 'GTiff', 
            '-n', str(nodata)]
    merge_command = merge_command + srtm_filenames
    subprocess.call(merge_command)

    mosaic_with_nodata_name = 'srtm.tif'
    set_nodata_command = ['gdalwarp', '-srcnodata', str(nodata),
            '-dstnodata', str(nodata), mosaic_name, mosaic_with_nodata_name]
    subprocess.call(set_nodata_command)
