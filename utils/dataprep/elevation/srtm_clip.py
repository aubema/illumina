#!/usr/bin/env python

import subprocess

output_srs = 'epsg:23030'
input_filename = 'srtm.tif'
reprojected_filename = 'srtm_utm30.tif'
corrected_bbox = [239981.58663178002, 4324182.311412051, 639981.58663178, 4624182.311412051] # from bbox2params
pixsiz = 50
clipped_filename = 'srtm_utm19_clipped.tif'

if __name__ == '__main__':
    reproject_cmd = ['gdalwarp', '-t_srs', output_srs, input_filename,
                reprojected_filename]
    subprocess.call(reproject_cmd)

    clip_cmd = ['gdalwarp', '-te'] + \
            [str(coord) for coord in corrected_bbox] + \
            ['-tr', str(pixsiz), str(pixsiz)] + \
                [reprojected_filename, clipped_filename]
    subprocess.call(clip_cmd)
