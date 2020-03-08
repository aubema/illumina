#!/usr/bin/env python3

import click
from glob import glob
import yaml, h5py
import hdftools
from subprocess import call
import os
import numpy as np
from PIL import Image

def OpenTIFF(path):
    return np.array(Image.open(path))

def warp_files(srcfiles, projection, extent):
    tmpfile = "tmp_warp.tiff"
    if os.path.isfile(tmpfile):
        os.remove(tmpfile)

    cmd  = ['gdalwarp']
    cmd += ['-t_srs',projection]
    cmd += ['-te',
        extent["xmin"],
        extent["ymin"],
        extent["xmax"],
        extent["ymax"]
    ]
    cmd += ['-tr',
        extent["pixel_size"],
        extent["pixel_size"]
    ]
    cmd += ['-r','cubicspline']
    cmd += ['-dstnodata','0.']
    cmd += srcfiles
    cmd += [tmpfile]
    cmd = list(map(str,cmd))
    print("EXECUTING :", ' '.join(cmd))
    call(cmd)

    return OpenTIFF(tmpfile)

def prep_shp(infile, projection, extent):
        cmd  = ['ogr2ogr']
        cmd += ['-spat',
            extent['xmin'],
            extent['ymin'],
            extent['xmax'],
            extent['ymax']
        ]
        cmd += ['-spat_srs','+init='+projection]
        cmd += ['-t_srs','+init='+projection]
        cmd += ['tmp_select.shp']
        cmd += ['/vsizip/'+os.path.abspath(infile)]
        cmd = list(map(str,cmd))
        print("EXECUTING :", ' '.join(cmd))
        call(cmd)

        cmd  = ['ogr2ogr']
        cmd += ['tmp_merge.shp']
        cmd += ['tmp_select.shp']
        cmd += ['-dialect','sqlite']
        cmd += ['-sql','SELECT ST_Union(geometry) AS geometry FROM tmp_select']
        print("EXECUTING :", ' '.join(cmd))
        call(cmd)

def rasterize(shpfile, projection, extent):
    tmpfile = "tmp_rasterize.tiff"
    if os.path.isfile(tmpfile):
        os.remove(tmpfile)

    cmd  = ['gdal_rasterize']
    cmd += ['-i','-at'] # inverted, all touched
    cmd += ['-burn','1']
    cmd += ['-ot','Byte']
    cmd += ['-a_srs',projection]
    cmd += ['-te',
        extent["xmin"],
        extent["ymin"],
        extent["xmax"],
        extent["ymax"]
    ]
    cmd += ['-tr',
        extent["pixel_size"],
        extent["pixel_size"]
    ]
    cmd += [shpfile]
    cmd += [tmpfile]
    cmd = list(map(str,cmd))
    print("EXECUTING :", ' '.join(cmd))
    call(cmd)

    return OpenTIFF(tmpfile)

def save(params, data, dstname, scale_factor=1.):
    scaled_data = [ d*scale_factor for d in data ]
    ds = hdftools.from_domain(params,scaled_data)
    ds.save(dstname)

@click.command()
def warp():
    """Warps the satellite imagery.

    Warps the satellite imagery based on the domain defined in
    'domain.ini'.

    \b
    Requires the folowing data:
        Unzipped SRTM data in a folder named 'SRTM'.
        If used, VIIRS data in a volder named 'VIIRS-DNB'.
        If VIIRS data is used, the 'hydropolys.zip' file.
    """
    with open(glob("*.ini")[0]) as f:
        params = yaml.safe_load(f)

    if os.path.isfile("GHSL.zip"):
        print("Found GHSL.zip file, processing.")
        data = [ warp_files(["/vsizip/GHSL.zip/GHSL.tif"], params['srs'], extent) \
            for extent in params['extents'] ]
        save(params, data, "obstf")
    else:
        print("WARNING: Could not find GHSL.zip file.")
        print("If you don't intend to use it, you can safely ignore this.")

    files = sorted(glob("SRTM/*.hgt"))
    if not len(files):
        print("ERROR: Could not find SRTM file(s), aborting.")
        raise SystemExit
    print("    ".join(map(str,files)))
    data = [ warp_files(files, params['srs'], extent) \
        for extent in params['extents'] ]
    save(params, data, "srtm")

    files = sorted(glob("VIIRS-DNB/*.tif"))
    if not len(files):
        print("WARNING: Did not find VIIRS file(s).")
        print("If you don't intend to use zones inventory, you can safely ignore this.")
    else:
        if not os.path.isfile("hydropolys.zip"):
            print("ERROR: Could not find hydropolys.zip file, aborting.")
            raise SystemExit

        print("    ".join(map(str,files)))
        data = [ warp_files(files, params['srs'], extent) \
            for extent in params['extents'] ]
        save(params, data, "stable_lights")

        prep_shp(
            "hydropolys.zip/hydropolys.shp",
            params['srs'],
            params['extents'][-1]
        )
        data = [ rasterize("tmp_merge.shp", params['srs'], extent) \
            for extent in params['extents'] ]
        save(params, data, "water_mask")

        for fname in glob("tmp*"):
            os.remove(fname)

        print("Done.")
