#!/usr/bin/env python3

from glob import glob
import yaml, h5py
import hdftools
from subprocess import call
import os
import numpy as np
from PIL import Image
import argparse

parser = argparse.ArgumentParser(description="""Warps sattelite imagery.

If no parameters are given, executes the standard pipeline.""")
parser.add_argument( "output_name", default=None, nargs='?', help="Name of the output file." )
parser.add_argument( "filename", nargs="*", help="Files to warp." )

p = parser.parse_args()

if p.output_name != None and len(p.filename) == 0:
    print("ERROR: If an output name is given, files to process must also be provided.")
    raise SystemExit

def OpenTIFF(path):
    return np.array(Image.open(path))

def warp(srcfiles, projection, extent):
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
        cmd += ['-spat_srs',projection]
        cmd += ['-t_srs',projection]
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

        if not os.path.isfile("tmp_merge.shp"):
            for fname in glob("tmp_select.*"):
                os.rename(fname,fname.replace("select","merge"))

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

def correction_filenames(srcfiles):
    return [ fname.replace(fname.split('_')[3],"zero_correction") \
        .replace("avg_rade9h.tif","csv") \
        for fname in srcfiles ]

def convert_correction_data(srcfiles):
    corr_files = np.unique(correction_filenames(srcfiles))

    data = np.nanmean([ np.loadtxt(fname,delimiter=',') for fname in corr_files ],0)
    data[np.isnan(data)] = -9999

    with open("VIIRS-DNB/correction.asc",'w') as f:
        f.write("NCOLS 72\n")
        f.write("NROWS 28\n")
        f.write("XLLCORNER -180\n")
        f.write("YLLCORNER -65\n")
        f.write("CELLSIZE 5\n")
        f.write("NODATA_VALUE -9999\n")
        np.savetxt(f, data)

    with open("VIIRS-DNB/correction.prj",'w') as f:
        f.write('GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]]')

with open(glob("*.ini")[0]) as f:
    params = yaml.safe_load(f)

if len(p.filename):
    data = [ warp(p.filename, params['srs'], extent) \
        for extent in params['extents'] ]
    save(params, data, p.output_name)

else:
    if os.path.isfile("GHSL.zip"):
        print("Found GHSL.zip file, processing.")
        data = [ warp(["/vsizip/GHSL.zip/GHSL.tif"], params['srs'], extent) \
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
    data = [ warp(files, params['srs'], extent) \
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

        correction = np.all([ os.path.isfile(fname) \
            for fname in correction_filenames(files) ])
        if not correction:
            print("WARNING: Could not find correction files that matched the VIIRS files.")
            print("If you indend to use it, please validate that you have the right ones.")
            print("Note that only the VCMCFG dataset from VIIRS can be corrected.")

        data = [ warp(files, params['srs'], extent) \
            for extent in params['extents'] ]

        if correction:
            convert_correction_data(files)
            corr = [ warp(["VIIRS-DNB/correction.asc"], params['srs'], extent) \
                for extent in params['extents'] ]
            save(params, corr, "VIIRS_background")
            for l in range(len(data)):
                data[l] -= corr[l]

        save(params, data, "stable_lights")

        prep_shp(
            "hydropolys.zip/hydropolys.shp",
            params['srs'],
            params['extents'][-1]
        )
        data = [ rasterize("tmp_merge.shp", params['srs'], extent) \
            for extent in params['extents'] ]
        save(params, data, "water_mask")

for fname in glob("tmp_*"):
    os.remove(fname)

print("Done.")
