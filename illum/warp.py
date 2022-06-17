#!/usr/bin/env python3

import os
from glob import glob
from subprocess import call

import click
import numpy as np
import yaml
from illum import MultiScaleData as MSD
from PIL import Image


def OpenTIFF(path):
    return np.array(Image.open(path))


def warp_files(srcfiles, projection, extent):
    tmpfile = "tmp_warp.tiff"
    if os.path.isfile(tmpfile):
        os.remove(tmpfile)

    cmd = [
        "gdalwarp",
        "-t_srs",
        projection,
        "-te",
        extent["xmin"],
        extent["ymin"],
        extent["xmax"],
        extent["ymax"],
        "-tr",
        extent["pixel_size"],
        extent["pixel_size"],
        "-r",
        "cubicspline",
        "-dstnodata",
        "0.",
        *srcfiles,
        tmpfile,
    ]
    cmd = list(map(str, cmd))
    print("EXECUTING :", " ".join(cmd))
    call(cmd)

    return OpenTIFF(tmpfile)


def prep_shp(infile, projection, extent):
    cmd = [
        "ogr2ogr",
        "-spat",
        extent["xmin"],
        extent["ymin"],
        extent["xmax"],
        extent["ymax"],
        "-spat_srs",
        projection,
        "-t_srs",
        projection,
        "tmp_select.shp",
        "/vsizip/" + os.path.abspath(infile),
    ]
    cmd = list(map(str, cmd))
    print("EXECUTING :", " ".join(cmd))
    call(cmd)

    cmd = [
        "ogr2ogr",
        "tmp_merge.shp",
        "tmp_select.shp",
        "-dialect",
        "sqlite",
        "-sql",
        "SELECT ST_Union(geometry) AS geometry FROM tmp_select",
    ]
    print("EXECUTING :", " ".join(cmd))
    call(cmd)

    if not os.path.isfile("tmp_merge.shp"):
        for fname in glob("tmp_select.*"):
            os.rename(fname, fname.replace("select", "merge"))


def rasterize(shpfile, projection, extent):
    tmpfile = "tmp_rasterize.tiff"
    if os.path.isfile(tmpfile):
        os.remove(tmpfile)

    cmd = [
        "gdal_rasterize",
        "-i",  # inverted
        "-at",  # all touched
        "-burn",
        "1",
        "-ot",
        "Byte",
        "-a_srs",
        projection,
        "-te",
        extent["xmin"],
        extent["ymin"],
        extent["xmax"],
        extent["ymax"],
        "-tr",
        extent["pixel_size"],
        extent["pixel_size"],
        shpfile,
        tmpfile,
    ]
    cmd = list(map(str, cmd))
    print("EXECUTING :", " ".join(cmd))
    call(cmd)

    return OpenTIFF(tmpfile)


def save(params, data, dstname, scale_factor=1.0):
    scaled_data = [d * scale_factor for d in data]
    ds = MSD.from_domain(params, scaled_data)
    ds.save(dstname)


def correction_filenames(srcfiles):
    try:
        return [
            fname.replace(fname.split("_")[3], "zero_correction").replace(
                "avg_rade9h.tif", "csv"
            )
            for fname in srcfiles
        ]
    except IndexError:
        return [fname.replace("avg_rade9h.tif", "csv") for fname in srcfiles]


def convert_correction_data(srcfiles):
    corr_files = np.unique(correction_filenames(srcfiles))

    data = np.nanmean(
        [np.loadtxt(fname, delimiter=",") for fname in corr_files], 0
    )
    data[np.isnan(data)] = -9999

    with open("VIIRS-DNB/correction.asc", "w") as f:
        f.write(
            "NCOLS 72\n"
            "NROWS 28\n"
            "XLLCORNER -180\n"
            "YLLCORNER -65\n"
            "CELLSIZE 5\n"
            "NODATA_VALUE -9999\n"
        )
        np.savetxt(f, data)

    with open("VIIRS-DNB/correction.prj", "w") as f:
        f.write(
            "GEOGCS["
            '"GCS_WGS_1984",'
            'DATUM["D_WGS_1984",'
            'SPHEROID["WGS_1984",6378137,298.257223563]],'
            'PRIMEM["Greenwich",0],'
            'UNIT["Degree",0.017453292519943295]'
            "]"
        )


@click.command(name="warp")
@click.argument("output_name", required=False)
@click.argument("infiles", required=False, nargs=-1)
def CLI_warp(output_name=None, infiles=[]):
    """Warps the satellite imagery.

    Warps the satellite imagery based on the domain defined in
    'domain.ini'.

    \b
    Requires the folowing data:
        Unzipped SRTM data in a folder named 'SRTM'.
        If used, VIIRS data in a volder named 'VIIRS-DNB'.
        If VIIRS data is used, the 'hydropolys.zip' file.

    Can also be used on specific files, in wich case an output name and a list
    of files to process must be given (the use of bash wildcards is encouraged)
    """
    warp(output_name, infiles)


def warp(output_name=None, infiles=[]):
    if output_name is not None and len(infiles) == 0:
        print(
            "ERROR: If an output name is given, "
            "files to process must also be provided."
        )
        raise SystemExit

    with open(glob("*.ini")[0]) as f:
        params = yaml.safe_load(f)

    if len(infiles):
        data = [
            warp_files(infiles, params["srs"], extent)
            for extent in params["extents"]
        ]
        save(params, data, output_name)

    else:
        if os.path.isfile("GHSL.zip"):
            print("Found GHSL.zip file, processing.")
            data = [
                warp_files(
                    ["/vsizip/GHSL.zip/GHSL.tif"], params["srs"], extent
                )
                for extent in params["extents"]
            ]
            save(params, data, "obstf")
        else:
            print("WARNING: Could not find GHSL.zip file.")
            print("If you don't intend to use it, you can safely ignore this.")

        files = sorted(glob("SRTM/*.hgt"))
        if not len(files):
            print("ERROR: Could not find SRTM file(s), aborting.")
            raise SystemExit
        print("    ".join(map(str, files)))
        data = [
            warp_files(files, params["srs"], extent)
            for extent in params["extents"]
        ]
        save(params, data, "srtm")

        files = sorted(glob("VIIRS-DNB/*.tif"))
        if not len(files):
            print("WARNING: Did not find VIIRS file(s).")
            print(
                "If you don't intend to use zones inventory, "
                "you can safely ignore this."
            )
        else:
            if not os.path.isfile("hydropolys.zip"):
                print("ERROR: Could not find hydropolys.zip file, aborting.")
                raise SystemExit

            print("    ".join(map(str, files)))

            correction = np.all(
                [
                    os.path.isfile(fname)
                    for fname in correction_filenames(files)
                ]
            )
            if not correction:
                print(
                    "WARNING: Could not find correction files "
                    "that matched the VIIRS files."
                )
                print(
                    "If you indend to use it, please validate "
                    "that you have the right ones."
                )
                print(
                    "Note that only the VCMCFG dataset from VIIRS "
                    "can be corrected."
                )

            data = [
                warp_files(files, params["srs"], extent)
                for extent in params["extents"]
            ]

            if correction:
                convert_correction_data(files)
                corr = [
                    warp_files(
                        ["VIIRS-DNB/correction.asc"], params["srs"], extent
                    )
                    for extent in params["extents"]
                ]
                save(params, corr, "VIIRS_background")
                for i, c in enumerate(corr):
                    data[i] -= c

            save(params, data, "stable_lights")

            prep_shp(
                "hydropolys.zip/hydropolys.shp",
                params["srs"],
                params["extents"][-1],
            )
            data = [
                rasterize("tmp_merge.shp", params["srs"], extent)
                for extent in params["extents"]
            ]
            save(params, data, "water_mask")

        for fname in glob("tmp_*"):
            os.remove(fname)

    print("Done.")
