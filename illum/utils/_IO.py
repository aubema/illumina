# Functions related to IO opperations
# Author: Alexandre Simoneau


import errno
import os

import astropy.io.fits
import astropy.wcs
import numpy as np
from osgeo import gdal, gdal_array

# Georeferenced FITS


def load_fits(filename, c=True):
    """Loads a FITS file."""
    hdu = astropy.io.fits.open(filename)[0]
    try:
        wcs = astropy.wcs.find_all_wcs(hdu.header, fix=False)[0]
    except IndexError:
        return hdu.data, None
    c = wcs.all_pix2world([[-0.5, -0.5], [0.5, 0.5]], 0)
    gt = (c[0, 0], c[1, 0] - c[0, 0], 0, c[0, 1], 0, c[1, 1] - c[0, 1])
    return hdu.data, gt


# GeoTiff


def load_geotiff(filename):
    """Open a georeferenced tiff image as a numpy array.

    Returns the data array, the projection and the geotransform."""
    rst = rasterio.open(filename)
    return rst.read(1), rst.crs, rst.get_transform()


def save_geotiff(filename, arr, **kwargs):
    """Saves a 2D numpy data array as a georeferenced tiff image."""
    if "dtype" not in kwargs:
        kwargs["dtype"] = str(arr.dtype)
    kwargs["height"], kwargs["width"] = arr.shape
    profile = rasterio.profiles.DefaultGTiffProfile(count=1, **kwargs)
    with rasterio.open(filename, "w", **profile) as f:
        f.write(arr.astype(kwargs["dtype"]), 1)


# Fortran binary files


def load_bin(filename, dtype=np.float32):
    """Load a ILLUMINA binary file.

    Returns the data as an array."""
    with open(filename) as f:
        shape = np.fromfile(f, dtype=np.uint32, count=4)[1:-1][::-1]
        data = np.fromfile(f, dtype=np.float32, count=-1)[1::3]
    return data.reshape(shape).astype(dtype)


def save_bin(filename, data):
    """Saves a numpy data array as an ILLUMINA binary file."""
    data = data.astype(np.float32)
    shape = data.shape[::-1]
    size = data.size
    data_flat = data.flatten()
    filler = np.ones(size, dtype=np.float32) * 5.6e-45

    head = np.array((8,) + shape + (8,)).astype(np.uint32)
    body = np.array([filler, data_flat, filler]).T.flatten()

    with open(filename, "w") as f:
        head.tofile(f)
    with open(filename, "a") as f:
        body.tofile(f)
