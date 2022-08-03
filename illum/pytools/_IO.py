# Functions related to IO opperations
# Author: Alexandre Simoneau


import errno
import os

import astropy.io.fits as fits
import numpy as np
from osgeo import gdal, gdal_array

# Georeferenced FITS


def load_fits(filename):
    """Loads a FITS file.

    Returns a 2-ple containing
      - A list of arrays defining each axis in order (x,y,z,...)
      - The contained data in transposed order (...,z,y,x)
    """
    hdu = fits.open(filename)[0]

    ax = [
        np.linspace(
            hdu.header["CRVAL%d" % (i + 1)],
            hdu.header["CRVAL%d" % (i + 1)]
            + hdu.header["CDELT%d" % (i + 1)]
            * (hdu.header["NAXIS%d" % (i + 1)] - 1),
            hdu.header["NAXIS%d" % (i + 1)],
        )
        for i in range(hdu.header["NAXIS"])
    ]

    return ax, hdu.data.T[:, ::-1].T


def save_fits(axis, data, filename):
    """Save an array to a fits file. Must be at least 2D.

    axis : a list of 2-tuple containing the base value and the increment for
           each axis.
    data : data array. The dimensions must be ordered (...,z,y,x)
    filename : name of the file to create
    """
    hdu = fits.PrimaryHDU()
    hdu.data = data.T[:, ::-1].T
    for i in range(len(axis)):
        hdu.header["CRPIX%d" % (i + 1)] = 1
        hdu.header["CRVAL%d" % (i + 1)] = axis[i][0]
        hdu.header["CDELT%d" % (i + 1)] = axis[i][1]
    hdu.writeto(filename, clobber=True)


# GeoTiff


def load_geotiff(filename):
    """Open a georeferenced tiff image as a numpy array.

    Returns the data array, the projection and the geotransform."""

    if not os.path.isfile(filename):
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT), filename
        )

    src = gdal.Open(filename, gdal.GA_Update)
    arr = src.GetRasterBand(1).ReadAsArray()

    return arr, src.GetProjection(), src.GetGeoTransform()


def save_geotiff(filename, arr, projection, geotransform, compression="NONE"):
    """Saves a 2D numpy data array as a georeferenced tiff image.

    Needs the projection (WKT string) and the geotransform coefficients.

    Compatible data types are:
        '(u)int[8,16,32]', 'float[32,64]', 'complex[64,128]'

    Different compressions algorithms can be used:
        'LZW', 'DEFLATE', 'PACKBITS', 'JPEG' or 'NONE'.
    'JPEG' compression shoudn't be used when accurate values are needed."""

    options = ["LZW", "PACKBITS", "DEFLATE", "JPEG", "NONE"]
    if compression not in options:
        raise ValueError(
            "Compression needs to be one of "
            + ", ".join("'" + opt + "'" for opt in options[:-1])
            + f" or '{options[-1]}'."
        )
    if compression == "JPEG":
        print(f"Warning: Lossy compression used when saving {filename}.tiff")

    shape = arr.shape[::-1] + (1,)
    driver = gdal.GetDriverByName("GTiff")

    type_code = gdal_array.NumericTypeCodeToGDALTypeCode(arr.dtype)
    if type_code is None:
        raise TypeError("Unsupported data type.")

    dst_dataset = driver.Create(
        filename + ".tiff", *shape, type_code, ["COMPRESS=" + compression]
    )
    dst_dataset.SetGeoTransform(geotransform)
    dst_dataset.SetProjection(projection)
    dst_dataset.GetRasterBand(1).WriteArray(arr)
    dst_dataset.FlushCache()


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
