# Mathematical functions
# Author: Alexandre Simoneau


import numpy as np


def deg2mx(deg, lat):
    """Converts a longitudinal measurement from degrees to meters."""
    a = 6378137.0
    b = 6356752.3
    r = a * ((b / a * np.tan(np.deg2rad(lat))) ** 2 + 1) ** -0.5
    return deg * r


def deg2my(deg):
    """Converts a latitudinal measurement from degrees to meters."""
    return 40007863 * deg / 360


def safe_divide(a, b):
    """Safely divide two arrays, with 0 as a result of a division by 0."""
    with np.errstate(divide="ignore", invalid="ignore"):
        c = np.true_divide(a, b)
        c[c == np.inf] = 0
        c = np.nan_to_num(c)
    return c


def geotransform(X, Y, GT):
    # From https://gdal.org/tutorials/geotransforms_tut.html
    X_geo = GT[0] + (X + 0.5) * GT[1] + (Y + 0.5) * GT[2]
    Y_geo = GT[3] + (X + 0.5) * GT[4] + (Y + 0.5) * GT[5]
    return X_geo, Y_geo
