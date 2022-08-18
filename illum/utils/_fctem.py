# Functions related to spectral and angular emission functions
# Author: Alexandre Simoneau

import numpy as np
import scipy.interpolate

import illum.utils


def LOP_norm(angles, x):
    """Normalises 'x' as a function of theta over the full sphere.
    Uses the two first elements of 'angles' as the integration step.
    `angles` must be in degrees."""
    a = np.deg2rad(angles)
    mids = np.concatenate([[a[0]], np.mean([a[1:], a[:-1]], 0), [a[-1]]])
    sinx = 2 * np.pi * (np.cos(mids[:-1]) - np.cos(mids[1:]))
    return illum.pytools.safe_divide(x, np.sum(x * sinx))


def SPD_norm(wav, norm_spct, x, factor=683.002):
    """Normalises a spectrum 'x' with a normalisation spectrum to 'factor'.
    Both arrays must be of the same lenght.
    Uses the two first elements of 'wav' as the integration step."""
    dlambda = wav[1] - wav[0]
    return illum.pytools.safe_divide(
        x, factor * np.sum(norm_spct * x) * dlambda
    )


def spct_norm(wav, x):
    """Normalises 'x' using the two first elements of 'wav'as the integration
    step."""
    dlambda = wav[1] - wav[0]
    return illum.pytools.safe_divide(x, np.sum(x) * dlambda)


def load_lop(angles, filename, interp="cubic"):
    """Load an LOP file interpolated to 'angles' and normalised.

      interp : Interpolation kind

    See 'scipy.interpolate.interp1d for interpolation kinds."""
    data = np.loadtxt(filename).T
    y = (
        data[0]
        if np.all(data[1] == angles)
        else scipy.interpolate.interp1d(
            data[1], data[0], kind=interp, bounds_error=False, fill_value=0.0
        )
    )
    return LOP_norm(angles, y)


def load_spct(wav, norm_spct, filename, interp="cubic", factor=683.002):
    """Load a spectrum file interpolated to 'wav' and normalised.

      interp : Interpolation kind
      factor : Normalisation factor

    See 'scipy.interpolate.interp1d for interpolation kinds."""
    data = np.loadtxt(filename, skiprows=1).T
    y = (
        data[1]
        if np.all(data[0] == wav)
        else scipy.interpolate.interp1d(
            data[0], data[1], kind=interp, bounds_error=False, fill_value=0.0
        )(wav)
    )
    return SPD_norm(wav, norm_spct, y, factor)
