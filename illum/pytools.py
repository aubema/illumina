#!/usr/bin/env python3
#
# Library of usefull function related to Illumina and PGM and FITS handling
#
# Author : Alexandre Simoneau
# unless noted otherwise
#
# March 2016

import astropy.io.fits as _fits
import matplotlib.colors as _colors
import matplotlib.pyplot as _plt
import numpy as _np
import scipy.interpolate as _I


def safe_divide(a, b):
    """Safely divide two arrays, with 0 as a result of a division by 0."""
    with _np.errstate(divide="ignore", invalid="ignore"):
        c = _np.true_divide(a, b)
        c[c == _np.inf] = 0
        c = _np.nan_to_num(c)
    return c


def LOP_norm(angles, x):
    """Normalises 'x' as a function of theta over the full sphere.
    Uses the two first elements of 'angles' as the integration step.
    `angles` must be in degrees."""
    a = _np.deg2rad(angles)
    mids = _np.concatenate([[a[0]], _np.mean([a[1:], a[:-1]], 0), [a[-1]]])
    sinx = 2 * _np.pi * (_np.cos(mids[:-1]) - _np.cos(mids[1:]))
    return safe_divide(x, _np.sum(x * sinx))


def SPD_norm(wav, norm_spct, x, factor=683.002):
    """Normalises a spectrum 'x' with a normalisation spectrum to 'factor'.
    Both arrays must be of the same lenght.
    Uses the two first elements of 'wav' as the integration step."""
    dlambda = wav[1] - wav[0]
    return safe_divide(x, factor * _np.sum(norm_spct * x) * dlambda)


def spct_norm(wav, x):
    """Normalises 'x' using the two first elements of 'wav'as the integration
    step."""
    dlambda = wav[1] - wav[0]
    return safe_divide(x, _np.sum(x) * dlambda)


def zon_norm(angles, wavelenght, zone):
    """Returns the normalisation factor of an ILLUMINA zone."""
    a = _np.deg2rad(angles)
    mids = _np.concatenate([[a[0]], _np.mean([a[1:], a[:-1]], 0), [a[-1]]])
    sinx = 2 * _np.pi * (_np.cos(mids[:-1]) - _np.cos(mids[1:]))
    dlambda = wavelenght[1] - wavelenght[0]
    return _np.sum(zone.T * sinx) * dlambda


def parse_inventory(filename, n=0):
    """Parse an inventory type file.
    Skips the first 'n' columns."""

    def lamp_norm(lampsData):
        trans = list(map(list, list(zip(*lampsData))))
        norm = sum(trans[0])
        if norm != 0.0:
            trans[0] = [n / norm for n in trans[0]]
        return list(map(list, list(zip(*trans))))

    with open(filename) as inv_file:
        zonData = strip_comments(inv_file.readlines())
    zonData = [s.split()[n:] for s in zonData]
    zonData = [[s.split("_") for s in i] for i in zonData]
    zonData = [[[float(s[0]), s[1], s[2]] for s in i] for i in zonData]
    zonData = list(map(lamp_norm, zonData))

    return zonData


def make_zones(theta, lop, wl, spct, ivtr, sources):
    """Returns an array of normalized zones.

    theta : angles used to define the Light Output Pattern (deg)
    lop : Light Output Pattern dictionnary
    wl : wavelength used to define the spectrum
    spct : Lamp spectrum dictionnary
    ivtr : Parsed inventory
    lops : List of LOPs in the inventory"""

    zones = _np.ones((len(ivtr), len(sources), len(theta), len(wl)))
    for i, zone in enumerate(ivtr):
        for j, s in enumerate(sources):
            zones[i, j] *= sum(
                c * spct[t].data * lop[l].vertical_profile()[::-1, None]
                for c, t, l in zone
                if l == s
            )

    for i, zone in enumerate(zones):
        zon = zone.sum(0)
        norm = zon_norm(theta, wl, zon)
        zones[i] = safe_divide(zone, norm)

    return zones


def load_pgm(filename):
    """Opens a PGM file.

    Returns the header dictionnary, the image parameters and the data as a
    3-ple. The image parameters are a list of the height, the width and the
    resolution."""
    with open(filename) as file:
        data = file.read().split("\n")[:-1]

    head = (s for s in data if s[0] == "#")
    head = {s.split()[1]: s.split()[2] for s in head}

    gain = float(head.setdefault("gain", 0))
    offset = float(head.setdefault("offset", 0))

    header = data[: len(head) + 2]
    data = _np.loadtxt(filename, skiprows=len(head) + 2, ndmin=2)
    data = data * gain + offset

    p = list(map(int, header[-1].split()))

    return head, p, data.reshape(p[1::-1])


def save_pgm(filename, head, p, data, offset=0.0):
    """Saves a PGM file.

    head : A dictionnary of headers
    p : The image parameters as returned by 'load_pgm'
    data : Data to be saved
    offset : Either a value or 'min', the minimum of the data."""

    offset = _np.min(data) if offset == "min" else offset
    new_gain = (_np.max(data) - offset) / float(p[2])
    if new_gain != 0:
        data = (data - offset) / new_gain
        data = _np.round(data).astype(int)
        data[data < 0] = 0

    head["gain"] = str(new_gain)
    head["offset"] = str(offset)

    headstring = (
        "P2\n"
        + "\n".join(["# %s %s" % s for s in iter(list(head.items()))])
        + "\n"
        + " ".join(map(str, p))
    )
    _np.savetxt(filename, data, fmt="%d", header=headstring, comments="")


def load_fits(filename):
    """Loads a FITS file.

    Returns a 2-ple containing
      - A list of arrays defining each axis in order (x,y,z,...)
      - The contained data in transposed order (...,z,y,x)
    """
    hdu = _fits.open(filename)[0]

    ax = [
        _np.linspace(
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
    hdu = _fits.PrimaryHDU()
    hdu.data = data.T[:, ::-1].T
    for i in range(len(axis)):
        hdu.header["CRPIX%d" % (i + 1)] = 1
        hdu.header["CRVAL%d" % (i + 1)] = axis[i][0]
        hdu.header["CDELT%d" % (i + 1)] = axis[i][1]
    hdu.writeto(filename, clobber=True)


def load_bin(filename, dtype=_np.float32):
    """Load a ILLUMINA binary file.

    Returns the data as an array."""
    with open(filename) as f:
        shape = _np.fromfile(f, dtype=_np.uint32, count=4)[1:-1][::-1]
        data = _np.fromfile(f, dtype=_np.float32, count=-1)[1::3]
    return data.reshape(shape).astype(dtype)


def save_bin(filename, data):
    """Saves a numpy data array as an ILLUMINA binary file."""
    data = data.astype(_np.float32)
    shape = data.shape[::-1]
    size = data.size
    data_flat = data.flatten()
    filler = _np.ones(size, dtype=_np.float32) * 5.6e-45

    head = _np.array((8,) + shape + (8,)).astype(_np.uint32)
    body = _np.array([filler, data_flat, filler]).T.flatten()

    with open(filename, "w") as f:
        head.tofile(f)
    with open(filename, "a") as f:
        body.tofile(f)


def strip_comments(item, token="#"):
    """Generator. Strips comments and whitespace from input lines.

    This generator strips comments, leading/trailing whitespace, and
    blank lines from its input.

    Arguments:
        item (obj):  Object to strip comments from.
        token (str, optional):  Comment delimiter.  Defaults to ``#``.

    Yields:
        str:  Next non-blank line from ``item`` with comments and
            leading/trailing whitespace removed.

    Credits: Doug R., StackOverflow
    """

    for line in item:
        s = line.split(token, 1)[0].strip()
        if s != "":
            yield s


def load_lop(angles, filename, interp="cubic"):
    """Load an LOP file interpolated to 'angles' and normalised.

      interp : Interpolation kind

    See 'scipy.interpolate.interp1d for interpolation kinds."""
    data = _np.loadtxt(filename).T
    y = (
        data[0]
        if _np.all(data[1] == angles)
        else _I.interp1d(
            data[1], data[0], kind=interp, bounds_error=False, fill_value=0.0
        )
    )
    return LOP_norm(angles, y)


def load_spct(wav, norm_spct, filename, interp="cubic", factor=683.002):
    """Load a spectrum file interpolated to 'wav' and normalised.

      interp : Interpolation kind
      factor : Normalisation factor

    See 'scipy.interpolate.interp1d for interpolation kinds."""
    data = _np.loadtxt(filename, skiprows=1).T
    y = (
        data[1]
        if _np.all(data[0] == wav)
        else _I.interp1d(
            data[0], data[1], kind=interp, bounds_error=False, fill_value=0.0
        )(wav)
    )
    return SPD_norm(wav, norm_spct, y, factor)


def plot_allsky(phi, r, data, n=100, **kwargs):
    """Plot all sky data.

    Parameters:
      phi      : Azimuthal angle (deg). Must be linearly spaced.
      r        : Elevation angle (deg). Must be in [0;90]
      data     : (NxM) or (NxMx3) data to plot with N the number of elevation
                 angles and M the number of azimuthal angles. If 3D the Data
                 will be represented as an RGB plot.
      n        : Supersampling. If 'interp' is not 'None', it will be applied
                 to both axis. (Used to make the plot round)
      interp   : Interpolation method. Can be 'None', 'linear', 'cubic' or
                 'quintic'
      autogain : If True, the data will be stretched to have a max of 1
                 (def. False)
      title    : Figure title
      clabel   : Colorbar label. If defined, will draw a colorbar (2D data
                 only)
      fname    : File name. If defined, will save the figure and close it.
      cmap     : Colormap to use (2D data only). See matplotlib doc for
                 options.
      labels   : Labels to add to the figure. {'label':position (deg)}
                 dictionnary.
      log      : Logaritmic color scale (def. False)
      vmin     : Minimal value of the colorscale (2D data only).
      vmax     : Maximal value of the colorscale (2D data only).
      showpts  : If set to True, will show datapoints.
    """
    if "showpts" in kwargs and kwargs["showpts"]:
        xv, yv = _np.meshgrid(_np.deg2rad(phi), 90 - _np.array(r))
    if "interp" in kwargs and kwargs["interp"] != "None":
        ndata = _np.concatenate([data, data[:, 0, ...][:, None]], axis=1)
        nphi = _np.linspace(phi[0], phi[0] + 360, n * len(phi), endpoint=False)
        nr = _np.linspace(r[0], r[-1], n * (len(r) - 1) + 1)
        if data.ndim == 3:
            data = _np.zeros((len(nr), len(nphi), data.shape[2]))
            for i in range(data.shape[2]):
                interp = _I.interp2d(
                    _np.concatenate([phi, [phi[0] + 360]]),
                    r,
                    ndata[:, :, i],
                    kind=kwargs["interp"],
                )
                data[:, :, i] = interp(nphi, nr)
        else:
            interp = _I.interp2d(
                _np.concatenate([phi, [phi[0] + 360]]),
                r,
                ndata,
                kind=kwargs["interp"],
            )
            data = interp(nphi, nr)
        phi = nphi
        r = nr
        n = 1
    else:
        data = _np.repeat(data, n, axis=1)

    r = _np.asarray(r).tolist()
    r[0:0] = [2 * r[0] - r[1]]
    r[-1:-1] = [r[-1]]
    r = 90 - _np.mean([r[:-1], r[1:]], 0)
    r[r < 0] = 0
    r[r > 90] = 90

    theta = _np.linspace(0, 2 * _np.pi, n * len(phi) + 1) - _np.mean(
        _np.radians(phi[:2])
    )

    Theta, R = _np.meshgrid(theta, r)
    if "autogain" in kwargs and kwargs["autogain"]:
        gain = _np.max(data)
        data /= gain
    if data.ndim == 3:
        color = data.reshape((-1, data.shape[2]))

    _plt.figure()
    ax = _plt.subplot(111, polar=True)
    ax.set_theta_zero_location("N")
    ax.xaxis.set_ticklabels(["N", "NE", "E", "SE", "S", "SW", "W", "NW"])
    if data.ndim == 3:
        m = _plt.pcolormesh(
            Theta, R, data[:, :, 0], color=color, linewidth=0, vmin=0, vmax=1
        )
        m.set_array(None)
    else:
        args = dict()
        if "cmap" in kwargs:
            args["cmap"] = kwargs["cmap"]
        if "vmin" in kwargs:
            args["vmin"] = kwargs["vmin"]
        if "vmax" in kwargs:
            args["vmax"] = kwargs["vmax"]
        if "log" in kwargs and kwargs["log"]:
            args["norm"] = _colors.LogNorm(
                vmin=kwargs["vmin"], vmax=kwargs["vmax"]
            )
        m = _plt.pcolormesh(Theta, R, data, linewidth=0, **args)

    if "showpts" in kwargs and kwargs["showpts"]:
        _plt.plot(xv, yv, "r.")
    title_str = ""
    if "title" in kwargs:
        title_str += kwargs["title"]
    if "autogain" in kwargs and kwargs["autogain"]:
        if title_str != "":
            title_str += "\n"
        title_str += "gain = %.4g" % gain
    if title_str != "":
        _plt.title(title_str)

    if "clabel" in kwargs:
        _plt.colorbar(label=kwargs["clabel"], pad=0.11)

    if "labels" in kwargs:
        for name, pos in list(kwargs["labels"].items()):
            _plt.annotate(
                name,
                xy=[_np.radians(pos), _np.max(r)],
                xytext=[_np.radians(pos), _np.max(r) + 7],
                verticalalignment="top" if (90 < pos < 270) else "bottom",
                horizontalalignment="right" if pos < 180 else "left",
                annotation_clip=False,
                arrowprops=dict(facecolor="black", width=0.1, headlength=0.1),
            )

    _plt.tight_layout()

    if "fname" in kwargs:
        _plt.savefig(kwargs["fname"])
        _plt.close()
