#!/usr/bin/env python3
#
# Inventory preprocessing for Illumina
#
# Author : Alexandre Simoneau
#
# December 2021

import numpy as np

import illum.pytools as pt
from illum import MultiScaleData as MSD


def from_lamps(
    dir_name,
    n_bins,
    params,
    out_name,
    x,
    lop,
    angles,
    wav,
    spct,
    viirs,
    refl,
    bool_array,
):

    print("Building inputs from discrete inventory.")

    # lamps distribution
    inv_name = params["lamps_inventory"]
    lampsData = np.loadtxt(inv_name, usecols=list(range(7)), ndmin=2)
    photometry = np.loadtxt(inv_name, usecols=[-2, -1], dtype=str, ndmin=2)
    domain = MSD.from_domain("domain.ini")

    sources = np.unique(photometry[:, 1])

    print("Classifying points.")

    points = dict()
    for layer in range(len(domain)):
        ysize, xsize = domain[layer].shape
        col, row = domain._get_col_row(lampsData[:, :2].T, layer)
        valid = (0 <= col) * (col < xsize) * (0 <= row) * (row < ysize) > 0
        ind = np.where(valid)[0]
        points[layer] = (col[ind], row[ind], ind)

    print("Calculating the generalized lamps.")

    for n in range(n_bins):
        for s in sources:
            profile = lop[s].vertical_profile()[::-1]
            np.savetxt(
                dir_name + f"fctem_wl_{x[n]:g}_lamp_{s}.dat",
                np.concatenate([profile, angles]).reshape((2, -1)).T,
            )

    with open(dir_name + "lamps.lst", "w") as zfile:
        zfile.write("\n".join(sources) + "\n")

    geometry = dict()
    for geo in ["obsth", "obstd", "obstf", "altlp", "lights"]:
        geometry[geo] = MSD.from_domain("domain.ini")

    lumlp = dict()
    for s in sources:
        for wl in x:
            lumlp[s, wl] = MSD.from_domain("domain.ini")

    for layer, pts in points.items():
        cols, rows, inds = pts
        if len(inds):
            for col, row in np.unique([cols, rows], axis=1).T:
                ind = inds[np.logical_and(cols == col, rows == row)]
                lumens = lampsData[:, 2][ind]

                for n, geo in zip(range(3, 7), ["obsth", "obstd", "obstf", "altlp"]):
                    geometry[geo][layer][row, col] = np.average(
                        lampsData[:, n][ind], weights=lumens
                    )
                geometry["lights"][layer][row, col] = 1

                local_sources = np.unique(photometry[ind][:, 1])
                for s in local_sources:
                    mask = photometry[:, 1][ind] == s
                    fctem = np.array(
                        [spct[type].data for type in photometry[:, 0][ind][mask]]
                    )
                    fctem = np.sum(fctem * lumens[mask, None], 0)

                    y = [np.mean(fctem[mask]) for mask in bool_array]

                    for i, wl in enumerate(x):
                        lumlp[s, wl][layer][row, col] = y[i]

    print("Saving data.")

    for geo, ds in geometry.items():
        ds.save(dir_name + out_name + "_" + geo)

    for key, ds in lumlp.items():
        s, wl = key
        ds.save(dir_name + f"{out_name}_{wl:g}_lumlp_{s}")


def from_zones(
    dir_name,
    inv_name,
    n_inv,
    n_bins,
    params,
    out_name,
    x,
    lop,
    angles,
    wav,
    spct,
    viirs,
    refl,
    bool_array,
):
    print("Building inputs from zones inventory.")

    # lamps distribution
    zonData = pt.parse_inventory(inv_name, n_inv)

    sources = np.unique([lamp[2] for zd in zonData for lamp in zd])

    for n in range(n_bins):
        for s in sources:
            profile = lop[s].vertical_profile()[::-1]
            np.savetxt(
                dir_name + f"fctem_wl_{x[n]:g}_lamp_{s}.dat",
                np.concatenate([profile, angles]).reshape((2, -1)).T,
            )

    with open(dir_name + "lamps.lst", "w") as zfile:
        zfile.write("\n".join(sources) + "\n")

    print("Making zone properties files.")

    circles = MSD.from_domain("domain.ini")  # Same geolocalisation
    zonfile = np.loadtxt(params["zones_inventory"], usecols=list(range(7)), ndmin=2)

    # zone number
    for i, dat in enumerate(zonfile, 1):
        circles.set_circle((dat[0], dat[1]), dat[2] * 1000, i)
    circles.save(dir_name + out_name + "_zone")

    weights = [sum(z[0] for z in zone) for zone in zonData]
    for w, dat in zip(weights, zonfile):
        circles.set_circle((dat[0], dat[1]), dat[2] * 1000, bool(w))
    circles.save(dir_name + "origin")

    for n, name in zip(range(3, 7), ["obsth", "obstd", "obstf", "altlp"]):
        for i, dat in enumerate(zonfile, 1):
            circles.set_circle((dat[0], dat[1]), dat[2] * 1000, dat[n])
        circles.save(dir_name + out_name + "_" + name)

    print("Inverting lamp intensity.")

    viirs_dat = MSD.Open("stable_lights.hdf5")
    for i in range(len(viirs_dat)):
        viirs_dat[i] *= 1e-5  # nW/cm^2/sr -> W/m^2/sr
        viirs_dat[i][viirs_dat[i] < 0] = 0.0

    water_mask = MSD.Open("water_mask.hdf5")
    for i, wm in enumerate(water_mask):
        viirs_dat[i][wm == 0] = 0.0

    circles = MSD.Open(dir_name + out_name + "_zone.hdf5")
    zon_mask = np.empty(len(circles), dtype=object)
    for i in range(len(zon_mask)):
        zon_mask[i] = np.arange(1, len(zonfile) + 1)[:, None, None] == circles[i]

    a = np.deg2rad(angles)
    mids = np.concatenate([[a[0]], np.mean([a[1:], a[:-1]], 0), [a[-1]]])
    sinx = 2 * np.pi * (np.cos(mids[:-1]) - np.cos(mids[1:]))

    # Pixel size in m^2
    S = np.array([viirs_dat.pixel_size(i) ** 2 for i in range(len(viirs_dat))])

    # Calculate zones lamps
    zones = pt.make_zones(angles, lop, wav, spct, zonData, sources).transpose(
        (0, 1, 3, 2)
    )

    # phie = DNB * S / int( R ( rho/pi Gdown + Gup ) ) dlambda
    Gdown = np.dot(zones[..., angles > 90], sinx[angles > 90])
    Gup = np.dot(zones[..., angles < 70], sinx[angles < 70]) / sinx[angles < 70].sum()
    integral = np.sum(viirs.data * (Gdown * refl / np.pi + Gup), (1, 2)) * (
        wav[1] - wav[0]
    )

    phie = [
        pt.safe_divide(
            viirs_dat[i] * S[i],
            np.sum(zon_mask[i] * integral[:, None, None], 0),
        )
        for i in range(len(S))
    ]

    ratio = [np.dot(zones[:, :, ind], sinx).mean(-1) for ind in bool_array]

    for n in range(n_bins):
        r = [
            np.sum(zon_mask[layer][:, None] * ratio[n][:, :, None, None], 0)
            for layer in range(len(phie))
        ]
        for i, s in enumerate(sources):
            new = MSD.from_domain("domain.ini")
            for layer in range(len(new)):
                new[layer] = phie[layer] * r[layer][i]
            new.save(dir_name + f"{out_name}_{x[n]:g}_lumlp_{s}")
