#!/usr/bin/env python3
#
# Preprocessing for Illumina
#
# Author : Alexandre Simoneau
# modified by Hector Linares for AOD definition
#
# March 2021

import os
import shutil
from glob import glob

import click
import illum
import numpy as np
import yaml
from illum import MultiScaleData as MSD
from illum import pytools as pt
from illum.inventory import from_lamps, from_zones
from illum.OPAC import OPAC
from scipy.interpolate import griddata


@click.command(name="inputs")
def CLI_inputs():
    """Prepares the executions inputs."""
    inputs()


def inputs():

    print("Preparing the inputs for the experiment.")

    dir_name = "Inputs/"
    shutil.rmtree(dir_name, True)
    os.makedirs(dir_name)
    shutil.copy("inputs_params.in", dir_name + "inputs_params.in")

    with open("inputs_params.in") as f:
        params = yaml.safe_load(f)

    if (
        params["zones_inventory"] is not None
        and params["lamps_inventory"] is not None
    ):

        print("Validating the inventories.")

        lamps = np.loadtxt(params["lamps_inventory"], usecols=[0, 1])
        zones = np.loadtxt(params["zones_inventory"], usecols=[0, 1, 2])
        zonData = pt.parse_inventory(params["zones_inventory"], 7)

        hasLights = [sum(x[0] for x in z) != 0 for z in zonData]

        circles = MSD.from_domain("domain.ini")
        for dat, b in zip(zones, hasLights):
            circles.set_circle((dat[0], dat[1]), dat[2] * 1000, b)

        zones_ind = MSD.from_domain("domain.ini")
        for i, dat in enumerate(zones, 1):
            zones_ind.set_circle((dat[0], dat[1]), dat[2] * 1000, i)

        failed = set()
        for j, coords in enumerate(lamps, 1):
            for i in range(len(circles)):
                try:
                    col, row = circles._get_col_row(coords, i)
                    if circles[i][row, col] and col >= 0 and row >= 0:
                        zon_ind = zones_ind[i][row, col]
                        failed.add((j, coords[0], coords[1], zon_ind))
                except IndexError:
                    continue

        if len(failed):
            for i, lat, lon, zon_ind in sorted(failed):
                print(
                    "WARNING: Lamp #%d (%.06g,%.06g) falls within non-null zone #%d"
                    % (i, lat, lon, zon_ind)
                )
            raise SystemExit()

    out_name = params["exp_name"]

    if params["road_orientation"]:
        print("Computing road orientation (Can be slow for large domains)")
        from illum.street_orientation import street_orientation

        with open("domain.ini") as f:
            domain_params = yaml.safe_load(f)
        srs = domain_params["srs"]
        lats, lons = MSD.from_domain("domain.ini").get_obs_pos()
        bearings = street_orientation(lats, lons, srs)
        np.savetxt(dir_name + "/brng.lst", bearings, fmt="%g")

    print("Loading photometry files.")

    # Angular distribution (normalised to 1)
    lop_files = glob("Lights/*.lop")
    angles = np.arange(181, dtype=float)
    lop = {
        os.path.basename(s)
        .rsplit(".", 1)[0]
        .split("_", 1)[0]: pt.load_lop(angles, s)
        for s in lop_files
    }

    # Spectral distribution (normalised with scotopric vision to 1)
    wav, viirs = np.loadtxt("Lights/viirs.dat", skiprows=1).T
    viirs /= np.max(viirs)
    wav, norm_spectrum = np.loadtxt("Lights/photopic.dat", skiprows=1).T
    norm_spectrum /= np.max(norm_spectrum)

    spct_files = glob("Lights/*.spct")
    spct = {
        os.path.basename(s)
        .rsplit(".", 1)[0]
        .split("_", 1)[0]: pt.load_spct(wav, norm_spectrum, s)
        for s in spct_files
    }

    print("Splitting in wavelengths bins.")

    if os.path.isfile("spectral_bands.dat"):
        bins = np.loadtxt("spectral_bands.dat", delimiter=",")
        n_bins = bins.shape[0]
    else:
        n_bins = params["nb_bins"]
        lmin = params["lambda_min"]
        lmax = params["lambda_max"]

        limits = np.linspace(lmin, lmax, n_bins + 1)
        bins = np.stack([limits[:-1], limits[1:]], axis=1)

    bool_array = (wav >= bins[:, 0:1]) & (wav < bins[:, 1:2])
    x = bins.mean(1)
    bw = bins[:, 1] - bins[:, 0]

    print("Interpolating reflectance.")

    aster_files = glob("Lights/*.aster")
    aster = {
        os.path.basename(s).split(".", 1)[0]: np.loadtxt(s)
        for s in aster_files
    }

    for type in aster:
        wl, refl = aster[type].T
        wl *= 1000.0
        refl /= 100.0
        aster[type] = np.interp(wav, wl, refl)

    sum_coeffs = sum(
        params["reflectance"][type] for type in params["reflectance"]
    )
    if sum_coeffs == 0:
        sum_coeffs = 1.0

    refl = sum(
        aster[type] * coeff / sum_coeffs
        for type, coeff in params["reflectance"].items()
    )

    reflect = [np.mean(refl[mask]) for mask in bool_array]

    with open(dir_name + "/refl.lst", "w") as zfile:
        zfile.write("\n".join(["%.06g" % n for n in reflect]) + "\n")

    print("Linking mie files.")

    illumpath = os.path.dirname(illum.__path__[0])
    shutil.copy2(
        os.path.abspath(illumpath + "/Molecular_optics/MolecularAbs.txt"),
        dir_name,
    )

    OPAC(x)

    shutil.copy("srtm.hdf5", dir_name)

    with open(dir_name + "/wav.lst", "w") as zfile:
        zfile.write("".join("%g %g\n" % (w, b) for w, b in zip(x, bw)))

    if params["zones_inventory"] is not None:
        dir_name = ".Inputs_zones/"
        inv_name = params["zones_inventory"]
        n_inv = 7
        shutil.rmtree(dir_name, True)
        os.makedirs(dir_name)
        from_zones(
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
        )

    if params["lamps_inventory"] is not None:
        dir_name = ".Inputs_lamps/"
        shutil.rmtree(dir_name, True)
        os.makedirs(dir_name)
        from_lamps(
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
        )
    dir_name = "Inputs/"

    print("Unifying inputs.")

    lfiles = {fname.split(os.sep)[-1] for fname in glob(".Inputs_lamps/*")}
    zfiles = {fname.split(os.sep)[-1] for fname in glob(".Inputs_zones/*")}
    for fname in lfiles - zfiles:
        shutil.move(os.path.join(".Inputs_lamps", fname), "Inputs")
    for fname in zfiles - lfiles:
        shutil.move(os.path.join(".Inputs_zones", fname), "Inputs")
    for fname in zfiles & lfiles:
        if "fctem" in fname:
            shutil.move(os.path.join(".Inputs_lamps", fname), "Inputs")
        elif fname.endswith(".lst"):
            with open(os.path.join(".Inputs_lamps", fname)) as f:
                ldat = f.readlines()
            with open(os.path.join(".Inputs_zones", fname)) as f:
                zdat = f.readlines()
            with open(os.path.join("Inputs", fname), "w") as f:
                f.write("".join(sorted(set(ldat + zdat))))
        elif fname.endswith(".hdf5"):
            ldat = MSD.Open(os.path.join(".Inputs_lamps", fname))
            zdat = MSD.Open(os.path.join(".Inputs_zones", fname))
            for i, dat in enumerate(ldat):
                zdat[i][dat != 0] = dat[dat != 0]
            zdat.save(os.path.join("Inputs", fname))
        else:
            print("WARNING: File %s not merged properly." % fname)
    if "origin.hdf5" not in zfiles:
        origin = MSD.from_domain("domain.ini")
        origin.save("Inputs/origin")
    shutil.rmtree(".Inputs_lamps", True)
    shutil.rmtree(".Inputs_zones", True)

    # Interpolation of the obstacles properties
    defined = MSD.Open(dir_name + "origin.hdf5")
    lights_file = dir_name + out_name + "_lights.hdf5"
    if os.path.isfile(lights_file):
        lights = MSD.Open(lights_file)
        for i, layer in enumerate(lights):
            defined[i] += layer

    for geo in ["obsth", "obstd", "obstf", "altlp"]:
        geometry = MSD.Open(dir_name + out_name + "_" + geo + ".hdf5")
        for i, mask in enumerate(defined):
            geometry[i] = (
                griddata(
                    points=np.where(mask),
                    values=geometry[i][mask.astype(bool)],
                    xi=tuple(np.ogrid[0 : mask.shape[0], 0 : mask.shape[1]]),
                    method="nearest",
                )
                if mask.any()
                else np.zeros_like(geometry[i])
            )
        geometry.save(dir_name + out_name + "_" + geo)

    print("Done.")
