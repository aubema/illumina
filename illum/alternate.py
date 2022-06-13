#!/usr/bin/env python3

import os
import shutil
from glob import glob

import click
import numpy as np
import yaml
from illum import MultiScaleData as MSD
from illum import pytools as pt
from illum.inventory import from_lamps, from_zones
from scipy.interpolate import interp1d as interp


@click.command()
@click.argument("name")
@click.option(
    "-z",
    "--zones",
    type=click.Path(exists=True),
    help="New zones inventory filename.",
)
@click.option(
    "-l",
    "--lights",
    type=click.Path(exists=True),
    help="New discrete lights inventory filename.",
)
def CLI_alternate(name, zones, lights):
    """Generates an alternate scenario at constant lumen.

    This scenatio will be based on the content of the `Inputs` folder and
    will be placed in a folder named `Inputs_NAME`.
    """
    alternate(name, zones, lights)


def alternate(name, zones, lights):
    if zones is None and lights is None:
        print("ERROR: At least one of 'zones' and 'lights' must be provided.")
        raise SystemExit

    dirname = "Inputs_%s/" % name

    if os.path.exists(dirname):
        shutil.rmtree(dirname)
    os.makedirs(dirname)

    with open("inputs_params.in") as f:
        params = yaml.safe_load(f)

    if zones is not None and lights is not None:

        print("Validating the inventories.")

        lamps = np.loadtxt(lights, usecols=[0, 1])
        zones = np.loadtxt(params["zones_inventory"], usecols=[0, 1, 2])
        zonData = pt.parse_inventory(zones, 0)

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

    shutil.copy("Inputs/inputs_params.in", dirname)

    print("\nLoading data...")

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
    viirs = pt.spct_norm(wav, viirs)
    scotopic = pt.load_spct(wav, np.ones(wav.shape), "Lights/scotopic.dat", 1)
    photopic = pt.load_spct(wav, np.ones(wav.shape), "Lights/photopic.dat", 1)

    ratio_ps = 1.0
    norm_spectrum = ratio_ps * photopic + (1 - ratio_ps) * scotopic

    spct_files = glob("Lights/*.spct")
    spct = {
        os.path.basename(s)
        .rsplit(".", 1)[0]
        .split("_", 1)[0]: pt.load_spct(wav, norm_spectrum, s)
        for s in spct_files
    }

    # Make bins
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
    x = bins.mean(1).tolist()

    out_name = params["exp_name"]

    asper_files = glob("Lights/*.aster")
    asper = {
        os.path.basename(s).split(".", 1)[0]: np.loadtxt(s)
        for s in asper_files
    }

    for type in asper:
        wl, refl = asper[type].T
        wl *= 1000.0
        refl /= 100.0
        asper[type] = interp(wl, refl, bounds_error=False, fill_value=0.0)(wav)

    sum_coeffs = sum(
        params["reflectance"][type] for type in params["reflectance"]
    )
    if sum_coeffs == 0:
        sum_coeffs = 1.0

    refl = sum(
        asper[type] * coeff / sum_coeffs
        for type, coeff in params["reflectance"].items()
    )

    reflect = [np.mean(refl[mask]) for mask in bool_array]

    with open(dirname + "/refl.lst", "w") as zfile:
        zfile.write("\n".join(["%.06g" % n for n in reflect]) + "\n")

    # Photopic/scotopic spectrum
    nspct = ratio_ps * photopic + (1 - ratio_ps) * scotopic
    nspct = nspct / np.sum(nspct)
    nspct = [np.mean(nspct[mask]) for mask in bool_array]

    for aero_file in glob("Inputs/*.txt"):
        shutil.copy(aero_file, aero_file.replace("Inputs", dirname))

    shutil.copy("srtm.hdf5", dirname)

    with open(dirname + "/wav.lst", "w") as zfile:
        zfile.write("\n".join(map(str, x)) + "\n")

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

        oldlumlp = MSD.from_domain("domain.ini")
        for fname in glob("Inputs/*lumlp*"):
            ds = MSD.Open(fname)
            wl = int(fname.split("_")[1])
            for i, dat in enumerate(ds):
                oldlumlp[i] += dat * nspct[x.index(wl)] * dl

        newlumlp = MSD.from_domain("domain.ini")
        for fname in glob(os.path.join(dir_name, "*lumlp*")):
            ds = MSD.Open(fname)
            wl = int(fname.split("_")[2])
            for i, dat in enumerate(ds):
                newlumlp[i] += dat * nspct[x.index(wl)] * dl

        ratio = MSD.from_domain("domain.ini")
        for i in range(len(ratio)):
            ratio[i] = pt.safe_divide(oldlumlp[i], newlumlp[i])

        for fname in glob(os.path.join(dir_name, "*lumlp*")):
            ds = MSD.Open(fname)
            for i, dat in enumerate(ratio):
                ds[i] *= dat
            ds.save(fname)

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

    print("Unifying inputs.")

    lfiles = {fname.split(os.sep)[-1] for fname in glob(".Inputs_lamps/*")}
    zfiles = {fname.split(os.sep)[-1] for fname in glob(".Inputs_zones/*")}
    for fname in lfiles - zfiles:
        shutil.move(os.path.join(".Inputs_lamps", fname), dirname)
    for fname in zfiles - lfiles:
        shutil.move(os.path.join(".Inputs_zones", fname), dirname)
    for fname in zfiles & lfiles:
        if "fctem" in fname:
            shutil.move(os.path.join(".Inputs_lamps", fname), dirname)
        elif fname.endswith(".lst"):
            with open(os.path.join(".Inputs_lamps", fname)) as f:
                ldat = f.readlines()
            with open(os.path.join(".Inputs_zones", fname)) as f:
                zdat = f.readlines()
            with open(os.path.join(dirname, fname), "w") as f:
                f.write("".join(sorted(set(ldat + zdat))))
        elif fname.endswith(".hdf5"):
            ldat = MSD.Open(os.path.join(".Inputs_lamps", fname))
            zdat = MSD.Open(os.path.join(".Inputs_zones", fname))
            for i, dat in enumerate(ldat):
                zdat[i][dat != 0] = dat[dat != 0]
            zdat.save(os.path.join(dirname, fname))
        else:
            print("WARNING: File %s not merged properly." % fname)
    if "origin.hdf5" not in zfiles:
        origin = MSD.from_domain("domain.ini")
        origin.save(dirname + "/origin")
    shutil.rmtree(".Inputs_lamps", True)
    shutil.rmtree(".Inputs_zones", True)

    print("Done.")
