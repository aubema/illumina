#!/usr/bin/env python3

from glob import glob

import numpy as np
import pandas as pd
import rasterio as rio
import yaml
from pykml import parser

import illum.pytools as pt

with open("iss_params.in") as f:
    p = yaml.safe_load(f)


exp_folds = glob("iss???e*/")
if not exp_folds:
    print("ERROR: Data folder not found. Please make sure it was not renamed.")
    exit()
if len(exp_folds) > 1:
    print(f"ERROR: Too many data folders found. {exp_folds}")
    exit()

exp_fold = exp_folds[0]
basename, ts = exp_fold.split("_")
print(f"Data `{basename}` ({ts[:4]}-{ts[4:6]}-{ts[6:8]}) found.")

with open(f"{exp_fold}/popular/{basename}RGBpos.kml") as f:
    obs_angle = float(
        str(parser.parse(f).getroot().getchildren()[0].description).split()[1]
    )

abc = pd.DataFrame(p["inventory"])
error = False
for spct in abc["spct"]:
    if not glob(f"Lights/{spct}_*.spct"):
        print(f"ERROR: `{spct}` spectral emission definition file not found.")
        error = True
for lop in abc["lop"]:
    if not glob(f"Lights/{lop}_*.lop"):
        print(f"ERROR: `{lop}` angular emission definition file not found.")
        error = True
if error:
    quit()

wav, norm_spectrum = np.loadtxt("Lights/photopic.dat", skiprows=1).T
norm_spectrum /= norm_spectrum.max()
wl, refl = np.loadtxt("Lights/asphalt.aster").T
refl = np.interp(wav, wl * 1000, refl / 100)
angles = np.arange(181, dtype=float)

# Load pre-normalize spcts & lops
spcts = np.array(
    [
        pt.load_spct(wav, norm_spectrum, glob(f"Lights/{spct}_*.spct")[0])
        for spct in abc["spct"]
    ]
)
lops = np.array(
    [pt.load_lop(angles, glob(f"Lights/{lop}_*.lop")[0]) for lop in abc["lop"]]
)

S = ((p["ISS_alt"] * p["pixel_size"]) / p["focal_dist"]) ** 2
a = np.deg2rad(angles)
mids = np.concatenate([[a[0]], np.mean([a[1:], a[:-1]], 0), [a[-1]]])
sinx = 2 * np.pi * (np.cos(mids[:-1]) - np.cos(mids[1:]))


# Calculating integral for each lop-spct combinaisons
# Atmospheric correction pre-calculated by Alejandro
# Intensity: (nw/sr/cm^2/angstrom) We suppose const value on photopic bandwidth
# Formula: I = DNB * S / integral( R(lambd) * T(lambd) *
# (1/pi* p(lambd)* F(lambd) + G(lambd))) dlambd)
integral = S / np.trapz(
    spcts.T
    * (
        (refl[:, None] / np.pi)
        * np.sum(lops[:, angles > 90] * sinx[angles > 90], axis=1)
        + lops[:, round(obs_angle)]
    )
    * norm_spectrum[:, None],
    x=wav,
    axis=0,
)

print("Creating discrete inventory..")
src = rio.open("tech.tiff")
im = src.read(1)
pts = pd.read_csv("peaks.csv")

Y, X = rio.transform.rowcol(src.transform, pts["lon"], pts["lat"])
pts["tech"] = im[Y, X].astype(np.uint8)
pts = pts[pts["tech"] > 0]
tech_idx = pts["tech"] - 1
pts["pow"] = integral[tech_idx] * pts.pop("flux") * p["bandwidth"] * 1e-4
pts["spct"] = abc["spct"][tech_idx].to_numpy()
pts["lop"] = abc["lop"][tech_idx].to_numpy()

pts.to_csv("sources.csv", index=False)
