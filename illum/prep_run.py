#!/usr/bin/env python3

import os
import shutil
import time
from glob import glob

import numpy as np
import pandas as pd
import rasterio as rio
import rasterio.merge
import yaml
from tqdm import tqdm

import illum


def ground_type(dist, bounds):
    ground_type = np.ones(dist.shape) * len(bounds)
    for i, lim in reversed(list(enumerate(sorted(bounds)))):
        ground_type[dist < lim] = i
    return ground_type


def warp_like(src, dst, **kwargs):
    if type(src) == str:
        srcs = [rio.open(f) for f in glob(src)]
        if not srcs:
            raise ValueError("No SRTM files found.")
        src_data, src_transform = rio.merge.merge(srcs)
        src_crs = srcs[0].crs
    elif type(src) == rio.io.DatasetReader:
        src_data = src.read()
        src_transform = src.transform
        src_crs = src.crs

    out = np.zeros_like(dst.read(1))
    return rio.warp.reproject(
        source=src_data[0],
        destination=out,
        src_transform=src_transform,
        src_crs=src_crs,
        dst_transform=dst.transform,
        dst_crs=dst.crs,
        **kwargs,
    )[0]


if __name__ == "__main__":
    with open("iss_params.in") as f:
        p = yaml.safe_load(f)

    if os.path.exists("run"):
        for i in range(5, 0, -1):
            print(
                f"\rWARNING: Directory already exists. It will be deleted in {i}s.",
                end="",
            )
            time.sleep(1)
        print("\rDirectory deleted." + " " * 42)
        shutil.rmtree("run")
    os.makedirs("run")

    edges = np.linspace(
        p["wavelengths"]["min"],
        p["wavelengths"]["max"],
        p["wavelengths"]["nb"] + 1,
    )
    wls = np.average([edges[:-1], edges[1:]], axis=0)

    reflectances = {
        os.path.basename(fname)
        .split(os.path.extsep)[0]: illum.SPD.from_aster(fname)
        .interpolate(wls)
        for fname in glob("Lights/*.aster")
    }

    photo = illum.SPD.from_txt("Lights/photopic.dat").normalize()
    spcts = {
        os.path.basename(fname)
        .split(os.path.extsep)[0]
        .split("_")[0]: illum.SPD.from_txt(fname)
        .normalize(photo)
        .interpolate(wls)
        for fname in glob("Lights/*.spct")
    }
    spcts = pd.DataFrame(spcts.values(), index=spcts.keys())["data"].apply(pd.Series).T

    lops = {
        os.path.basename(fname)
        .split(os.path.extsep)[0]
        .split("_")[0]: illum.APD.from_txt(fname)
        .normalize()
        .interpolate(step=1)
        for fname in glob("Lights/*.lop")
    }

    for i, lop in enumerate(pd.DataFrame(p["inventory"])["lop"].values, 1):
        np.savetxt(
            f"run/fctem_{i}.dat",
            lops[str(lop)].data.T[180:, ::-1],
            fmt="%.15f",
            delimiter="\n",
        )

    dst_dist = rio.open("distance.tiff")
    dist = dst_dist.read(1)
    dst_tech = rio.open("tech.tiff")

    pts = pd.read_csv("sources.csv")
    coords = rio.transform.rowcol(dst_tech.transform, pts["lon"], pts["lat"])

    tech = np.zeros_like(dist, dtype=np.uint8)
    tech[coords] = pts["tech"] + 1

    illum.utils.save_bin("run/topogra.bin", warp_like("SRTM/*.hgt", dst_dist))
    illum.utils.save_bin("run/obsth.bin", np.full(dist.shape, p["height"]["obstacles"]))
    illum.utils.save_bin("run/altlp.bin", np.full(dist.shape, p["height"]["lamps"]))
    illum.utils.save_bin("run/azimu.bin", rio.open("bearing.tiff").read(1))
    illum.utils.save_bin("run/lmpty.bin", tech)
    illum.utils.save_bin("run/gndty.bin", ground_type(dist, p["ground_distances"]))

    with open("batch", "w") as f:
        pass

    for i, wl in tqdm(enumerate(wls), total=len(wls)):
        fold = f"run/band_{i+1:0{len(f'{len(wls)+1:d}')}d}"
        os.makedirs(fold)

        refls = (
            str(reflectances[p["refl"][r]].data[i])
            for r in ["road", "frontyard", "frontwall", "roof", "backwall", "backyard"]
        )

        with open(f"{fold}/illum-health.in", "w") as f:
            f.write("# ILLUMINA HEALTH INPUT PARAMETERS\n")
            f.write(f"{p['resolution']},{p['resolution']}\n")
            f.write(f"{dist.shape[1]},{dist.shape[0]}\n")
            f.write(
                f"{p['aerosols']['aod']},"
                f"{p['aerosols']['alpha']},"
                f"{p['aerosols']['scale']}\n"
            )
            f.write(f"{str(wl)},{wls[1]-wls[0]}\n")
            f.write(f"{','.join(refls)}\n")
            f.write(f"{p['aerosols']['pressure']}\n")
            f.write(f"{p['height']['observer']}\n")

        for fname in glob("run/*.bin") + glob("run/*.dat"):
            os.symlink(
                os.path.relpath(fname, start=fold),
                os.path.join(fold, os.path.basename(fname)),
            )

        os.symlink(
            f"{illum.path}/data/Molecular_optics/MolecularAbs.txt",
            f"{fold}/MolecularAbs.txt",
        )

        lumlp = np.zeros_like(dist)
        lumlp[coords] = pts["pow"] * spcts.iloc[i][pts["spct"]].to_numpy()
        illum.utils.save_bin(f"{fold}/lumlp.bin", lumlp)

        with open("batch", "a") as f:
            f.write(f"cd {os.path.abspath(fold)}\n")
            f.write(os.path.abspath(f"{illum.path}/../bin/illum-health") + "\n")
