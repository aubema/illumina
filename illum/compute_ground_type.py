#!/usr/bin/env python3

import numpy as np
import pandas as pd
import yaml
from osgeo import gdal


def open_tiff(filename, dtype=np.float32):
    # Load file, and access the band and get a NumPy array
    src = gdal.Open(filename, gdal.GA_Update)
    band = src.GetRasterBand(1)
    ar = band.ReadAsArray()
    return src, ar


def save_geotiff(filename, data, src):
    nband = 1
    nrow, ncol = data.shape
    driver = gdal.GetDriverByName("GTiff")
    dst_dataset = driver.Create(
        filename + ".tiff", ncol, nrow, nband, gdal.GDT_Float32
    )
    # sets same geotransform as input
    dst_dataset.SetGeoTransform(src.GetGeoTransform())
    # sets same projection as input
    dst_dataset.SetProjection(src.GetProjection())
    dst_dataset.GetRasterBand(1).WriteArray(data.astype(float))
    dst_dataset = None


with open("iss_params.in") as f:
    p = yaml.safe_load(f)

src, im = open_tiff(p["wd"] + "/domain.tiff")
df = pd.read_csv(p["wd"] + "/obs.csv")

dists = (10, 15, 25, 30)
ground_type = np.ones(im.shape) * len(dists)
for i, lim in reversed(list(enumerate(dists))):
    mask = df.distance < lim
    ground_type[df.Y[mask], df.X[mask]] = i

save_geotiff(p["wd"] + "/ground_type", ground_type, src)
