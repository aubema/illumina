#!/usr/bin/env python

import argparse
import re

import numpy as np
import pyproj
import yaml
from matplotlib.path import Path

from illum.pytools import load_bin

parser = argparse.ArgumentParser(
    description="Integrates Illumina binary file over a polygon."
)
parser.add_argument("domain", help="Domain characteristics file [domain.ini].")
parser.add_argument("bin", help="Binary file to integrate.")
parser.add_argument(
    "kml",
    nargs="+",
    help="KML file or files that defines the area over wich to integrate.",
)

p = parser.parse_args()

bin = load_bin(p.bin)

# Load zones

regex = re.compile(r"<coordinates>\s*(.*)\s*<\/coordinates>")

zones = dict()
for kml_file in p.kml:
    with open(kml_file) as f:
        kml_data = f.read()

    coords_txt = re.findall(regex, kml_data)[0]
    coords_txt = coords_txt.strip().replace(" ", ";")
    coords = np.asarray(np.matrix(coords_txt))

    zones[kml_file.split(".", 1)[0]] = coords[:, :2]

# Project zone coordinates

with open(p.domain) as f:
    domain = yaml.load(f)

domain["xmin"], domain["ymin"], domain["xmax"], domain["ymax"] = list(
    map(float, domain["bbox"].split())
)

p1 = pyproj.CRS.from_epsg(4326)  # WGS84
p2 = pyproj.CRS.from_user_input(domain["srs"])
transform = pyproj.Transformer.from_crs(p1, p2, always_xy=True).transform

for zone, data in zones.items():
    lat, lon = data.T

    x, y = pyproj.transform(lon, lat)

    data[:, 0] = (x - domain["xmin"]) / domain["pixsize"] + 1
    data[:, 1] = (y - domain["ymin"]) / domain["pixsize"] + 1

# Define masks
nr, nc = bin.shape
ygrid, xgrid = np.mgrid[:nr, :nc]
xypix = np.vstack((xgrid.ravel(), ygrid.ravel())).T

for zone, data in zones.items():
    pth = Path(data, closed=False)
    mask = pth.contains_points(xypix)
    mask = mask.reshape(bin.shape)
    print(zone, np.sum(bin[mask]))
