#!/usr/bin/env python3
# ******************************************************************************
#        Determine street orientation from list of coordinates
#
# Authors : Julien-Pierre Houle & Alexandre Simoneau
# Date :    May 2021
# ******************************************************************************

import numpy as np
import osmnx as ox
import pyproj
import rasterio
from rasterio import features
from scipy.ndimage import distance_transform_edt

from .pytools import geotransform, load_geotiff


# https://gist.github.com/calebrob6/5039fd409921606aa5843f8fec038c03
def download_roads(domain):
    arr, proj, gt = load_geotiff(domain)
    proj_xy = pyproj.Proj(proj)
    proj_ll = pyproj.Proj("epsg:4326")
    xy2ll = pyproj.Transformer.from_proj(proj_xy, proj_ll, always_xy=True)

    def func(x, y):
        return xy2ll.transform(*geotransform(x, y, gt))

    H, W = arr.shape
    n = func(np.arange(W), 0)[1]
    s = func(np.arange(W), H - 1)[1]
    y = np.concatenate([n, s])
    e = func(W - 1, np.arange(H))[0]
    w = func(0, np.arange(H))[0]
    x = np.concatenate([e, w])

    coords = (y.max(), y.min(), x.max(), x.min())
    print("Downloading road network")
    Graph = ox.graph_from_bbox(
        *coords,
        network_type="drive",
        simplify=False,
        retain_all=True,
        truncate_by_edge=True,
        clean_periphery=True,
    )
    return ox.graph_to_gdfs(Graph, nodes=False)


# https://gis.stackexchange.com/a/151861/92556
def rasterize_roads(domain, roads):
    rst = rasterio.open(domain)
    meta = rst.meta.copy()
    meta["compress"] = "lzw"
    print("Rasterizing roads")
    roads = roads.to_crs(meta["crs"])
    burned = features.rasterize(
        roads.geometry,
        out_shape=rst.shape,
        fill=1,
        default_value=0,
        all_touched=True,
        transform=rst.transform,
    )
    return burned


def dist_ang(domain, roads):
    arr, proj, gt = load_geotiff(domain)
    print("Computing geometrical parameters")
    dist, i = distance_transform_edt(
        roads, return_indices=True, sampling=(abs(gt[5]), gt[1])
    )
    Y, X = np.indices(roads.shape, sparse=True)
    return dist, np.arctan2((i[1] - X) * gt[1], (i[0] - Y) * gt[5])


def roads_analysis(domain):
    return dist_ang(domain, rasterize_roads(domain, download_roads(domain)))
