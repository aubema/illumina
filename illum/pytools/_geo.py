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


# https://gdal.org/tutorials/geotransforms_tut.html
def geotransform(X, Y, GT):
    X_geo = GT[0] + (X + 0.5) * GT[1] + (Y + 0.5) * GT[2]
    Y_geo = GT[3] + (X + 0.5) * GT[4] + (Y + 0.5) * GT[5]
    return X_geo, Y_geo


# https://gist.github.com/calebrob6/5039fd409921606aa5843f8fec038c03
def download_roads(rst):
    proj_xy = pyproj.Proj(rst.crs)
    proj_ll = pyproj.Proj("epsg:4326")
    xy2ll = pyproj.Transformer.from_proj(proj_xy, proj_ll, always_xy=True)

    def func(x, y):
        return xy2ll.transform(*geotransform(x, y, rst.get_transform()))

    n = func(np.arange(rst.width), 0)[1]
    s = func(np.arange(rst.width), rst.height - 1)[1]
    y = np.concatenate([n, s])
    e = func(rst.width - 1, np.arange(rst.height))[0]
    w = func(0, np.arange(rst.height))[0]
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
def rasterize_roads(rst, roads):
    meta = rst.meta.copy()
    meta["compress"] = "lzw"
    print("Rasterizing roads")
    roads = roads.to_crs(rst.crs)
    burned = features.rasterize(
        roads.geometry,
        out_shape=rst.shape,
        fill=1,
        default_value=0,
        all_touched=True,
        transform=rst.transform,
    )
    return burned


def dist_ang(gt, roads):
    print("Computing geometrical parameters")
    dist, i = distance_transform_edt(
        roads, return_indices=True, sampling=(abs(gt[5]), gt[1])
    )
    Y, X = np.indices(roads.shape, sparse=True)
    return dist, np.arctan2((i[1] - X) * gt[1], (i[0] - Y) * gt[5])


def roads_analysis(domain):
    rst = rasterio.open(domain)
    network = download_roads(rst)
    burned = rasterize_roads(rst, network)
    return dist_ang(rst.get_transform(), burned)


def compute_ground_type(dist, bounds):
    ground_type = np.ones(dist.shape) * len(bounds)
    for i, lim in reversed(list(enumerate(sorted(bounds)))):
        ground_type[dist < lim] = i
