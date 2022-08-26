#!/usr/bin/env python3
# ******************************************************************************
#        Determine street orientation from list of coordinates
#
# Authors : Julien-Pierre Houle & Alexandre Simoneau
# Date :    May 2021
# ******************************************************************************

from math import ceil

import numpy as np
import osmnx as ox
import pyproj
import rasterio
import rasterio.features
import scipy.ndimage
import shapely.geometry
import shapely.ops


def transform(s_crs=4326, t_crs=4326):
    return pyproj.Transformer.from_crs(s_crs, t_crs, always_xy=True).transform


# https://gdal.org/tutorials/geotransforms_tut.html
def geotransform(X, Y, GT):
    print("WARNING: Depreciated GDAL geotransform")
    X_geo = GT[0] + (X + 0.5) * GT[1] + (Y + 0.5) * GT[2]
    Y_geo = GT[3] + (X + 0.5) * GT[4] + (Y + 0.5) * GT[5]
    return X_geo, Y_geo


def estimate_utm_epsg(lon, lat):
    return pyproj.database.query_utm_crs_info(
        datum_name="WGS84",
        area_of_interest=pyproj.aoi.AreaOfInterest(
            min(lon), min(lat), max(lon), max(lat)
        ),
    )[0].code


def fishnet(bounds, xres=None, yres=None, cols=None, rows=None, buffer=200):
    if xres is not None and cols is not None:
        print("WARNING: 'xres' and 'cols' provided. Ignoring 'cols'.")
    if yres is not None and rows is not None:
        print("WARNING: 'yres' and 'rows' provided. Ignoring 'rows'.")

    xmin, ymin, xmax, ymax = bounds.bounds
    if xres is None:
        xres = (xmax - xmin) / cols
    cols = ceil((xmax - xmin) / xres)
    if yres is None:
        yres = (ymax - ymin) / rows
    rows = ceil((ymax - ymin) / yres)

    geoms = {}
    for row in range(rows):
        for col in range(cols):
            box = shapely.geometry.box(
                xmin + xres * col,
                ymax - yres * (row + 1),
                xmin + xres * (col + 1),
                ymax - yres * row,
            )
            geo = bounds.intersection(box).buffer(buffer, join_style=3)
            bbox = box.buffer(buffer, join_style=3)
            if not geo.is_empty:
                geoms[col, row] = box, geo, bbox

    return geoms


# https://gist.github.com/calebrob6/5039fd409921606aa5843f8fec038c03
def download_roads(bounds, crs):
    bounds_lonlat = shapely.ops.transform(transform(s_crs=crs), bounds)
    Graph = ox.graph_from_bbox(
        *bounds_lonlat.bounds,
        network_type="drive",
        simplify=False,
        retain_all=True,
        truncate_by_edge=True,
    )
    return ox.graph_to_gdfs(Graph, nodes=False)


# https://gis.stackexchange.com/a/151861/92556
def rasterize_roads(roads, bounds, crs, xres, yres):
    roads = roads.to_crs(crs)
    xmin, ymin, xmax, ymax = bounds.bouds
    width = (xmax - xmin) / xres
    height = (ymax - ymin) / yres
    transform = rasterio.transform.from_bounds(
        xmin, ymin, xmax, ymax, width, height
    )
    burned = rasterio.features.rasterize(
        roads.geometry,
        out_shape=(height, width),
        fill=1,
        default_value=0,
        all_touched=True,
        transform=transform,
    )
    return burned


def dist_ang(roads, xres, yres):
    dist, i = scipy.ndimage.distance_transform_edt(
        roads, return_indices=True, sampling=(abs(yres), xres)
    )
    Y, X = np.indices(roads.shape, sparse=True)
    return dist, np.arctan2((i[1] - X) * xres, (i[0] - Y) * yres)


def roads_analysis(domain):
    rst = rasterio.open(domain)
    coords = rasterio.transform.xy(rst.transform, *np.where(rst.read(1)))
    epsg = estimate_utm_epsg(*coords)
    bounds = shapely.geometry.MultiPoint(
        tuple(zip(*transform(t_srs=epsg)(*coords)))
    ).minimum_rotated_rectangle

    fishnet(bounds, xres=10000, yres=10000, buffer=0.1)

    xres, yres = 1, 1
    for poly in bounds.values():
        network = download_roads(poly, epsg)
        burned = rasterize_roads(network, poly, epsg, xres, yres)
        dist_ang(burned, xres, yres)


def compute_ground_type(dist, bounds):
    ground_type = np.ones(dist.shape) * len(bounds)
    for i, lim in reversed(list(enumerate(sorted(bounds)))):
        ground_type[dist < lim] = i
