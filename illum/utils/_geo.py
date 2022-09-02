#!/usr/bin/env python3
# ******************************************************************************
#        Determine street orientation from list of coordinates
#
# Authors : Julien-Pierre Houle & Alexandre Simoneau
# Date :    May 2021
# ******************************************************************************

from math import ceil

import illum.utils as u
import numpy as np
import osmnx as ox
import pyproj
import rasterio
import rasterio.features
import scipy.ndimage
import shapely.geometry
import shapely.ops
from progressbar import progressbar


def transform(s_crs=4326, t_crs=4326):
    s_crs = pyproj.CRS.from_epsg(s_crs)
    t_crs = pyproj.CRS.from_epsg(t_crs)
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
            np.min(lon), np.min(lat), np.max(lon), np.max(lat)
        ),
    )[0].code


def fishnet(bounds, res, xsize=None, ysize=None, cols=None, rows=None):
    if xsize is not None and cols is not None:
        print("WARNING: 'xsize' and 'cols' provided. Ignoring 'cols'.")
    if ysize is not None and rows is not None:
        print("WARNING: 'ysize' and 'rows' provided. Ignoring 'rows'.")

    xmin, ymin, xmax, ymax = bounds.bounds
    if cols is None:
        cols = round((xmax - xmin) / xsize)
    xsize = ceil((xmax - xmin) / cols / res) * res
    if rows is None:
        rows = round((ymax - ymin) / ysize)
    ysize = ceil((ymax - ymin) / rows / res) * res

    geoms = {}
    for row in range(rows):
        for col in range(cols):
            box = shapely.geometry.box(
                xmin + xsize * col,
                ymax - ysize * (row + 1),
                xmin + xsize * (col + 1),
                ymax - ysize * row,
            )
            geo = bounds.intersection(box)
            if not geo.is_empty:
                geoms[row, col] = [box, geo]

    return geoms, (int(ysize // res), int(xsize // res))


# https://gist.github.com/calebrob6/5039fd409921606aa5843f8fec038c03
def download_roads(bounds, crs):
    bounds_lonlat = shapely.ops.transform(transform(s_crs=crs), bounds)
    Graph = ox.graph_from_polygon(
        bounds_lonlat,
        network_type="drive",
        simplify=False,
        retain_all=True,
        truncate_by_edge=True,
    )
    return ox.graph_to_gdfs(Graph, nodes=False)


# https://gis.stackexchange.com/a/151861/92556
def rasterize_roads(roads, bounds, epsg, xres, yres):
    crs = pyproj.CRS.from_epsg(epsg)
    roads = roads.to_crs(crs)
    xmin, ymin, xmax, ymax = bounds.bounds
    width = int((xmax - xmin) / xres)
    height = int((ymax - ymin) / yres)
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


def roads_analysis(domain, res=1, xsize=10000, ysize=10000, buffer=200):
    print("Analyzing domain.")
    rst = rasterio.open(domain)
    coords = rasterio.transform.xy(rst.transform, *np.where(rst.read(1)))
    epsg = estimate_utm_epsg(*coords)
    print(f"Using epsg:{epsg} projection.")
    bounds = shapely.geometry.MultiPoint(
        tuple(zip(*transform(t_crs=epsg)(*coords)))
    ).minimum_rotated_rectangle

    print("Splitting domain.")
    res = 1
    buffer = ceil(200 / res) * res
    tiles, shape = fishnet(bounds, res=res, xsize=10000, ysize=10000)
    ny, nx = shape
    nb = int(buffer // res)
    mask = slice(nb, -nb), slice(nb, -nb)
    Ty, Tx = np.max(list(tiles.keys()), 0) + 1
    print(f"Domain split in {len(tiles)} ({Ty}x{Tx})")

    for idx, tile in progressbar(tiles.items(), redirect_stdout=True):
        y, x = idx
        box, geo = tile
        boxb = box.buffer(buffer, join_style=2)
        geob = geo.buffer(buffer, join_style=2)

        print(f"Processing tile {idx}.")
        # print("Downloading road network.")
        network = download_roads(geob, epsg)
        # print("Rasterizing road network.")
        burned = rasterize_roads(network, boxb, epsg, res, res)

        # print("Computing distance and angle to roads.")
        dist, ang = dist_ang(burned, res, res)

        profile = dict(
            crs=pyproj.CRS.from_epsg(epsg),
            transform=rasterio.transform.from_bounds(*box.bounds, nx, ny),
        )
        u.save_geotiff(f"road_distance_x{x}_y{y}.tif", dist[mask], **profile)
        u.save_geotiff(f"angle_to_road_x{x}_y{y}.tif", ang[mask], **profile)


def compute_ground_type(dist, bounds):
    ground_type = np.ones(dist.shape) * len(bounds)
    for i, lim in reversed(list(enumerate(sorted(bounds)))):
        ground_type[dist < lim] = i
