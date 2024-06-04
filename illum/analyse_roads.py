#!/usr/bin/env python3

import multiprocessing
import os
import shutil

import geopandas as gpd
import numpy as np
import osmnx
import pandas as pd
import rasterio as rio
import rasterio.features
import scipy
import shapely
import yaml
from tqdm import tqdm


def pmap(func, it):
    with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
        return list(tqdm(pool.imap_unordered(func, it), total=len(it)))


def fishnet(bounds, size):
    xmin, ymin, xmax, ymax = bounds.bounds
    cols = round((xmax - xmin) / size)
    rows = round((ymax - ymin) / size)

    geoms = {}
    for row in range(rows):
        for col in range(cols):
            box = shapely.geometry.box(
                xmin + size * col,
                ymax - size * (row + 1),
                xmin + size * (col + 1),
                ymax - size * row,
            )
            geo = bounds.intersection(box)
            if not geo.is_empty:
                geoms[row, col] = (box, geo)

    return geoms


def download_network(limits, crs):
    bounds_lonlat = gpd.GeoSeries([limits], crs=crs).to_crs(4326)[0]
    # print("Downloading road network.")
    graph = osmnx.graph_from_polygon(
        bounds_lonlat,
        network_type="drive",
        simplify=False,
        retain_all=True,
        truncate_by_edge=True,
    )
    # print("Processing.")
    return osmnx.projection.project_graph(graph, to_crs=crs)


def get_bounds(bounds, resolution):
    xmin, ymin, xmax, ymax = bounds.bounds
    xmin = xmin // res * res
    ymin = ymin // res * res
    xmax = xmax // res * res
    ymax = ymax // res * res

    height = int((ymax - ymin) // res)
    width = int((xmax - xmin) // res)

    return xmin, ymin, xmax, ymax, width, height


def analyse_roads(bounds, graph, crs, X, Y):
    roads = osmnx.graph_to_gdfs(graph, nodes=False).filter(["geometry"])

    # print("Finding closest road.")
    nearest_roads = roads.loc[
        osmnx.distance.nearest_edges(graph, X.flat, Y.flat)
    ].reset_index(drop=True)
    del roads

    # print("Computing bearing.")
    pts = gpd.GeoSeries(gpd.points_from_xy(X.flat, Y.flat, crs=crs))

    distance = pts.distance(nearest_roads).to_numpy().reshape(X.shape)
    short = shapely.get_coordinates(
        pts.shortest_line(nearest_roads).to_crs(4326)
    ).reshape(-1, 4)
    del pts
    del nearest_roads

    bearing = osmnx.bearing.calculate_bearing(
        short[:, 1], short[:, 0], short[:, 3], short[:, 2]
    ).reshape(X.shape)

    return distance, bearing


if __name__ == "__main__":
    with open("iss_params.in") as f:
        p = yaml.safe_load(f)

    src = rio.open("domain.tiff")
    arr = src.read(1)
    arr[5:-5, 5:-5] = 0

    center_lon = (src.bounds.left + src.bounds.right) / 2
    center_lat = (src.bounds.top + src.bounds.bottom) / 2
    crs = "32" + ("6" if center_lat >= 0 else "7") + "%02d" % (center_lon / 6 + 31)
    domain = shapely.geometry.MultiPoint(
        gpd.points_from_xy(
            *rio.transform.xy(src.transform, *np.where(arr)), crs=4326
        ).to_crs(crs)
    ).minimum_rotated_rectangle
    buff_domain = domain.buffer(p["buffer"] * p["resolution"], join_style="mitre")

    grid = fishnet(domain, p["size"])

    def process(elem):
        box, geo = elem

        box_buff = box.buffer(p["buffer"] * p["resolution"], join_style="mitre")
        geo_buff = geo.buffer(
            p["buffer"] * p["resolution"], join_style="mitre"
        ).intersection(buff_domain)

        x0, y0, x1, y1, w, h = get_bounds(box_buff, p["resolution"])
        xs = np.linspace(x0 + resolution / 2, x1, w)
        ys = np.linspace(y1 - resolution / 2, y0, h)
        X, Y = np.meshgrid(xs, ys)

        try:
            graph = download_network(geo_buff, crs)
        except ValueError:
            distance = angle = np.zeros_like(X)
        else:
            distance, angle = analyse_roads(box_buff, graph, crs, X, Y)

        return geo, box_buff, distance, angle

    results = pmap(process, grid.values())

    xmin, ymin, xmax, ymax, width, height = get_bounds(
        domain.buffer(buffer), p["resolution"]
    )
    transform = rio.transform.from_bounds(xmin, ymin, xmax, ymax, width, height)

    distance = np.full((height, width), np.nan)
    bearing = np.full((height, width), np.nan)

    for geo, box, dist, brng in results:
        Mask = rio.features.geometry_mask(
            [geo], (height, width), transform, invert=True
        )
        x0, y0, x1, y1, w, h = get_bounds(box, p["resolution"])
        trans = rio.transform.from_bounds(x0, y0, x1, y1, w, h)
        mask = rio.features.geometry_mask([geo], (h, w), trans, invert=True)
        distance[Mask] = dist[mask]
        bearing[Mask] = brng[mask]

    kwargs = dict(
        mode="w",
        driver="GTiff",
        height=height,
        width=width,
        count=1,
        dtype=distance.dtype,
        crs=rio.crs.CRS.from_epsg(crs),
        transform=transform,
    )

    with rio.open("distance.tiff", **kwargs) as dfile:
        dfile.write(distance, 1)

    with rio.open("bearing.tiff", **kwargs) as bfile:
        bfile.write(bearing, 1)
