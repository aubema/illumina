#!/usr/bin/env python3

import geopandas as gpd
import numpy as np
import osmnx as ox
import pandas as pd
import shapely.geometry as geometry
import yaml
from osgeo import gdal
from pyproj import Transformer

with open("iss_params.in") as f:
    p = yaml.safe_load(f)


def open_tiff(filename, dtype=np.float32):
    # Load file, and access the band and get a NumPy array
    src = gdal.Open(filename, gdal.GA_Update)
    band = src.GetRasterBand(1)
    ar = band.ReadAsArray()
    return src, ar


def geotransform(X, Y, GT):
    # From https://gdal.org/tutorials/geotransforms_tut.html
    X_geo = GT[0] + (X + 0.5) * GT[1] + (Y + 0.5) * GT[2]
    Y_geo = GT[3] + (X + 0.5) * GT[4] + (Y + 0.5) * GT[5]
    return X_geo, Y_geo


src, im = open_tiff(f"{p['wd']}/Vrad.tiff")
GT = src.GetGeoTransform()

src, tech = open_tiff(f"{p['wd']}/tech.tiff")
src, domain = open_tiff(f"{p['wd']}/domain.tiff")

df = pd.DataFrame()
Y, X = np.where((im > 0) & (tech > 0))
df["val"] = im[Y, X]
df["tech"] = tech[Y, X]
df["lons"], df["lats"] = geotransform(X, Y, GT)
df.to_pickle(f"{p['wd']}/xyz.pickle")

df = pd.DataFrame()
Y, X = np.where(domain)
df["lons"], df["lats"] = geotransform(X, Y, GT)
df["X"] = X
df["Y"] = Y

outProj = (
    "epsg:32"
    + ("6" if np.mean(df["lats"]) >= 0 else "7")
    + "%02d" % (np.mean(df["lons"]) / 6 + 31)
)  # WGS84/UTM

print("Build convex hull")
corner_mask = (
    (df["lons"] < (0.975 * np.min(df["lons"]) + 0.025 * np.max(df["lons"])))
    | (df["lons"] > (0.025 * np.min(df["lons"]) + 0.975 * np.max(df["lons"])))
    | (df["lats"] < (0.975 * np.min(df["lats"]) + 0.025 * np.max(df["lats"])))
    | (df["lats"] > (0.025 * np.min(df["lats"]) + 0.975 * np.max(df["lats"])))
)
corner_pts = df[corner_mask][["lons", "lats"]].to_numpy()
convex_hull = geometry.MultiPoint(corner_pts).convex_hull


print("Load graph")
Graph = ox.graph_from_polygon(
    convex_hull,
    network_type="drive",
    simplify=False,
    retain_all=True,
    truncate_by_edge=True,
    clean_periphery=True,
)

print("Process graph")
Graph = ox.utils_graph.get_undirected(Graph)
Graph = ox.bearing.add_edge_bearings(Graph, precision=0)
Graph = ox.projection.project_graph(Graph, to_crs=outProj)
nodes, edges = ox.graph_to_gdfs(Graph)
df_routes = edges.filter(["name", "bearing", "geometry"], axis=1)
del nodes
del edges

inProj = "epsg:4326"
X, Y = Transformer.from_crs(inProj, outProj, always_xy=True).transform(
    df["lons"], df["lats"]
)

print("Get nearest edges")
nearest_edges = df_routes.loc[
    map(tuple, ox.distance.nearest_edges(Graph, X, Y))
]
del Graph
del Y
del X

print("Compute distance")
df["distance"] = nearest_edges.distance(
    gpd.GeoSeries(
        gpd.points_from_xy(df["lons"], df["lats"], crs=inProj)
    ).to_crs(outProj),
    align=False,
).to_numpy()

df.to_csv(
    f"{p['wd']}/obs.csv",
    columns=["lons", "lats", "X", "Y", "distance"],
    header=True,
    index=False,
)
