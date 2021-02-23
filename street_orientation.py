#!/usr/bin/env python3

# ******************************************************************************
#        Determine street orientation from list of coordinates
#
# Authors : Julien-Pierre Houle & Alexandre Simoneau
# Date :    February 2021
# ******************************************************************************

import numpy as np
import osmnx as ox

def get_bearing(lat1, lon1, lat2, lon2):
    x = np.cos(np.deg2rad(lat2)) * np.sin(np.deg2rad(lon2-lon1))
    y = np.cos(np.deg2rad(lat1)) * np.sin(np.deg2rad(lat2)) - \
        np.sin(np.deg2rad(lat1)) * np.cos(np.deg2rad(lat2)) * np.cos(np.deg2rad(lon2-lon1))
    return np.rad2deg(np.arctan2(x,y))

def street_orientation(lats, lons):
    print("    Loading graph")
    Graph = ox.graph_from_bbox(
        north=max(lats)+0.01,
        south=min(lats)-0.01,
        east=max(lons)+0.01,
        west=min(lons)-0.01,
        network_type='drive',
        simplify=False,
        retain_all=True,
        truncate_by_edge=True,
        clean_periphery=True
    )
    Graph = ox.utils_graph.get_undirected(Graph)
    Graph = ox.bearing.add_edge_bearings(Graph, precision=0)
    nodes, edges = ox.graph_to_gdfs(Graph)
    df_routes = edges.filter(['name','bearing','geometry'], axis=1)

    print("    Get nearest edges")
    edges_ID = ox.distance.get_nearest_edges(Graph,lons,lats,method="kdtree")
    nearest_edges = df_routes.loc[map(tuple,edges_ID)]

    print("    Compute bearings")
    coords = np.array([ x.coords.xy for x in nearest_edges['geometry'] ])

    bearing = nearest_edges['bearing'].to_numpy()
    bearing_AC = get_bearing(coords[:,1,0],coords[:,0,0],lats,lons)
    bearing[(bearing_AC - bearing) % 360 > 180] += 180
    bearing -= 90 # point towards road

    return bearing % 360
