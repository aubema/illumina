#!/usr/bin/env python3

import math

import click
import numpy as np
import pyproj
import yaml


def eng_format(x, unit=""):
    # Credit: 200_success on StackOverflow
    # https://codereview.stackexchange.com/a/50971
    #
    # U+03BC is Greek lowercase mu
    UNITS = (
        [" ", " k", " M", " G"] + ([None] * 10) + [" f", " p", " n", " \u03bc", " m"]
    )

    power_of_1000 = int(math.floor(math.log10(x) // 3))
    exponent = 3 * power_of_1000
    prefix = UNITS[power_of_1000]
    if prefix is None:
        prefix = "*10^%d " % exponent

    significand = x * 10 ** (-exponent)
    return f"{significand:.2f}{prefix}{unit}"


def round_odd(n):
    return int(n - n % 2 + 1)


@click.command(name="domain")
def CLI_domain():
    """Defines the simulation domain.

    Reads domain definition parameters from 'domain_params.in'.

    Outputs the definition in 'domain.ini'.
    """
    domain()


def domain():
    with open("domain_params.in") as f:
        domain = yaml.safe_load(f)

    obs_lat = domain.pop("latitude")
    obs_lon = domain.pop("longitude")

    try:
        if len(obs_lat) != len(obs_lon):
            print("ERROR: Latitude and Longitude must have the same length.")
            exit()
    except TypeError:  # lat and lon not lists
        obs_lat = [obs_lat]
        obs_lon = [obs_lon]

    center_lat = (max(obs_lat) + min(obs_lat)) / 2.0
    center_lon = (max(obs_lon) + min(obs_lon)) / 2.0

    # Define projection
    if domain["srs"] == "auto":
        default_srs = (
            "epsg:32"
            + ("6" if center_lat >= 0 else "7")
            + "%02d" % (center_lon / 6 + 31)
        )  # WGS84/UTM
        domain["srs"] = default_srs

    wgs84 = pyproj.CRS.from_epsg(4326)
    proj = pyproj.CRS.from_user_input(domain["srs"])
    transform = pyproj.Transformer.from_crs(wgs84, proj, always_xy=True).transform

    x0, y0 = transform(center_lon, center_lat)

    obs_x, obs_y = list(
        zip(*(transform(lon, lat) for lat, lon in zip(obs_lat, obs_lon)))
    )

    domain["observers"] = [
        {"latitude": lat, "longitude": lon, "x": x, "y": y}
        for lat, lon, x, y in zip(obs_lat, obs_lon, obs_x, obs_y)
    ]

    obs_x = np.array(obs_x)
    obs_y = np.array(obs_y)

    obs_size_x = 2 * np.max(np.abs(obs_x - x0))
    obs_size_y = 2 * np.max(np.abs(obs_y - y0))

    R = int(domain["nb_pixels"] / 2)
    r = int((domain["nb_pixels"] / float(domain.pop("scale_factor"))) / 2)
    scale = (R + 0.5) / (r + 0.5)

    domain["nb_pixels"] = R
    domain["nb_core"] = r
    domain["extents"] = list()

    for i in range(domain["nb_layers"]):
        psize = domain["scale_min"] * scale**i
        buff = min(255 - R, domain["buffer"] * 1e3 // psize)

        print("Layer", i)
        print("Pixel size:", eng_format(psize, "m"))
        print("Domain size:", eng_format(psize * (2 * R + 1), "m"))
        print("")

        n_obs_x = max(1.0, math.ceil(obs_size_x / psize))
        n_obs_y = max(1.0, math.ceil(obs_size_y / psize))

        layer_half_size_x = psize * (buff + R + n_obs_x / 2.0)
        layer_half_size_y = psize * (buff + R + n_obs_y / 2.0)

        xmin = x0 - layer_half_size_x
        xmax = x0 + layer_half_size_x
        ymin = y0 - layer_half_size_y
        ymax = y0 + layer_half_size_y

        extent = dict()
        extent["layer"] = i
        extent["pixel_size"] = psize
        extent["buffer"] = int(buff)
        extent["observer_size_x"] = int(n_obs_x)
        extent["observer_size_y"] = int(n_obs_y)
        bbox = dict(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
        extent.update(bbox)
        domain["extents"].append(extent)

    domain.pop("buffer")

    with open("domain.ini", "w") as f:
        yaml.safe_dump(domain, f, default_flow_style=False)

    # print lon/lat bbox formatted for earthdata
    transform = pyproj.Transformer.from_crs(proj, wgs84, always_xy=True).transform
    SE = transform(xmax, ymin)
    SW = transform(xmin, ymin)
    NE = transform(xmax, ymax)
    NW = transform(xmin, ymax)

    N = max(NE[1], NW[1])
    S = min(SE[1], SW[1])
    E = max(NE[0], SE[0])
    W = min(NW[0], SW[0])

    print("Bounding box:")
    print(f"  SW: {S:f},{W:f}")
    print(f"  NE: {N:f},{E:f}")
