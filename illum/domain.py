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
        [" ", " k", " M", " G"]
        + ([None] * 10)
        + [" f", " p", " n", " \u03bc", " m"]
    )

    power_of_1000 = int(math.floor(math.log10(x) // 3))
    exponent = 3 * power_of_1000
    prefix = UNITS[power_of_1000]
    if prefix is None:
        prefix = "*10^%d " % exponent

    significand = x * 10 ** (-exponent)
    return "%.2f%s%s" % (significand, prefix, unit)


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

    wgs84 = pyproj.Proj("epsg:4326")
    proj = pyproj.Proj(domain["srs"])

    x0, y0 = pyproj.transform(
        wgs84, proj, center_lon, center_lat, always_xy=True
    )

    obs_x, obs_y = list(
        zip(
            *(
                pyproj.transform(wgs84, proj, lon, lat, always_xy=True)
                for lat, lon in zip(obs_lat, obs_lon)
            )
        )
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
        psize = domain["scale_min"] * scale ** i
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
    SE = pyproj.transform(proj, wgs84, xmax, ymin, always_xy=True)
    SW = pyproj.transform(proj, wgs84, xmin, ymin, always_xy=True)
    NE = pyproj.transform(proj, wgs84, xmax, ymax, always_xy=True)
    NW = pyproj.transform(proj, wgs84, xmin, ymax, always_xy=True)

    N = max(NE[1], NW[1])
    S = min(SE[1], SW[1])
    E = max(NE[0], SE[0])
    W = min(NW[0], SW[0])

    print("Bounding box:")
    print("  SW: %f,%f" % (S, W))
    print("  NE: %f,%f" % (N, E))
