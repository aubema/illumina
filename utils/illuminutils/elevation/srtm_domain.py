# -*- coding: utf-8 -*-
"""Helper functions for SRTM data download

Input:
    Illumina Parameters ini file (hardcoded)

:Author: Jean-Denis Giguère
:Licence: GPL v3 or later
:Date: 2013-02-10
"""
def latitude_tiles(prefix, t_range):
    tiles = []
    for t in range(t_range[0], t_range[1] + 1):
        tile_part = "%s%02d" % (prefix, t)
        tiles.append(tile_part)

    return tiles

def longitude_tiles(prefix, t_range):
    tiles = []
    for t in range(t_range[0], t_range[1] + 1):
        tile_part = "%s%03d" % (prefix, t)
        tiles.append(tile_part)

    return tiles


def all_tiles(lat_tiles, lon_tiles):
    tiles = []
    for lat in lat_tiles:
        for lon in lon_tiles:
            tiles.append("%s%s" % (lat, lon))

    return tiles

latitude_tiles_prefix = 'N'
latitude_tiles_range = (33, 34)

longitude_tiles_prefix = 'W'
longitude_tiles_range = (116,117)

tiles = all_tiles(latitude_tiles(latitude_tiles_prefix,latitude_tiles_range),
    longitude_tiles(longitude_tiles_prefix, longitude_tiles_range))
