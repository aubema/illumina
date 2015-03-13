# -*- coding: utf-8 -*-
""""Build mosaic from tiles

:Author: Jean-Denis Giguère
:Licence: GPL v3 or later
:Date: 2013-02-10
"""

import logging
import subprocess
import sys

class Mosaic(object):
    def __init__(self):
        """Tiles to merge"""
        self._nodata = None
        self._tiles = None
        self._filename = None

    @property
    def tiles(self):
        "List of every tiles to be merged"
        return self._tiles

    @tiles.setter
    def tiles(self,  tiles_list):
        "Set tiles based on an array"
        self._tiles = tiles_list

    @property
    def nodata(self):
        "Numeric value for nodata in tiles"
        return self._nodata

    @nodata.setter
    def nodata(self, nodataval):
        self._nodata = nodataval

    @property
    def filename(self):
        "Filename of the resulting mosaic"
        return self._filename

    @filename.setter
    def filename(self,  filename):
        self._filename = filename

    def check_mosaic_params(self):
        "Make sure that mosaic params are coherents"
        if self.filename is None:
            logging.error("Mosaic has no name")
            sys.exit(1)
        try:
            if len(self.tiles) < 1:
                logging.error("At least on tile is required")
                sys.exit(1)
        except TypeError:
            logging.error("Invalid tiles list")
            sys.exit(1)

    def build_cmd(self):
        "Return a string reprensentation of the build command"
        self.check_mosaic_params()
        merge_command = ['gdal_merge.py', '-o', self.filename,  '-of', 'GTiff']
        if self.nodata is not None:
            merge_command += ['-n', str(self.nodata),  '-a_nodata', str(self.nodata)]
        merge_command += self.tiles
        self._merge_command_list = merge_command
        return ' '.join(merge_command)

    def build(self):
        "Build the mosaic"
        self.build_cmd()
        subprocess.call(self._merge_command_list)
