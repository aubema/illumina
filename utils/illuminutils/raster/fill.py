# -*- coding: utf-8 -*-
""""Fill nodata in raster

:Author: Jean-Denis Giguère
:Licence: GPL v3 or later
:Date: 2013-02-16
"""

import logging
import subprocess
import sys

class Fill(object):
    """A raster in which null will be filled and the transformation parameter"""

    def __init__(self):
        self._srcfile = None
        self._dstfile = None
        self._nodata = None

    @property
    def srcfile(self):
        "Raster to reproject"
        return self._srcfile

    @srcfile.setter
    def srcfile(self,  fname):
        self._srcfile = fname

    @property
    def dstfile(self):
        "Name of reprojected raster"
        return self._dstfile

    @dstfile.setter
    def dstfile(self,  fname):
        self._dstfile = fname

    @property
    def nodata(self):
        "Numeric value for nodata in tiles"
        return self._nodata

    @nodata.setter
    def nodata(self, nodataval):
        self._nodata = nodataval

    def check_transform_params(self):
        "Make sure that reprojection params are coherents"
        if self.srcfile is None:
            logging.error("Reproject has no source name")
            sys.exit(1)
        if self.dstfile is None:
            logging.error("Reproject has no destination name")
            sys.exit(1)

    def fill_cmd(self):
        "Return a string reprensentation of the fill command"
        self.check_transform_params()
        fill_cmd = ['gdal_fillnodata.py', ]
        fill_cmd += [self.srcfile,  self.dstfile]
        self._fill_command_list = fill_cmd
        return ' '.join(fill_cmd)

    def fill(self):
        "Fill the raster"
        self.fill_cmd()
        subprocess.call(self._fill_command_list)
