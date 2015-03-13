# -*- coding: utf-8 -*-
""""Reproject raster

:Author: Jean-Denis Giguère
:Licence: GPL v3 or later
:Date: 2013-02-16
"""

import logging
import subprocess
import sys

class Reproject(object):
    """A raster to reproject and the transformation parameter"""

    def __init__(self):
        self._srcfile = None
        self._dstfile = None
        self._nodata = None
        self._resampling_method = 'bilinear'

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
    def proj4string(self):
        "Proj4 string representing the spatial reference system"
        return self._proj4string

    @proj4string.setter
    def proj4string(self,  proj4str):
        self._proj4string = proj4str

    @property
    def nodata(self):
        "Numeric value for nodata in tiles"
        return self._nodata

    @nodata.setter
    def nodata(self, nodataval):
        self._nodata = nodataval

    @property
    def resampling_method(self):
        return self._resampling_method

    @resampling_method.setter
    def resampling_method(self,  method_name):
        self._resampling_method = method_name

    def check_transform_params(self):
        "Make sure that reprojection params are coherents"
        if self.srcfile is None:
            logging.error("Reproject has no source name")
            sys.exit(1)
        if self.dstfile is None:
            logging.error("Reproject has no destination name")
            sys.exit(1)

    def reproject_cmd(self):
        "Return a string reprensentation of the reproject command"
        self.check_transform_params()
        reproject_cmd = ['gdalwarp', ]
        reproject_cmd += ['-t_srs',  self.proj4string]
        if self.nodata is not None:
            reproject_cmd += ['-srcnodata',  str(self.nodata),  '-dstnodata',  str(self.nodata)]
        reproject_cmd += ['-r',  self.resampling_method]
        reproject_cmd += [self.srcfile,  self.dstfile]
        self._reproject_command_list = reproject_cmd
        return ' '.join(reproject_cmd)

    def reproject(self):
        "Reproject the raster"
        self.reproject_cmd()
        subprocess.call(self._reproject_command_list)
