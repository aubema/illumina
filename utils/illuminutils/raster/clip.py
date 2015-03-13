# -*- coding: utf-8 -*-
""""Clip raster and change pixels size

:Author: Jean-Denis Giguère
:Licence: GPL v3 or later
:Date: 2013-02-16
"""

import logging
import subprocess
import sys

class Clip(object):
    """A raster to clip and the transformation parameter"""
    def __init__(self):
        self._srcfile = None
        self._dstfile = None
        self._nodata = None
        self._resampling_method = 'bilinear'
        self._pixsize = None
        self._bbox = None

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

    @property
    def pixsize(self):
        "Size of pixels [meters]"
        return self._pixsize

    @pixsize.setter
    def pixsize(self,  size):
        self._pixsize = size

    @property
    def bbox(self):
        "Bounding box as array [minx, miny, maxx, maxy]"
        return self._bbox

    @bbox.setter
    def bbox(self,  bbox):
        self._bbox = bbox

    @property
    def resampling_method(self):
        return self._resampling_method

    @resampling_method.setter
    def resampling_method(self,  method_name):
        self._resampling_method = method_name

    def check_clip_params(self):
        "Make sure that reprojection params are coherents"
        if self.srcfile is None:
            logging.error("Clip has no source name")
            sys.exit(1)
        if self.dstfile is None:
            logging.error("Clip has no destination name")
            sys.exit(1)

    def clip_cmd(self):
        "Clip command as as string to execute in shell"
        self.check_clip_params()
        clip_cmd = ['gdalwarp', ]
        clip_cmd += ['-te', ] + [str(b) for b in self.bbox]
        clip_cmd += ['-tr', str(self.pixsize),  str(self.pixsize)]
        if self.nodata is not None:
            clip_cmd += ['-srcnodata',  str(self.nodata),  '-dstnodata',  str(self.nodata)]
        clip_cmd += ['-r',  self.resampling_method]
        clip_cmd += [self.srcfile,  self.dstfile]
        self._clip_command_list = clip_cmd
        return ' '.join(clip_cmd)

    def clip(self):
        "Clip the raster"
        self.clip_cmd()
        subprocess.call(self._clip_command_list)


