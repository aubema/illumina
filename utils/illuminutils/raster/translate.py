# -*- coding: utf-8 -*-
""""Translate raster between two formats

:Author: Jean-Denis Giguère
:Licence: GPL v3 or later
:Date: 2013-02-19
"""

import logging
import subprocess
import sys

class Translate(object):
    """A raster dataset to translate and its transformation parameters"""

    def __init__(self):
        self._dataset = None
        self._dstfile = None
        self._nodata = None

    @property
    def dataset(self):
        "Dataset from which extract the raster"
        return self._dataset

    @dataset.setter
    def dataset(self,  fname):
        self._dataset = fname

    @property
    def dstfile(self):
        "Name of translated raster"
        return self._dstfile

    @dstfile.setter
    def dstfile(self,  fname):
        self._dstfile = fname

    @property
    def nodata(self):
        "Numeric value for nodata"
        return self._nodata

    @nodata.setter
    def nodata(self, nodataval):
        self._nodata = nodataval

    def check_translation_params(self):
        "Make sure that translation params are coherents"
        if self.dataset is None:
            logging.error("Translate has no dataset name")
            sys.exit(1)
        if self.dstfile is None:
            logging.error("Translate has no destination name")
            sys.exit(1)

    def translation_cmd(self):
        "Return a string reprensentation of the translation command"
        self.check_translation_params()
        translate_cmd = ['gdal_translate', ]
        #TODO: Check if no data management is required for translation
        #if self.nodata is not None:
        #    translate_cmd += ['-a_nodata',  str(self.nodata)]
        translate_cmd += [self.dataset,  self.dstfile]
        self._translate_command_list = translate_cmd
        return ' '.join(translate_cmd)

    def translate(self):
        "Translate the raster dataset"
        self.translation_cmd()
        subprocess.call(self._translate_command_list, shell=True)
