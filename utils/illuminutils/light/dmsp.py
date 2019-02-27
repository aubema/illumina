#!/usr/bin/env python2
# -*- coding: utf-8 -*-
""""DMPS-OLS stable light

:Author: Jean-Denis Giguère
:Licence: GPL v3 or later
:Date: 2013-02-20
"""

import logging
import os.path

from illuminutils import VERBOSE_LVL, start_logging
from illuminutils.raster.mosaic import Mosaic
from illuminutils.raster.reproject import Reproject
from illuminutils.raster.clip import Clip
from illuminutils.raster.fill import Fill
from illuminutils.raster.pgm import Pgm
from illuminutils.experiments.illuminini import IlluminaParameters

class Dmsp(object):
    def __init__(self):
        #TODO: Ini file is hardcoded
        self.params = IlluminaParameters('madrid.ini')
        
        #TODO: Domain is static in strm_domain.
        #      Build it based on bounding box
        
        self._filename = os.path.join("DMSP-OLS", "F182010.v4c_web.stable_lights_coarse_clip.tif")
        self._reprojected_name = os.path.join("DMSP-OLS",  "stable_lights_%s.tif" % (self.params.srs_suffixname))
        self._clipped_name = os.path.join("DMSP-OLS", "stable_lights_%s_clipped.tif" % (self.params.srs_suffixname) )
        self._filled_name = os.path.join("DMSP-OLS", "stable_lights_%s_filled.tif" % (self.params.srs_suffixname) )
        self._pgm_name = os.path.join("DMSP-OLS","stable_lights.pgm")
        self._reprojected = True
        self._clipped = True
        self._filled = True
        self._saved_as_pgm = False
        
    @property
    def filename(self):
        """DMSP-OLS stable lights filename"""
        return self._filename
        
    @property
    def reprojected_name(self):
        "Name of the reprojected file"
        return self._reprojected_name
        
    @property
    def clipped_name(self):
        "Name of the clipped file"
        return self._clipped_name
        
    @property
    def filled_name(self):
        "Name of the filled file"
        return self._filled_name
        
    @property
    def  pgm_name(self):
        "Name of the pgm"
        return self._pgm_name
        
    @property
    def reprojected(self):
        "If dataset is reprojected, we will not reproject it again"
        return self._reprojected
        
    @property
    def clipped(self):
        "If dataset is clipped, we will not clip it again"
        return self._clipped
        
    @property
    def filled(self):
        "If dataset is filled, we will not fill it again"
        return self._filled
        
    @property
    def saved_as_pgm(self):
        "If dataset as already been saved as pgm, we will not save it again"
        return self._saved_as_pgm
        
        
    def reproject(self):
        reproject = Reproject()
        reproject.srcfile = self.filename
        reproject.dstfile = self.reprojected_name
        reproject.proj4string = self.params.proj4string
        logging.info("Reprojectiong DMSP-OLS mosaic")
        reproject.reproject()
        
    def clip(self):
        clip = Clip()
        clip.srcfile = self.reprojected_name
        clip.dstfile = self.clipped_name
        clip.pixsize = self.params.pixsize
        clip.bbox = self.params.bbox
        clip.clip()
        
    def fill(self):
        fill = Fill()
        fill.srcfile = self.clipped_name
        fill.dstfile = self.filled_name
        fill.fill()
        
    def save_pgm(self):
        logging.info("Saving DMSP-OLS stable light as PGM")
        pgm = Pgm()
        pgm.srcfile = self.clipped_name
        pgm.dstfile = self.pgm_name
        pgm.write()
        
if __name__ == '__main__':
    start_logging()
    dmsp_dataset = Dmsp()
    if not dmsp_dataset.reprojected:
        dmsp_dataset.reproject()
    if not dmsp_dataset.clipped:
        dmsp_dataset.clip()
    #if not dmsp_dataset.filled:
    #    dmsp_dataset.fill()
    if not dmsp_dataset.saved_as_pgm:
        dmsp_dataset.save_pgm()
    logging.info("DMSP-OLS stable lights dataset is ready!")
