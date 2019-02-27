#!/usr/bin/env python
# -*- coding: utf-8 -*-
""""SRTM for elevation data

Input:
    Illumina Parameters ini file (hardcoded)
    Hardcoded run status

Output:
    Final elevation PGM and intermediate files

:Author: Jean-Denis Giguère
:Licence: GPL v3 or later
:Date: 2013-02-10
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

import srtm_domain


class SRTM(object):
    def __init__(self):
        self.nodata = -32768
        #TODO: Ini file is hardcoded
        self.params = IlluminaParameters('palomar.ini')
	

        #TODO: Domain is static in strm_domain.
        #      Build it based on bounding box
        self._tiles = [os.path.join("srtm",  "%s.hgt" % (tname, )) for tname in srtm_domain.tiles]
        self._mosaic_name = os.path.join("srtm","srtm.tif")
        self._reprojected_name = os.path.join("srtm",  "srtm_%s.tif" % (self.params.srs_suffixname))
        self._clipped_name = os.path.join("srtm", "srtm_%s_clipped.tif" % (self.params.srs_suffixname) )
        self._filled_name = os.path.join("srtm", "srtm_%s_filled.tif" % (self.params.srs_suffixname) )
        self._pgm_name = os.path.join("srtm","srtm.pgm")
        self._mosaic_built = False
        self._mosaic_reprojected = False
        self._mosaic_clipped = False
        self._mosaic_filled = False
        self._mosaic_saved_as_pgm = False

    @property
    def tiles(self):
        """Tiles list"""
        return self._tiles

    @property
    def mosaic_name(self):
        "Mosaic output filaname"
        logging.debug("Getting SRTM mosaic_name: %s" % (self._mosaic_name))
        return self._mosaic_name

    @property
    def mosaic_built(self):
        "If mosaic is built, we will not build it again"
        return self._mosaic_built

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
    def mosaic_reprojected(self):
        "If mosaic is reprojected, we will not reproject it again"
        return self._mosaic_reprojected

    @property
    def mosaic_clipped(self):
        "If mosaic is clipped, we will not clip it again"
        return self._mosaic_clipped

    @property
    def mosaic_filled(self):
        "If mosaic is filled, we will not fill it again"
        return self._mosaic_filled

    @property
    def mosaic_saved_as_pgm(self):
        "If mosaic as already been saved as pgm, we will not save it again"
        return self._mosaic_saved_as_pgm

    def download(self):
        "TODO: Download required SRTM tiles"
        pass

    def mosaic(self):
        logging.debug("SRTM mosaic configuration...")
        srtm_mosaic = Mosaic()
        srtm_mosaic.nodata = self.nodata
        srtm_mosaic.tiles = self.tiles
        logging.debug("Setting SRTM mosaic name to %s" % (self.mosaic_name, ))
        srtm_mosaic.filename = self.mosaic_name
        logging.info("Building using SRTM mosaic")
        #logging.log(VERBOSE_LVL,  srtm_mosaic.build_cmd)
        srtm_mosaic.build()

    def reproject(self):
        srtm_reproject = Reproject()
        srtm_reproject.nodata = self.nodata
        srtm_reproject.srcfile = self.mosaic_name
        srtm_reproject.dstfile = self.reprojected_name
        srtm_reproject.proj4string = self.params.proj4string
        logging.info("Reprojectiong SRTM mosaic")
        srtm_reproject.reproject()

    def clip(self):
        srtm_clip = Clip()
        srtm_clip.nodata = self.nodata
        srtm_clip.srcfile = self.reprojected_name
        srtm_clip.dstfile = self.clipped_name
        srtm_clip.pixsize = self.params.pixsize
        srtm_clip.bbox = self.params.bbox
        srtm_clip.clip()

    def fill(self):
        srtm_fill = Fill()
        srtm_fill.nodata = self.nodata
        srtm_fill.srcfile = self.clipped_name
        srtm_fill.dstfile = self.filled_name
        srtm_fill.fill()

    def save_pgm(self):
        logging.info("Saving SRTM as PGM")
        srtm_pgm = Pgm()
        srtm_pgm.srcfile = self.filled_name
        srtm_pgm.dstfile = self.pgm_name
        srtm_pgm.write()

if __name__ == '__main__':
    start_logging()
    srtm_dataset = SRTM()
    if not srtm_dataset.mosaic_built:
        srtm_dataset.mosaic()
    if not srtm_dataset.mosaic_reprojected:
        srtm_dataset.reproject()
    if not srtm_dataset.mosaic_clipped:
        srtm_dataset.clip()
    if not srtm_dataset.mosaic_filled:
        srtm_dataset.fill()
    if not srtm_dataset.mosaic_saved_as_pgm:
        srtm_dataset.save_pgm()
    logging.info("SRTM elevation dataset is ready!")
