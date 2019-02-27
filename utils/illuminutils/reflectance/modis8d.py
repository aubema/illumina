#!/usr/bin/env python2
# -*- coding: utf-8 -*-
""""Modis MYD09A1 reflectance data

:Author: Jean-Denis Giguère
:Licence: GPL v3 or later
:Date: 2013-02-10
"""

import logging
import os.path

from illuminutils import VERBOSE_LVL, start_logging
from illuminutils.raster.translate import Translate
from illuminutils.raster.reproject import Reproject
from illuminutils.raster.mosaic import Mosaic
from illuminutils.raster.clip import Clip
from illuminutils.raster.fill import Fill
from illuminutils.raster.pgm import Pgm
from illuminutils.experiments.illuminini import IlluminaParameters

class Band(object):
    def __init__(self):
        self.nodata = -28672
        #TODO: SRS suffix is encoded for modis data
        self.srs_suffix = 'utm30'
        #TODO: Ini file is hardcoded
        self.params = IlluminaParameters('palomar.ini')
        self._dataset_filenames = None
        self._suffix = None
        self._prefix = None
        self._translated = False
        self._mosaiced = False
        self._reprojected = False
        self._clipped = False
        self._filled = False
        self._saved_as_pgm = False

    @property
    def dataset_filenames(self):
        "HDF dataset containing this band"
        return self._dataset_filenames

    @dataset_filenames.setter
    def dataset_filenames(self, filenames):
        self._dataset_filenames = filenames

    @property
    def suffix(self):
        "Suffix of this specific band"
        return self._suffix

    @suffix.setter
    def suffix(self,  suffix):
        self._suffix = suffix

    @property
    def prefix(self):
        "Prefix of this dataset"
        return self._prefix

    @prefix.setter
    def prefix(self,  prefix):
        self._prefix = prefix

    @property
    def subdataset_names(self):
        """HDF subdataset name"""
        __subdataset_names = []
        for dts_name in self.dataset_filenames:
            sdts_name = 'HDF4_EOS:EOS_GRID:"%s":MOD_Grid_500m_Surface_Reflectance:%s' % (dts_name, self.suffix)
            __subdataset_names.append(sdts_name)

        logging.debug("Subdataset names are %s" % (__subdataset_names, ))
        return __subdataset_names

    @property
    def subdataset_prefixes(self):
        "Prefix for this band of the  based on dataset filename"
        sd_common_prefixes_list = [__dts_fname.split('.')[:5] for __dts_fname in self.dataset_filenames]
        sd_common_prefixes = ['.'.join(__list) for __list in sd_common_prefixes_list]
        sd_prefixes = ["%s_%s" % (_sd_c_prefix,  self.suffix) for _sd_c_prefix in sd_common_prefixes]
        logging.debug("Subdataset prefixes : %s" % (sd_prefixes, ) )
        return sd_prefixes

    @property
    def geotiff_names(self):
        "HDF band export as Geotiff  filaname"
        _geotiff_names = ["%s.tif" % (__prefix, ) for __prefix in self.subdataset_prefixes]
        logging.debug("Geotiff names are %s" % (_geotiff_names, ))
        return _geotiff_names


    @property
    def translated(self):
        "If the dataset is translated, we will not translate it again"
        return self._translated

    @translated.setter
    def translated(self, status):
        self._translated = status

    @property
    def mosaiced(self):
        "If mosaic is built, we will not build it again"
        return self._mosaiced

    @mosaiced.setter
    def mosaiced(self,  status):
        self._mosaiced = status

    @property
    def mosaic_name(self):
        return os.path.join( "MYD09A1","%s_%s_mosaic.tif" % (self.prefix, self.suffix, ))

    @property
    def reprojected_name(self):
        "Name of the reprojected file"
        return os.path.join( "MYD09A1","%s_%s_%s.tif" % (self.prefix, self.suffix, self.srs_suffix, ))

    @property
    def clipped_name(self):
        "Name of the clipped file"
        return os.path.join( "MYD09A1","%s_%s_%s_clipped.tif" % (self.prefix, self.suffix, self.srs_suffix, ))

    @property
    def filled_name(self):
        "Name of the filled file"
        return os.path.join( "MYD09A1","%s_%s_%s_filled.tif" % (self.prefix, self.suffix, self.srs_suffix, ))

    @property
    def  pgm_name(self):
        "Name of the pgm"
        return os.path.join( "MYD09A1","%s_%s_%s.pgm" % (self.prefix, self.suffix, self.srs_suffix,  ))

    @property
    def reprojected(self):
        "If data is reprojected, we will not reproject it again"
        return self._reprojected

    @reprojected.setter
    def reprojected(self,  status):
        self._reprojected = status

    @property
    def clipped(self):
        "If data is clipped, we will not clip it again"
        return self._clipped

    @clipped.setter
    def clipped(self, status):
        self._clipped = status

    @property
    def filled(self):
        "If data is filled, we will not fill it again"
        return self._filled

    @filled.setter
    def filled(self,  status):
        self._filled = status

    @property
    def saved_as_pgm(self):
        "If data as already been saved as pgm, we will not save it again"
        return self._saved_as_pgm

    @saved_as_pgm.setter
    def saved_as_pgm(self,  status):
        self._saved_as_pgm = status

    def download(self):
        "TODO: Download required Modis reflectance data"
        pass

    def translate(self):
        logging.debug("Modis extraction configuration...")
        for (__dataset,  __dstfile) in zip(self.subdataset_names, self.geotiff_names):
            logging.info("Translating %s..." % (__dataset, ))
            logging.info("...to %s" % (__dstfile, ))
            band_translate = Translate()
            band_translate.nodata = self.nodata
            band_translate.dataset = __dataset
            band_translate.dstfile = __dstfile
            band_translate.translate()

    def mosaic(self):
        logging.debug("Modis mosaic configuration...")
        modis_mosaic = Mosaic()
        modis_mosaic.nodata = self.nodata
        modis_mosaic.tiles = self.geotiff_names
        logging.debug("Setting Modis mosaic name to %s" % (self.mosaic_name, ))
        modis_mosaic.filename = self.mosaic_name
        logging.info("Building using Modis mosaic")
        #logging.log(VERBOSE_LVL,  srtm_mosaic.build_cmd)
        modis_mosaic.build()

    def reproject(self):
        modis_reproject = Reproject()
        modis_reproject.nodata = self.nodata
        modis_reproject.srcfile = self.mosaic_name
        modis_reproject.dstfile = self.reprojected_name
        modis_reproject.proj4string = self.params.proj4string
        logging.info("Reprojectiong Modis mosaic")
        modis_reproject.reproject()

    def clip(self):
        modis_clip = Clip()
        modis_clip.nodata = self.nodata
        modis_clip.srcfile = self.reprojected_name
        modis_clip.dstfile = self.clipped_name
        modis_clip.pixsize = self.params.pixsize
        modis_clip.bbox = self.params.bbox
        modis_clip.clip()

    def fill(self):
        modis_fill = Fill()
        modis_fill.nodata = self.nodata
        modis_fill.srcfile = self.clipped_name
        modis_fill.dstfile = self.filled_name
        modis_fill.fill()

    def save_pgm(self):
        logging.info("Saving SRTM as PGM")
        modis_pgm = Pgm()
        modis_pgm.extra_scale_factor = 0.0001
        modis_pgm.extra_offset = 0.0
        modis_pgm.srcfile = self.filled_name
        modis_pgm.dstfile = self.pgm_name
        modis_pgm.write()

    def prepare(self):
        "Complete remaining operations on band dataset"
        if not self.translated:
            logging.info("Translating MODIS band %s" % (self.suffix, ))
            self.translate()
        if not self.mosaiced:
            logging.info("Mosaicinq MODIS band %s" % (self.suffix, ))
            self.mosaic()
        if not self.reprojected:
            logging.info("Reprojecting MODIS band %s" % (self.suffix, ))
            self.reproject()
        if not self.clipped:
            logging.info("Clipping MODIS band %s" % (self.suffix, ))
            self.clip()
        if not self.filled:
            logging.info("Fillinq MODIS band %s" % (self.suffix, ))
            self.fill()
        if not self.saved_as_pgm:
            logging.info("Saving MODIS band %s as pgm" % (self.suffix, ))
            self.save_pgm()
        logging.info("Band %s is ready!" % (self.suffix, ))
        
def water_state(state_flag):
    
    water_flag = bin(state_flag)[-6:-3]
    
    # shallow ocean
    if water_flag == '000':
        water = 0
    # land
    elif water_flag == '001':
        water = 1
    # ocean coastlines and lake shorelines
    elif water_flag == '010':
        water = 2
    # shallow inland water
    elif water_flag == '011':
        water = 3
    # ephemeral water
    elif water_flag == '100':
        water = 4
    # deep inland water
    elif water_flag == '101':
        water = 5
    # continental/moderate ocean
    elif water_flag == '110':
        water = 6
    # deep ocean
    elif water_flag == '111':
        water = 7
    else:
        raise ValueError('Bad water state flag')
        
    return water
    
# version vectorisee pour operation matricielle
water_state_v = numpy.vectorize(water_state, otypes=[numpy.int])
        

class MYD09A1(object):
    """MODIS 8 days reflectance product"""

    def __init__(self):
        self._bands = []
        self.prefix = None
        self.translated = False
        self.mosaiced = False
        self.reprojected = False
        self.clipped = False
        self.filled = False
        self.saved_as_pgm = False

    @property
    def bands(self):
        return self._bands

    @bands.setter
    def bands(self,  band_suffixes):
        for b in band_suffixes:
            __c_band = Band()
            __c_band.prefix = self.prefix
            __c_band.suffix = b
            __c_band.dataset_filenames = self.filenames
            logging.debug("Is band %s translated?: %s" % (b,  str(self.translated)))
            __c_band.translated = self.translated
            __c_band.mosaiced = self.mosaiced
            __c_band.clipped = self.clipped
            __c_band.reprojected = self.reprojected
            __c_band.filled = self.filled
            __c_band.saved_as_pgm = self.saved_as_pgm
            self._bands.append(__c_band)

    def prepare(self):
        for b in self.bands:
            b.prepare()


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    feb_dataset = MYD09A1()
    feb_dataset.prefix = "february"
    feb_dataset.filenames = [os.path.join("MYD09A1", "MYD09A1.A2010041.h17v04.005.2010051022235.hdf"),
                               os.path.join("MYD09A1","MYD09A1.A2010041.h17v05.005.2010050230516.hdf")]
    feb_dataset.translated = False
    feb_dataset.mosaiced = False
    feb_dataset.reprojected = False
    feb_dataset.clipped = False
    feb_dataset.filled = False
    feb_dataset.saved_as_pgm = False
    #TODO: Bands must be added at the end...
    feb_dataset.bands = ["sur_refl_b01", "sur_refl_b03", "sur_refl_b04"]
    feb_dataset.prepare()

    aug_dataset = MYD09A1()
    aug_dataset.prefix = "august"
    aug_dataset.filenames = [os.path.join("MYD09A1","MYD09A1.A2010217.h17v04.005.2010244182438.hdf"),
                             os.path.join("MYD09A1","MYD09A1.A2010217.h17v05.005.2010244102916.hdf")]
    aug_dataset.translated = False
    aug_dataset.mosaiced = False
    aug_dataset.reprojected = False
    aug_dataset.clipped = False
    aug_dataset.filled = False
    aug_dataset.saved_as_pgm = False
    #TODO: Bands must be added at the end...
    aug_dataset.bands = ["sur_refl_b01", "sur_refl_b03", "sur_refl_b04"]
    aug_dataset.prepare()

