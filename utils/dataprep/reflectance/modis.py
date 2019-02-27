#!/usr/bin/env python2
# -*- coding: utf-8 -*-
""""Modis reflectance data extraction

:Author: Jean-Denis Giguère
:License: GPLv3 or later
"""

import logging
import subprocess

class ReflectanceDataset:
        
    def __init__(self,  tilesnames):
        self.__subdataset_identifier = ['sur_refl_b01',  'sur_refl_b02',  
                                        'sur_refl_b03', 'sur_refl_b04',  
                                        'sur_refl_b05',  'sur_refl_b06', 
                                        'sur_refl_b07',  'sur_refl_qc_500m',
                                        'sur_refl_state_500m']
        self.tiles = tilesnames
        
    def subdataset_names(self,  filename):
        sd_names = []
        for sd_id in self.__subdataset_identifier:
            sd_name = 'HDF4_EOS:EOS_GRID:"%s":MOD_Grid_500m_Surface_Reflectance:%s' % (filename,  sd_id)
            sd_names.append(sd_name)
        return sd_names
        
    def subdataset_prefixes(self, filename):
        sd_prefixes = []
        sd_common_prefix_part =  filename.split('.')
        sd_common_prefix = '.'.join(sd_common_prefix_part[:5])
        for sd_id in self.__subdataset_identifier:
            sd_prefix = '%s_%s' % (sd_common_prefix,  sd_id)
            sd_prefixes.append(sd_prefix)
        return sd_prefixes
        
    def subdataset_tif_names(self,  filename):
        return [sdp + '.tif' for sdp in self.subdataset_prefixes(filename)]
        
    def subdataset_tif_names_all(self):
        tif_names_by_fname = []
        for tile in self.tiles:
            tile_tifs = self.subdataset_tif_names(tile)
            tif_names_by_fname.append(tile_tifs)
        return tif_names_by_fname
        
    def subdataset_mosaic_names(self):
        fname1 = self.tiles[0]
        fname1_part = fname1.split('.')
        mosaic_prefix = '.'.join(fname1_part[0:2])
        return ["%s_%s.tif" % (mosaic_prefix , sdi ) for sdi in self.__subdataset_identifier]
        
    def reprojected_names(self):
        return ["%s_%s.tif" % (m.strip('.tif'),  self.srs_suffix) for m in self.subdataset_mosaic_names()]
        
    def clipped_names(self):
        return ["%s_%s_clipped.tif" % (m.strip('.tif'),  self.srs_suffix) for m in self.subdataset_mosaic_names()]
        
    def extract_band_tif(self):
        for tile in self.tiles:
            for (sd_name,  sd_tif) in zip(self.subdataset_names(tile),  self.subdataset_tif_names(tile)):
                logging.info("Extracting %s to %s" % (sd_name,  sd_tif))
                extract_command = ['gdal_translate',  sd_name,  sd_tif]
                subprocess.call(extract_command)
                
    def mosaic(self):
        for mosaic_file_args in zip(self.subdataset_mosaic_names(),  *self.subdataset_tif_names_all()):
            
            _mosaic_name = mosaic_file_args[0]
            _pre_mosaic_name = 'pre_' + _mosaic_name
            merge_command = ['gdal_merge.py', '-o', _pre_mosaic_name,  '-of', 'GTiff', 
            '-n', str(self.nodata)]
            merge_command = merge_command + list(mosaic_file_args[1:])
            subprocess.call(merge_command)
            
            set_nodata_command = ['gdalwarp', '-srcnodata', str(self.nodata),
                    '-dstnodata', str(self.nodata), _pre_mosaic_name, _mosaic_name]
            subprocess.call(set_nodata_command)
    def reproject(self):
        for (_mosaic,  _reproj_mosaic) in zip(self.subdataset_mosaic_names(),  self.reprojected_names()):
            reproject_cmd = ['gdalwarp', '-t_srs', self.srs, _mosaic, _reproj_mosaic]
            subprocess.call(reproject_cmd)

    def clip(self):
        for (_reproj_mosaic,  _clipped_mosaic) in zip(self.reprojected_names(), self.clipped_names()):
            clip_cmd = ['gdalwarp', '-te'] + \
                       [str(coord) for coord in self.bbox] + \
                       ['-tr', str(self.pixsiz), str(self.pixsiz)] + \
                       [_reproj_mosaic, _clipped_mosaic]
            subprocess.call(clip_cmd)

            
 

    
if __name__ == '__main__':
    logging.basicConfig(filename='reflectance_extractor.log',  level=logging.INFO)
    modis_tiles = ['MYD09A1.A2010033.h17v04.005.2010043020424.hdf', 'MYD09A1.A2010033.h17v05.005.2010043002711.hdf']
    dataset = ReflectanceDataset(modis_tiles)
    dataset.nodata = -28672 # From MYD09GA documentation
    dataset.bbox = [239981.58663178002, 4324182.311412051, 639981.58663178, 4624182.311412051]
    dataset.srs = 'epsg:23030'
    dataset.srs_suffix = 'utm30'
    dataset.pixsiz = 1000.0
    #dataset.extract_band_tif()
    #dataset.mosaic()
    #dataset.reproject()
    dataset.clip()

