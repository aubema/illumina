"""Spatial domain related tests

:Author: Jean-Denis Giguere
"""

import unittest
import sys

sys.path.insert(1, '..')


import gdal2pgm

class BBoxSpatialParamsConversionTestCase(unittest.TestCase):
    """Bounding box to and from Spatial parameters conversion

    There are two ways to represent spatial domain:
        * Bounding box : minimum and maximum coordonates of spatial
        in form minx,miny,maxx,maxy and pixel size
        * Spatial parameters : Latitude and longitude of
        south-west corner of domain, width and height of raster
        and pixel size
    """

    def setUp(self):
        self.bbox = (740305, 5016940, 745742, 5022037)
        self.srs = '+init=epsg:2959'
        self.pixsiz = 30
        self.lat0 = 45.26487
        self.lon0 = -71.936879
        self.width = 182
        self.height = 170
        self.bbox_alt = (-734391,443604,-659531,452892)
        self.srs_alt = '+init=epsg:32198'
        self.pixsiz_alt = 150
        self.lat0_alt=47.53526
        self.lon0_alt=-78.311929
        self.width_alt = 500
        self.height_alt = 62

    def testBBoxToSpatialParams(self):
        self.bboxSpatialParams = gdal2pgm.SpatialParams(
            bbox=self.bbox, srs=self.srs, pixsiz=self.pixsiz)
        self.computedSpatialParams = self.bboxSpatialParams.compute()
        self.assertAlmostEqual(self.computedSpatialParams.lat0, self.lat0,5)
        self.assertAlmostEqual(self.computedSpatialParams.lon0, self.lon0,5)
        self.spatialParams_alt = gdal2pgm.SpatialParams(
            bbox=self.bbox_alt, srs=self.srs_alt, pixsiz=self.pixsiz_alt)
        self.computedSpatialParams_alt = self.spatialParams_alt.compute()
        self.assertAlmostEqual(self.computedSpatialParams_alt.lat0, 
                               self.lat0_alt,5)
        self.assertAlmostEqual(self.computedSpatialParams_alt.lon0, 
                               self.lon0_alt,5)

    def testSpatialParamsToBBox(self):
        self.defSpatialParams = gdal2pgm.SpatialParams(
            srs=self.srs, lat0=self.lat0, lon0=self.lon0, pixsiz=self.pixsiz,
            width=self.width, height=self.height)
        self.defSpatialParamsComputed = self.defSpatialParams.compute()
        self.assertEqual(len(self.defSpatialParamsComputed.bbox),4)
        (self.minx, self.miny, self.maxx, self.maxy) = \
                self.defSpatialParamsComputed.bbox
        self.assertAlmostEqual(self.minx,self.bbox[0],1)
        self.defSpatialParams_alt = gdal2pgm.SpatialParams(
            srs=self.srs_alt, lat0=self.lat0_alt, lon0=self.lon0_alt, 
            pixsiz=self.pixsiz_alt, width=self.width, height=self.height)
        self.defSpatialParamsComputed_alt = \
                self.defSpatialParams_alt.compute()
        (self.minx_alt, self.miny_alt, self.maxx_alt, self.max_alt) = \
                self.defSpatialParamsComputed_alt.bbox
        self.assertAlmostEqual(self.minx_alt,self.bbox_alt[0],1)


if __name__ == '__main__':
    unittest.main()
