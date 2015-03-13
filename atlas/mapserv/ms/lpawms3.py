#!/usr/bin/env python
"""lpawms3.py - WMS-like web mapping service for light polluation atlas

lpawms3.py is OSGeo Mapserver mapscript python script which emulate
a WMS service for the light polluation atlas (http://galileo.graphycs.cegepsherbrooke.qc.ca/atlas/)

It allows on the fly layer creation based on the illumina naming convention.

lpawms3.py also build on the symbology based on linear color scale.

This script was originally written by Jean-Denis Giguère.

"""
import mapscript

from osgeo import gdal
from osgeo.gdalconst import *

#  MapScript Wrappers for WxS Services
req = mapscript.OWSRequest()
req.loadParams()
layer_name = req.getValueByName('LAYERS') or req.getValueByName('LAYER')

# The empty class is 0, then there are 7 non-empty classes
empty_class_id = 0
no_empty_class_id = range(1,8)

# Based on the naming convention of illumina, we open the right geotiff
# We will compute the min and max value for the symbology
d = gdal.Open('../data2013/%s.tif' % (layer_name,), GA_ReadOnly)
b = d.GetRasterBand(1)
(b_min,b_max) = b.ComputeRasterMinMax(0)


# We compute the class size and find the center of each class
buckets_size = (b_max - b_min) / 7.0
buckets_bounds = [b_min + i * buckets_size for i in range(8)]
buckets_mids = [b_min + (0.5 + i) * buckets_size for i in range(7)]

#The minimal value on the map
min_map = b_min + (b_max - b_min) / 100.0

#We open a template mapscript which contain layer definition
map = mapscript.mapObj( '../map/1layer3.map' )
l = map.getLayer(0)


#We need to override the data based on illumina experiment naming convention
l.data = layer_name + '.tif'
#We also override layer name to be kind
l.name = layer_name

for class_id in range(8):
    c_class = l.getClass(class_id)
    # We change class name (for labels)
    if class_id == 1:
        c_class.name = "%.2e" % (buckets_mids[class_id-1],)
    elif class_id == 7:
        c_class.name = "%.2e" % (buckets_mids[class_id-1],)
    # We change class expression
    if class_id == 0:
        c_class.setExpression("([pixel] <= %f)" % (min_map,))
    else:
        c_class.setExpression("([pixel] <= %f)" % (buckets_bounds[class_id],))
        
# MapScript Wrappers for WxS Services
map.OWSDispatch( req )



