#!/usr/bin/env python
import mapscript

from osgeo import gdal
from osgeo.gdalconst import *

req = mapscript.OWSRequest()
req.loadParams()


layer_name = req.getValueByName('LAYERS') or req.getValueByName('LAYER')


empty_class_id = 0
no_empty_class_id = range(1,8)

d = gdal.Open('../data/%s.tif' % (layer_name,), GA_ReadOnly)
b = d.GetRasterBand(1)
(b_min,b_max) = b.ComputeRasterMinMax(0)

buckets_size = (b_max - b_min) / 7.0
buckets_bounds = [b_min + i * buckets_size for i in range(8)]
buckets_mids = [b_min + (0.5 + i) * buckets_size for i in range(7)]

min_map = b_min + (b_max - b_min) / 100.0

map = mapscript.mapObj( '../map/1layer.map' )
l = map.getLayer(0)


l.data = layer_name + '.tif'
l.name = layer_name
#l.addProcessing("SCALE=%f,%f" % (b_min, b_max))
#l.addProcessing("SCALE_BUCKET=8")

for class_id in range(8):
    c_class = l.getClass(class_id)
    if class_id == 1:
        c_class.name = "%.2e" % (buckets_mids[class_id-1],)
    elif class_id == 7:
        c_class.name = "%.2e" % (buckets_mids[class_id-1],)
    if class_id == 0:
        c_class.setExpression("([pixel] <= %f)" % (min_map,))
    else:
        c_class.setExpression("([pixel] <= %f)" % (buckets_bounds[class_id],))
map.OWSDispatch( req )



