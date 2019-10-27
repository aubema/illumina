#!/usr/bin/env python

import MultiScaleData as MSD
import gdal, osr
import argparse
import geopandas as gpd
import numpy as np

parser = argparse.ArgumentParser(
    description="Convert an Illumina HDF file to a georeferenced format."
)
parser.add_argument( '-f', '--format', default='vector',
    choices=['vector','raster'], help="Output type." )
parser.add_argument( 'filename', help="Input file name. Must be HDF format." )
parser.add_argument( 'outname', help="Output file name. The extension is added automatically.")

p = parser.parse_args()

hdf = MSD.Open(p.filename)

if p.format == 'raster':
    driver = gdal.GetDriverByName("GTiff")
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(int(hdf._attrs['srs'].split(':')[1]))

    for l,data in enumerate(hdf):
        xmin = hdf._attrs['layers'][l]['xmin']
        ymax = hdf._attrs['layers'][l]['ymax']
        pix_size = hdf._attrs['layers'][l]['pixel_size']
        ds = driver.Create(
            p.outname+"_%d.tif" % l,
            data.shape[1],
            data.shape[0],
            1,
            gdal.GDT_Float64
        )
        ds.SetProjection(srs.ExportToWkt())
        ds.SetGeoTransform((xmin,pix_size,0,ymax,0,-pix_size))
        ds.GetRasterBand(1).WriteArray(data[::-1])

elif p.format == 'vector':
    points = {'x':[],'y':[],'val':[]}
    for l,data in enumerate(hdf):
        xmin = hdf._attrs['layers'][l]['xmin']
        ymin = hdf._attrs['layers'][l]['ymin']
        pix_size = hdf._attrs['layers'][l]['pixel_size']

        pts = np.where(hdf[l])
        points['x'].extend( (pts[1]+0.5)*pix_size + xmin )
        points['y'].extend( (pts[0]+0.5)*pix_size + ymin )
        points['val'].extend(hdf[l][hdf[l]!=0])

    gdf = gpd.GeoDataFrame(
        points,
        crs={'init':hdf._attrs['srs']},
        geometry=gpd.points_from_xy(points['x'],points['y'])
    )

    gdf.to_file(p.outname+'.geojson', driver='GeoJSON')
