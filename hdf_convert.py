#!/usr/bin/env python

import MultiScaleData as MSD
import argparse
import geopandas as gpd
import numpy as np

parser = argparse.ArgumentParser(
    description="Convert an Illumina HDF file to a georeferenced format."
)
parser.add_argument( '-f', '--format', default='vector',
    choices=['vector','raster'], help="Output type." )
parser.add_argument( '-log', action='store_true', help="If set, outputs the log10 of the data.")
parser.add_argument( '-area', action='store_true', help="If set, outputs the data in units per area [km^2]." )
parser.add_argument( 'filename', help="Input file name. Must be HDF format." )
parser.add_argument( 'outname', help="Output file name. The extension is added automatically.")

p = parser.parse_args()

hdf = MSD.Open(p.filename)
hdf.set_buffer(-1)
hdf.set_overlap(-1)

if p.format == 'raster':
    import gdal, osr
    driver = gdal.GetDriverByName("GTiff")
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(int(hdf._attrs['srs'].split(':')[1]))

    for l,data in enumerate(hdf):
        pix_size = hdf._attrs['layers'][l]['pixel_size']
        b = hdf._attrs['layers'][l]['buffer']
        xmin = hdf._attrs['layers'][l]['xmin'] + b*pix_size
        ymax = hdf._attrs['layers'][l]['ymax'] - b*pix_size
        data = data[b:-b,b:-b]
        if p.area:
            data /= (hdf.pixel_size(l)/1000.)**2
        if p.log:
            data = np.log10(data)
        ds = driver.Create(
            p.outname+"_%d.tif" % l,
            data.shape[1],
            data.shape[0],
            1,
            gdal.GDT_Float64
        )
        ds.SetProjection(srs.ExportToWkt())
        ds.SetGeoTransform((xmin,pix_size,0,ymax,0,-pix_size))
        ds.GetRasterBand(1).WriteArray(data)
        ds.GetRasterBand(1).SetNoDataValue(-np.inf if p.log else -1)

elif p.format == 'vector':
    points = {'x':[],'y':[],'val':[]}
    for l,data in enumerate(hdf):
        xmin = hdf._attrs['layers'][l]['xmin']
        ymin = hdf._attrs['layers'][l]['ymin']
        pix_size = hdf._attrs['layers'][l]['pixel_size']

        pts = np.where(hdf[l]!=-1)
        points['x'].extend( (pts[1]+0.5)*pix_size + xmin )
        points['y'].extend( (pts[0]+0.5)*pix_size + ymin )
        data = hdf[l][hdf[l]!=-1]
        if p.area:
            data /= (hdf.pixel_size(l)/1000.)**2
        if p.log:
            data = np.log10(data)
        points['val'].extend(data)

    gdf = gpd.GeoDataFrame(
        points,
        crs=hdf._attrs['srs'],
        geometry=gpd.points_from_xy(points['x'],points['y'])
    )

    gdf.to_file(p.outname+'.geojson', driver='GeoJSON')
