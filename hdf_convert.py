#!/usr/bin/env python

import MultiScaleData as MSD
import click
import geopandas as gpd
import numpy as np

@click.command()
@click.option("--vector/--raster","-v/-r",default=True,help="Output type.")
@click.option("-log",is_flag=True,help="Logarithmic scale (base 10)")
@click.option("-area",is_flag=True,help="Normalized by pixel area (in km²)")
@click.argument("filename",type=click.Path(exists=True))
@click.argument("outname")
def convert(filename,outname,vector,log,area):
    """Convert an Illumina HDF file to a georeferenced format.

    Converts FILENAME to OUTNAME.EXT where ext is defined based on the output format.
    The output format is either vector (GeoJSON) or raster (Tiff).
    """
    hdf = MSD.Open(filename)
    hdf.set_buffer(-1)
    hdf.set_overlap(-1)

    if not vector:
        from osgeo import gdal, osr
        driver = gdal.GetDriverByName("GTiff")
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(int(hdf._attrs['srs'].split(':')[1]))

        for l,data in enumerate(hdf):
            pix_size = hdf._attrs['layers'][l]['pixel_size']
            b = hdf._attrs['layers'][l]['buffer']
            xmin = hdf._attrs['layers'][l]['xmin'] + b*pix_size
            ymax = hdf._attrs['layers'][l]['ymax'] - b*pix_size
            if b != 0:
                data = data[b:-b,b:-b]
            if area:
                data /= (hdf.pixel_size(l)/1000.)**2
            if log:
                data = np.log10(data)
            ds = driver.Create(
                outname+"_%d.tif" % l,
                data.shape[1],
                data.shape[0],
                1,
                gdal.GDT_Float64
            )
            ds.SetProjection(srs.ExportToWkt())
            ds.SetGeoTransform((xmin,pix_size,0,ymax,0,-pix_size))
            ds.GetRasterBand(1).WriteArray(data)
            ds.GetRasterBand(1).SetNoDataValue(-np.inf if log else -1)

    elif vector:
        points = {'x':[],'y':[],'val':[]}
        for l,data in enumerate(hdf):
            xmin = hdf._attrs['layers'][l]['xmin']
            ymax = hdf._attrs['layers'][l]['ymax']
            pix_size = hdf._attrs['layers'][l]['pixel_size']

            pts = np.where(hdf[l]!=-1)
            points['x'].extend( (pts[1]+0.5)*pix_size + xmin )
            points['y'].extend( ymax - (pts[0]+0.5)*pix_size )
            data = hdf[l][hdf[l]!=-1]
            if area:
                data /= (hdf.pixel_size(l)/1000.)**2
            if log:
                data = np.log10(data)
            points['val'].extend(data)

        gdf = gpd.GeoDataFrame(
            points,
            crs=hdf._attrs['srs'],
            geometry=gpd.points_from_xy(points['x'],points['y'])
        )

        gdf.to_file(outname+'.geojson', driver='GeoJSON')
