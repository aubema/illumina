#!/opt/python264/2.7.1/bin/python
import logging
import os
import sys
from glob import glob
import datetime
import numpy
from matplotlib import pyplot, image

from pysal.esda.mapclassify import Natural_Breaks
from colorbrewer import YlOrRd

import osr
import sys
import gdal
from gdalconst import *

from pyproj import Proj

from argparse import ArgumentParser

VERSION='1.1 build 2'

def parse_header_lines(header_lines):
    header_info = {}
    header_info['lat0'] = float(header_lines[1].split()[2])
    header_info['lon0'] = float(header_lines[2].split()[2])
    header_info['pixsize'] = float(header_lines[3].split()[2])
    header_info['gain'] = float(header_lines[4].split()[2])
    header_info['offset'] = float(header_lines[5].split()[2])
    #print header_lines[6].split()
    #print header_lines[6].split()
    header_info['width'] = int(header_lines[6].split()[0])
    header_info['height'] = int(header_lines[6].split()[1])
    header_info['greylevel'] = int(header_lines[6].split()[2])
    
    
    return header_info
    
def gen_array(pgm_lines):
    pgm_header_info = parse_header_lines(pgm_lines[0:7])
    pgm_data_lines = pgm_lines[7:]
    pgm_data = []
    for l in pgm_data_lines:
        v = [int(gl) for gl in l.split()]
        pgm_data = pgm_data + v
    data_array = numpy.array(pgm_data)
    data_array = data_array.reshape(pgm_header_info['height'], 
        pgm_header_info['width'])
        
    return data_array

    
def plot_array(pgm_data_array, filename):
    pyplot.imshow(pgm_data_array)
    pyplot.savefig(filename)
    
def gen_GTiff(pgm_lines, gtiff_fname, proj4string):
    pgm_header_info = parse_header_lines(pgm_lines[0:7])
    pgm_data_array = gen_array(pgm_lines)
    pgm_shape = pgm_data_array.shape
    val_array = pgm_data_array * float(pgm_header_info['gain']) + float(pgm_header_info['offset'])
    color_idx_array = Natural_Breaks(val_array.ravel(),8, 0).yb
    color_idx_array.shape = pgm_shape
    #Select driver
    driver = gdal.GetDriverByName('GTiff')
    
    #Create file
    dataset = driver.Create( gtiff_fname, pgm_header_info['width'], 
        pgm_header_info['height'], 1, gdal.GDT_Byte)
    epsg = Proj(proj4string)
    logging.debug('Lower left coordinate (long, lat) : (%f, %f)' % (pgm_header_info['lon0'], pgm_header_info['lat0']))
    x_bl, y_bl = epsg(pgm_header_info['lon0'], 
            pgm_header_info['lat0'])
    logging.debug('Lower left coordinates (x, y): (%f, %f)' % (x_bl, y_bl))
    x_ul = x_bl
    y_ul = y_bl + pgm_header_info['pixsize'] * pgm_header_info['height']
    x0 = x_ul
    # Check reference
    y0 = y_ul
    
    dataset.SetGeoTransform([x0, pgm_header_info['pixsize'], 0,
    y0, 0, -1 * pgm_header_info['pixsize'] ])
    srs = osr.SpatialReference()
    srs.ImportFromProj4(proj4string)
    
    dataset.SetProjection(srs.ExportToWkt())
    dataset.GetRasterBand(1).WriteArray(color_idx_array)
    dataset.GetRasterBand(1).SetColorInterpretation(gdal.GCI_PaletteIndex)
    
    color_table = gdal.ColorTable( gdal.GPI_RGB)
    ce = gdal.ColorEntry()
    ce.__dict__['this'] = (255, 255, 255, 255)
    color_table.SetColorEntry(0, ce.__dict__['this'])
    for i in range(7):
        ce = gdal.ColorEntry()
        ce.__dict__['this'] = YlOrRd[7][i]
        color_table.SetColorEntry(i+1, ce.__dict__['this'])
    dataset.GetRasterBand(1).SetColorTable(color_table)
    
    dataset = None
    
    
def test_main():
    logging.basicConfig(level=logging.DEBUG)
    #print parse_header_lines(header_lines)
    pgm_fname = 'PCL-OT-before_midnight-rd4000-ta0.200-wl615-el90-az0.pgm'
    pgm_file = open(pgm_fname, 'r')
    pgm_lines = pgm_file.readlines()
    header_lines = pgm_lines[0:7]
    pgm_header_info = parse_header_lines(header_lines)
    pgm_data_array = gen_array(pgm_lines, 'test.tif')
    #plot_array(pgm_data_array, 'test.png')
    gen_GTiff(pgm_lines)
    
def main():
    print "This is lpa_pgm2tif.py version %s" % (VERSION,)
    parser = ArgumentParser()
    parser.add_argument("proj4string")
    parser.add_argument('-f', help="Name of the pgm file to transform")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    logging.basicConfig(level=logging.INFO)
    if args.f is None:
        pgms = glob('*.pgm')
        for pgm in pgms:
             tif = pgm.rstrip('.pgm') + '.tif'
             if not os.path.exists(tif):
                 transform_pgm(pgm, args.proj4string)
    else:
        pgm = args.f
        tif = pgm.rstrip('.pgm') + '.tif'
        if not os.path.exists(tif):
            transform_pgm(pgm, args.proj4string)
    
    
def transform_pgm(pgm_fname, proj4string):
    logging.info('Conversion of %s ... starting at %s' % (pgm_fname, str(datetime.datetime.now())))
    pgm_file = open(pgm_fname, 'r')
    pgm_lines = pgm_file.readlines()
    header_lines = pgm_lines[0:7]
    pgm_header_info = parse_header_lines(header_lines)
    pgm_data_array = gen_array(pgm_lines)
    pgm_fname_prefix = pgm_fname[:pgm_fname.rfind('.')]
    gtiff_fname = pgm_fname_prefix + '.tif'
    gen_GTiff(pgm_lines, gtiff_fname, proj4string)
    logging.info('Conversion of %s ... ends at %s' % (pgm_fname, str(datetime.datetime.now())))
    
if __name__ == '__main__':
    main()
   
    
    
    
    
