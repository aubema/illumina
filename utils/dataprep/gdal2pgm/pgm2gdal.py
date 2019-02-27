#!/usr/bin/env python2
import logging

import sys

import numpy
from matplotlib import pyplot, image

import osr
import gdal
from gdalconst import *

from pyproj import Proj

def parse_header_lines(header_lines):
    header_info = {}
    # Warning, latitude and longitude are mixed up!
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
    
def gen_GTiff(pgm_lines, gtiff_fname):
    pgm_header_info = parse_header_lines(pgm_lines[0:7])
    pgm_data_array = gen_array(pgm_lines)
    
    #Select driver
    driver = gdal.GetDriverByName('GTiff')
    
    #Create file
    dataset = driver.Create( gtiff_fname, pgm_header_info['width'], 
        pgm_header_info['height'], 1, gdal.GDT_Int32)
    
    epsg4083 = Proj(init='epsg:4083')
    logging.debug('Lower left coordinate (long, lat) : (%f, %f)' % (pgm_header_info['lon0'], pgm_header_info['lat0']))
    x_bl, y_bl = epsg4083(pgm_header_info['lon0'], 
            pgm_header_info['lat0'])
    logging.debug('Lower left coordinates (x, y): (%f, %f)' % (x_bl, y_bl))
    x0 = x_bl
    # Check reference
    y0 = y_bl
    
    dataset.SetGeoTransform([x0, pgm_header_info['pixsize'], 0,
    y0, 0, -1 * pgm_header_info['pixsize'] ])
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4083)
    
    dataset.SetProjection(srs.ExportToWkt())
    dataset.GetRasterBand(1).WriteArray(pgm_data_array)
    
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
    
if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    pgm_fname = sys.argv[1]
    logging.info('Converting %s ...' % (pgm_fname,))
    pgm_file = open(pgm_fname, 'r')
    pgm_lines = pgm_file.readlines()
    header_lines = pgm_lines[0:7]
    pgm_header_info = parse_header_lines(header_lines)
    pgm_data_array = gen_array(pgm_lines)
    pgm_fname_prefix = pgm_fname[:pgm_fname.rfind('.')]
    gtiff_fname = pgm_fname_prefix + '.tif'
    gen_GTiff(pgm_lines, gtiff_fname)
    
    
