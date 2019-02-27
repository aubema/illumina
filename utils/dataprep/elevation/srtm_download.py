#!/usr/bin/env python2
# -*- coding: utf-8 -*-
""""Téléchargement des données SRTM3

La définition du domaine se fait dans le module srtm_domain
SRTM tiles are archived by world region, ftp_path might need update

:Author: Jean-Denis Giguère
:License: GPLv3 or later
"""
import os
import os.path
import urllib
import zipfile
import logging

from srtm_domain import tiles

homedir = os.path.expanduser('~')
workdir = os.path.join(homedir, 'tmp', 'illumina-data') # You can use alternative directory

workdir = os.path.join(homedir, 'illumina-data',  'madrid')
srtmdir = os.path.join(workdir, 'srtm')
try:
    os.chdir(srtmdir)
except OSError:
    os.makedirs(srtmdir)
    os.chdir(srtmdir)

def tile_already_downloaded(filename):
    logging.info()
    fullpath = os.path.join(srtmdir,  filename)
    return os.path.isfile(fullpath)



# Change to world region may be required
ftp_path = 'http://dds.cr.usgs.gov/srtm/version2_1/SRTM3/Eurasia/'
srtm_suffix = '.hgt.zip'
    
def main():
    logging.basicConfig(level=logging.DEBUG)
    for tile in tiles:
        filename = tile + srtm_suffix
        if not tile_already_downloaded(filename):
            ftp_link = ftp_path + filename
            urllib.urlretrieve(ftp_link, filename)
            zip_file = zipfile.ZipFile(filename)
            srtm_file = zip_file.namelist()[0] # There is one file by archive
            srtm_file_object = open(srtm_file, 'wb')
            srtm_file_object.write(zip_file.read(srtm_file))
            srtm_file_object.close()
        
if __name__ == '__main__':
    main()
