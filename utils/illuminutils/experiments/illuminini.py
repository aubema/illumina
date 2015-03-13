# -*- coding: utf-8 -*-
""""Illumina parameter from ini config file

:Author: Jean-Denis Giguère
:Licence: GPL v3 or later
:Date: 2013-02-10
"""

import ConfigParser
import logging

class IlluminaParameters(object):
    """Illumina parameters for an experiments
     
    """
    def __init__(self, ini):
        logging.debug('Creating new Illumina parameters...')
        self.config = ConfigParser.RawConfigParser(allow_no_value=True)
        try:
            self.config.read(ini)
        except TypeError,  AttributeError:
            self.config.read(ini)
        except ConfigParser.MissingSectionHeaderError:
            logging.error("Misformed ini file: Doesn't have Section header")
        except:
            logging.error('Cannot read parameters from ini files %s' % (ini, ))
        logging.debug('Ini has sections %s' % (str(self.config.sections())))
        if self.config.has_section('domain'):
            bbox_str =self.config.get('domain', 'bbox').split(',')
            self._bbox = [float(crd) for crd in bbox_str]
            self._pixsize = float(self.config.get('domain',  'pixsize'))
            self._srssuffix = self.config.get('domain', 'srssuffix')
            self._proj4string = self.config.get('domain', 'proj4string')
            
    @property
    def bbox(self):
        """Bounding box as an array [minx, miny, maxx, maxy]"""
        return self._bbox
        
    @property
    def pixsize(self):
        "Pixel size (square) [meters]"
        return self._pixsize
    
    @property
    def srs_suffixname(self):
        "String to append to identify SRS"
        return self._srssuffix
        
    @property
    def proj4string(self):
        "Proj4 string to represent spatial reference system"
        return self._proj4string




        
