#!/usr/bin/env python
"""Generate mapfile with an entry for every geotif in the directory
"""

import glob

layer_template = """LAYER
    NAME "%(title)s"
    METADATA
      "wms_title"  "%(title)s"
    END
    TYPE RASTER
    DATA %(fname)s
    STATUS ON
    PROJECTION
      "init=epsg:900913"
    END
    CLASSITEM "[pixel]"
    CLASS
      EXPRESSION "0"
      STYLE
        COLOR 255 255 255
      END
    END
    CLASS
      EXPRESSION ([pixel] > 0 AND [pixel] < 65536)
      STYLE
        COLORRANGE 255 255 255 255 0 0
        DATARANGE 1 65535
      END
    END
  END
"""

def dict_from_fname(fname):
    fname_dict = {}
    fname_dict['title'] = fname[:fname.rfind('.')]
    fname_dict['fname'] = fname
    return(fname_dict)
    
if __name__ == '__main__':

    mapfile = open('layers.map', 'w')
    for fname in glob.glob('*.tif'):
        f_dict = dict_from_fname(fname)
        f_str = layer_template % f_dict
        mapfile.writelines(f_str)
    mapfile.close()
    

