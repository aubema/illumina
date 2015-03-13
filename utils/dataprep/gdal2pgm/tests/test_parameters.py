"""Basic parameters testing

:Author: Jean-Denis Giguere


"""

import unittest

import subprocess

gdal2pgm = ['../gdal2pgm.py',]

class ParameterTestCase(unittest.TestCase):
    def testNoParams(self):
        """Calling gdal2pgm without arguments"""
        self.retcode = subprocess.call(gdal2pgm)
        self.assertEqual(self.retcode, 101, 
                        'Wrong return code for empty argument list')

    def testNoPixsiz(self):
        """Calling gdal2pgm without pixsiz argument"""
        self.args = ['--lat0=49.0','--lon0=-72.1']
        self.retcode = subprocess.call(gdal2pgm + self.args)
        self.assertEqual(self.retcode, 102,
                         'Pixel size is mandatory')

if __name__ == '__main__':
    unittest.main()
