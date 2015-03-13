# -*- coding: utf-8 -*-
import unittest
import sys
import io
import datetime

# Ajout du répertoire lib
sys.path.append("../../lib")
from parseIntrantsIn import Intrants

config_file_A = """
[domaine]
date: 2012-02-15
latMin: 41.583 ;latitude minimale du domaine [degré décimal]
latMax: 44.267 ;latitude maximale du domaine [degré décimal]
lonMin: -2.266 ;longitude minimale du domaine [degré décimal]
lonMax: 2.53 ;longitude maximale du domaine [degré décimal]
taillePixel: 1000; [mètre]
"""


class TestParseDomain(unittest.TestCase):
    def setUp(self):
        
        self.intrants_A = Intrants(io.BytesIO(config_file_A))
        
    def test_date_format(self):
        self.assertIs(type(self.intrants_A.domaine.date), datetime.datetime) 
    
    def test_date_format(self):
        self._direct_date = datetime.datetime(2012,2,15)
        self.assertEqual(self.intrants_A.domaine.date, self._direct_date)
    
    def test_reversed_lon(self):
        pass
        
if __name__ == '__main__':
    unittest.main(verbosity=2)
