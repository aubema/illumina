import unittest
import modis

class FilenameTest(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None
        self.dts = modis.ReflectanceDataset(['MYD09A1.A2010033.h17v04.005.2010043020424.hdf', 
                                             'MYD09A1.A2010033.h17v05.005.2010043002711.hdf'])
        
    def test_tif_names(self):
        self.a_tif_names = self.dts.subdataset_tif_names('MYD09A1.A2010033.h17v04.005.2010043020424.hdf')
        self.e_tif_names = ['MYD09A1.A2010033.h17v04.005.2010043020424_sur_refl_b01.tif',
                            'MYD09A1.A2010033.h17v04.005.2010043020424_sur_refl_b02.tif',
                            'MYD09A1.A2010033.h17v04.005.2010043020424_sur_refl_b03.tif',
                            'MYD09A1.A2010033.h17v04.005.2010043020424_sur_refl_b04.tif',
                            'MYD09A1.A2010033.h17v04.005.2010043020424_sur_refl_b05.tif',
                            'MYD09A1.A2010033.h17v04.005.2010043020424_sur_refl_b06.tif',
                            'MYD09A1.A2010033.h17v04.005.2010043020424_sur_refl_b07.tif',
                            'MYD09A1.A2010033.h17v04.005.2010043020424_sur_refl_qc_500m.tif',
                            'MYD09A1.A2010033.h17v04.005.2010043020424_sur_refl_state_500m.tif',
                            ]
        self.assertItemsEqual(self.a_tif_names,  self.e_tif_names)
        
        
    def test_mosaic_names(self):
        self.a_mosaic_names = self.dts.subdataset_mosaic_names()
        self.e_mosaic_names = ['MYD09A1.A2010033_sur_refl_b01.tif',
                               'MYD09A1.A2010033_sur_refl_b02.tif',
                               'MYD09A1.A2010033_sur_refl_b03.tif',
                               'MYD09A1.A2010033_sur_refl_b04.tif',
                               'MYD09A1.A2010033_sur_refl_b05.tif',
                               'MYD09A1.A2010033_sur_refl_b06.tif',
                               'MYD09A1.A2010033_sur_refl_b07.tif',
                               'MYD09A1.A2010033_sur_refl_qc_500m.tif',
                               'MYD09A1.A2010033_sur_refl_state_500m.tif',
                            ]
        
if __name__ == '__main__':
    unittest.main()
