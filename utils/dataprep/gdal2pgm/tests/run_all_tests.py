"""run_all_tests.py - Run every test cases

:Author: Jean-Denis Giguere
"""

import unittest
import sys

import test_spatial
import test_parameters

if __name__ == '__main__':
    allTests = unittest.TestSuite()
    spatialTests = unittest.defaultTestLoader.loadTestsFromModule(
        test_spatial)
    parameterTests = unittest.defaultTestLoader.loadTestsFromModule(
        test_parameters)
    allTests.addTests([spatialTests,
                       parameterTests,
                      ])
    unittest.main(defaultTest='allTests')
