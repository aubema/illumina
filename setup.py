from subprocess import check_output

import illum
from setuptools import setup

gdal_version = (
    check_output(["gdalinfo", "--version"]).decode().split(",")[0].split()[1]
)

setup(
    name="illum",
    version=illum.__version__,
    packages=["illum"],
    install_requires=[
        "astropy",
        "Click<8",
        "fiona",
        "gdal==" + gdal_version,
        "geopandas",
        "GitPython",
        "h5py",
        "matplotlib",
        "numpy",
        "osmnx",
        "pandas",
        "pillow",
        "progressbar2",
        "pyproj",
        "PyQt5",
        "pyyaml",
        "scipy",
    ],
    entry_points="""
        [console_scripts]
        illum=illum.main:illum
    """,
)
