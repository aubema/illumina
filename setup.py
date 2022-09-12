#!/usr/bin/env python3

from subprocess import check_output

from setuptools import setup

with open("illum/__init__.py") as f:
    info = {}
    for line in f:
        if line.startswith("__version__"):
            exec(line, info)
            break

gdal_version = (
    check_output(["gdalinfo", "--version"]).decode().split(",")[0].split()[1]
)

setup(
    name="illum",
    version=info["__version__"],
    packages=["illum"],
    install_requires=[
        "astropy",
        "Click",
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
        "pyyaml",
        "scipy",
        "xmltodict",
    ],
    extras_require={
        "dev": [
            "black",
            "flake8",
            "ipython",
            "isort",
        ],
    },
    entry_points="""
        [console_scripts]
        illum=illum.main:illum
    """,
)
