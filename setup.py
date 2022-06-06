import illum
from setuptools import setup

setup(
    name="illum",
    version=illum.__version__,
    packages=["illum"],
    install_requires=[
        "astropy",
        "Click<8",
        "fiona",
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
