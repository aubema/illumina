from setuptools import setup
import illum

setup(
    name='illum',
    version=illum.__version__,
    packages=["illum"],
    install_requires=[
        'Click<8',
        'progressbar2',
        'pyproj',
        'pyyaml',
        'numpy',
        'h5py',
        'pillow',
        'matplotlib',
        'scipy',
        'astropy',
        'pandas',
        'geopandas',
        'fiona',
        'osmnx',
        'GitPython'
    ],
    entry_points='''
        [console_scripts]
        illum=illum.main:illum
    ''',
)
