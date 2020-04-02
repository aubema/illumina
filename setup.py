from setuptools import setup

setup(
    name='illum',
    version='2.5',
    py_modules=['main','defineDomain','Illuminutils','make_inputs',
    'make_zones','make_lamps','makeBATCH','find-failed-runs','extract-output-data',
    'init_run','alt_scenario_maker','hdf_convert','pytools','MultiScaleData'],
    install_requires=[
        'Click',
        'pyproj',
        'pyyaml',
        'numpy',
        'h5py',
        'pillow',
        'matplotlib',
        'scipy',
        'astropy',
        'pandas',
        'geopandas'
    ],
    entry_points='''
        [console_scripts]
        illum=main:illum
    ''',
)
