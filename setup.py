from setuptools import setup

setup(
    name='illum',
    version='2.1',
    py_modules=['main'],
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
    ],
    entry_points='''
        [console_scripts]
        illum=main:illum
    ''',
)
