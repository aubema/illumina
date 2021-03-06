#!/usr/bin/env python3

import shutil
import os
from glob import glob
import click

@click.command()
def init():
    """Initialize an execution folder."""
    print("Initializing Illumina execution folder.")
    ppath = os.environ['PATH'].split(os.pathsep)
    illumpath = [s for s in ppath if "illumina" in s and "bin" not in s][0]

    if not os.path.exists("Lights"):
        shutil.copytree(illumpath+"/Example/Lights",'Lights')

    example_files = glob(illumpath+"/Example/*.*")
    files = [s for s in example_files if not(s.endswith("hdf5") or s.endswith("csv"))]

    for filename in files:
        shutil.copy2(filename,'.')
