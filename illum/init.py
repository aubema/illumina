#!/usr/bin/env python3

import os
import shutil
from glob import glob
from importlib.resources import path

import click


@click.command(name="init")
def CLI_init():
    """Initialize an execution folder."""
    init()


def init():
    print("Initializing Illumina execution folder.")
    illumpath = path("illum", "data").as_posix()

    if not os.path.exists("Lights"):
        shutil.copytree(illumpath + "/Example/Lights", "Lights")

    example_files = glob(illumpath + "/Example/*.*")
    files = [
        s
        for s in example_files
        if not (s.endswith("hdf5") or s.endswith("csv"))
    ]

    for filename in files:
        shutil.copy2(filename, ".")
