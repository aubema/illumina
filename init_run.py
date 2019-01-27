#!/usr/bin/env python

import shutil
import os
from glob import glob

ppath = os.environ['PATH'].split(os.pathsep)
illumpath = filter(
    lambda s: "illumina" in s and "bin" not in s,
    ppath)[0]

if not os.path.exists("Lights"):
    shutil.copytree(illumpath+"/Example/Lights",'Lights')

example_files = glob(illumpath+"/Example/*.*")
files = filter(
    lambda s: not(s.endswith("bin") or s.endswith("csv")),
    example_files)

for filename in files:
    shutil.copy2(filename,'.')
