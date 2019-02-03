#!/usr/bin/env python
#
# Illumina output extract
#
# Author : Alexandre Simoneau
# unless noted otherwise
#
# February 2019

import os,sys

for dirpath,dirnames,filenames in os.walk(sys.argv[1]):
    out_names = filter(lambda fname: fname.endswith(".out") and \
        not fname.endswith(".mie.out"), filenames)
    if len(out_names) == 0:
        continue
    for oname in out_names:
        with open(os.sep.join([dirpath,oname])) as f:
            lines = f.readlines()

        print dirpath.split(os.sep,2)[-1].replace(os.sep,'-'), float(lines[-4])
