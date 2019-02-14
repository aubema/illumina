#!/usr/bin/env python
#
# Illumina output extract
#
# Author : Alexandre Simoneau
# unless noted otherwise
#
# February 2019

import os
import argparse

parser = argparse.ArgumentParser(description="Extract Illumina output.")
parser.add_argument( "exec_dir", help="Execution directory." )
parser.add_argument( "-d, --domain", dest="params_filename",
    help="Domain definition file [domain.ini]. If present, extract contribution maps." )

p = parser.parse_args()

failed = []
for dirpath,dirnames,filenames in os.walk(p.exec_dir):
    out_names = filter(lambda fname: fname.endswith(".out") and \
        not fname.endswith(".mie.out"), filenames)
    if len(out_names) == 0:
        continue
    for oname in out_names:
        with open(os.sep.join([dirpath,oname])) as f:
            lines = f.readlines()

        path = dirpath.split("exec")[-1][1:]
        try:
            val = float(lines[-4])
        except:
            failed.append(path)
        else:
            print path, val

print "\nFailed:"

for path in failed:
    print path
