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
parser.add_argument( "exec_dir", default='.', nargs='?', help="Execution directory." )
parser.add_argument( '-d', '--domain', dest="params_filename",
    help="Domain definition file [domain.ini]. If present, extract contribution maps." )
parser.add_argument( '-p', '--param', action='append', nargs=2, default=[],
    metavar=('NAME','VALUE(S)'), help="Values of the parameter NAME to extract. Multiple values must be separated by commas.")

p = parser.parse_args()

for dirpath,dirnames,filenames in os.walk(p.exec_dir):
    out_names = filter(lambda fname: fname.endswith(".out") and \
        not fname.endswith(".mie.out"), filenames)
    if len(out_names) == 0:
        continue
    try:
        for pname,vals in p.param:
            if pname not in dirpath:
                print "ERROR: Parameter '%s' not found." % pname
                exit()
            for val in vals.split(','):
                if "%s_%s" % (pname,val) in dirpath:
                    break
            else:
                raise ValueError()
    except ValueError:
        continue

    for oname in out_names:
        with open(os.sep.join([dirpath,oname])) as f:
            lines = f.readlines()

        path = dirpath.split("exec")[-1][1:]
        val = float(lines[-4])
        print path, val
