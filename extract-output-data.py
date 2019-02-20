#!/usr/bin/env python
#
# Illumina output extract
#
# Author : Alexandre Simoneau
# unless noted otherwise
#
# February 2019

import os, re
import argparse
from collections import defaultdict as ddict
from hdftools import MSD, from_domain
from functools import partial
from pytools import load_bin

parser = argparse.ArgumentParser(description="Extract Illumina output.")
parser.add_argument( "exec_dir", default='.', nargs='?', help="Execution directory." )
parser.add_argument( '-d', '--domain', dest="params_filename",
    help="Domain definition file [domain.ini]. If present, extract contribution maps." )
parser.add_argument( '-p', '--param', action='append', nargs=2, default=[],
    metavar=('NAME','VALUE(S)'), help="Values of the parameter NAME to extract. Multiple values must be separated by commas.")

p = parser.parse_args()
regex = re.compile(r'/layer_(\d+)')

data = ddict(float) if p.params_filename is None else \
    ddict(partial(from_domain,p.params_filename))

for dirpath,dirnames,filenames in os.walk(p.exec_dir):
    out_names = filter(lambda fname: fname.endswith(".out") and \
        not fname.endswith(".mie.out"), filenames)
    if len(out_names) == 0:
        continue
    try:
        for pname,pvals in p.param:
            if pname not in dirpath:
                print "ERROR: Parameter '%s' not found." % pname
                exit()
            for pval in pvals.split(','):
                if "%s_%s" % (pname,pval) in dirpath:
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
        if p.params_filename is None:
            data[regex.sub('',path)] += val
        else:
            n_layer = int(regex.search(dirpath).groups()[0])
            pcl_name = filter(lambda s: "pcl.bin" in s, filenames)[0]
            pcl_path = os.path.join(dirpath,pcl_name)
            pcl_data = load_bin(pcl_path)
            pcl_data *= val / pcl_data.sum()
            data[regex.sub('',path)][n_layer] = pcl_data

for key,val in data.iteritems():
    if p.params_filename is None:
        print key,val
    else:
        val.save(key.replace(os.sep,'-'))
        print key,sum( v.sum() for v in val )
