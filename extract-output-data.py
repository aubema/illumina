#!/usr/bin/env python3
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
from glob import glob
from copy import deepcopy as copy
from functools import partial
import numpy as np

parser = argparse.ArgumentParser(description="Extract Illumina output.")
parser.add_argument( "exec_dir", default='.', nargs='?',
    help="Path to the input folder from where `makeBATCH` was executed." )
parser.add_argument( '-c', '--contrib', action="store_true",
    help="If present, extract contribution maps." )
parser.add_argument( '-p', '--param', action='append', nargs=2, default=[],
    metavar=('NAME','VALUE(S)'), help="Values of the parameter NAME to extract. "\
        "Multiple values must be separated by commas.")
parser.add_argument( '-f', '--full', action="store_true",
    help="If present, will extract all available outputs.")

p = parser.parse_args()
regex_layer = re.compile(r'-layer_(\d+)')
regex_coords = re.compile(r'observer_coordinates_(-?\d+\.\d+_-?\d+\.\d+)')

def MSDOpen(filename,cached={}):
    if filename in cached:
        return cached[filename]
    ds = MSD.Open(filename)
    cached[filename] = ds
    return ds

skyglow = ddict(float)
if p.full:
    outputs = ddict(partial(np.zeros,4))
if p.contrib:
    contrib = dict()

for dirpath,dirnames,filenames in os.walk(p.exec_dir):
    if not os.path.isfile(os.path.join(dirpath,'illumina.in')):
        continue

    with open(os.path.join(dirpath,'illumina.in')) as f:
        basename = f.readlines()[1].split()[0]

    out_names = [fname for fname in filenames if fname.endswith(".out") and \
        fname.startswith(basename+'_')]
    if not out_names:
        out_names = [ basename + '.out' ]

    for oname in out_names:
        if oname == basename + ".out":
            params = dirpath.split('exec'+os.sep)[1].replace(os.sep,'-')
        else:
            params = oname[len(basename)+1:-4]

        try:
            for pname,pvals in p.param:
                if pname not in params:
                    print("ERROR: Parameter '%s' not found." % pname)
                    exit()
                for pval in pvals.split(','):
                    if "%s_%s" % (pname,pval) in params:
                        break
                else:
                    raise ValueError()
        except ValueError:
            continue

        with open(os.sep.join([dirpath,oname])) as f:
            lines = f.readlines()
        idx_results = np.where([ "==" in l for l in lines ])[0][-1] + 1

        val = float(lines[-1])
        if p.full:
            vals =  np.array([ float(l) for l in lines[idx_results+1::2] ])
            outputs[regex_layer.sub('',params)] += vals
        else:
            skyglow[regex_layer.sub('',params)] += val

        if p.contrib:
            try:
                n_layer = int(regex_layer.search(params).groups()[0])
            except AttributeError:
                # No match, only 1 layer
                n_layer = 0

            key = regex_layer.sub('',params)
            if key not in contrib:
                try:
                    coords = re.match(regex_coords,params).group(1)
                    blank = dirpath.split('exec')[0]+"/obs_data/%s/blank.hdf5" % coords
                except AttributeError:
                    # No match, only 1 coord
                    blank = glob(dirpath.split('exec')[0]+"/obs_data/*/blank.hdf5")[0]

                contrib[key] = copy(MSDOpen(blank))

            pix_size = ( contrib[key].pixel_size(n_layer) / 1000. ) ** 2 # in km^2
            if oname == basename + ".out":
                pcl_name = [s for s in filenames if "pcl.bin" in s][0]
            else:
                pcl_name = '_'.join([ basename,'pcl',params+".bin" ])
            pcl_path = os.path.join(dirpath,pcl_name)
            pcl_data = load_bin(pcl_path)
            pcl_data *= val / pcl_data.sum()
            b = (pcl_data.shape[0] - contrib[key][n_layer].shape[0]) // 2
            contrib[key][n_layer] = pcl_data[b:-b,b:-b] if b else pcl_data

if p.full:
    results_names = ["Case"] + [ s[:s.index('(')] for s in lines[idx_results::2] ]
    print('\t'.join(results_names))
    for key,vals in outputs.items():
        print(key,*vals,sep="\t")
        if p.contrib:
            contrib[key].save(key)
else:
    for key,val in skyglow.items():
        print(key,val,sep="\t")
        if p.contrib:
            contrib[key].save(key)
