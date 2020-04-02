#!/usr/bin/env python3
#
# Illumina output extract
#
# Author : Alexandre Simoneau
# unless noted otherwise
#
# February 2019

import click
import os, re
from collections import defaultdict as ddict
import MultiScaleData as MSD
from functools import partial
from pytools import load_bin
from glob import glob
from copy import deepcopy as copy

def MSDOpen(filename,cached={}):
    if filename in cached:
        return cached[filename]
    ds = MSD.Open(filename)
    cached[filename] = ds
    return ds

@click.command()
@click.argument("exec_dir",default='.',type=click.Path(exists=True))
@click.option("-c","--contrib",is_flag=True,help="If present, extract contribution maps.")
@click.option("-p","--params",multiple=True,nargs=2,
    help="Parameter name,value pair to extract. Can be provided more than once.")
def extract(exec_dir,contrib,params):
    """Extract Illumina outputs.

    Will walk the EXEC_DIR to locate and extract illumina outputs.
    May fail if some execution failed, so one should validate the completude
    of the runs with 'illum failed' before using this.

    If not given, EXEC_DIR will default to the current directory.
    """
    regex_layer = re.compile(r'-layer_(\d+)')
    regex_coords = re.compile(r'observer_coordinates_(-?\d+\.\d+_-?\d+\.\d+)')

    skyglow = ddict(float)
    if contrib:
        contributions = dict()

    for dirpath,dirnames,filenames in os.walk(exec_dir):
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
                params_name = dirpath.split('exec'+os.sep)[1].replace(os.sep,'-')
            else:
                params_name = oname[len(basename)+1:-4]

            try:
                for pname,pvals in params:
                    if pname not in params_name:
                        print("ERROR: Parameter '%s' not found." % pname)
                        exit()
                    for pval in pvals.split(','):
                        if "%s_%s" % (pname,pval) in params_name:
                            break
                    else:
                        raise ValueError()
            except ValueError:
                continue

            with open(os.sep.join([dirpath,oname])) as f:
                lines = f.readlines()

            val = float(lines[-1])
            skyglow[regex_layer.sub('',params_name)] += val
            if contrib:
                try:
                    n_layer = int(regex_layer.search(params_name).groups()[0])
                except AttributeError:
                    # No match, only 1 layer
                    n_layer = 0

                key = regex_layer.sub('',params_name)
                if key not in contributions:
                    try:
                        coords = re.match(regex_coords,params_name).group(1)
                        blank = dirpath.split('exec')[0]+"/obs_data/%s/blank.hdf5" % coords
                    except AttributeError:
                        # No match, only 1 coord
                        blank = glob(dirpath.split('exec')[0]+"/obs_data/*/blank.hdf5")[0]

                    contributions[key] = copy(MSDOpen(blank))

                pix_size = ( contributions[key].pixel_size(n_layer) / 1000. ) ** 2 # in km^2
                if oname == basename + ".out":
                    pcl_name = [s for s in filenames if "pcl.bin" in s][0]
                else:
                    pcl_name = '_'.join([ basename,'pcl',params_name+".bin" ])
                pcl_path = os.path.join(dirpath,pcl_name)
                pcl_data = load_bin(pcl_path)
                pcl_data *= val / pcl_data.sum()
                b = (pcl_data.shape[0] - contributions[key][n_layer].shape[0]) // 2
                contributions[key][n_layer] = pcl_data[b:-b,b:-b] if b else pcl_data

    for key,val in skyglow.items():
        print(key,val)
        if contrib:
            contributions[key].save(key)
