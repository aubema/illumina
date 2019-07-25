#!/usr/bin/env python2
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
    metavar=('NAME','VALUE(S)'), help="Values of the parameter NAME to extract. "\
        "Multiple values must be separated by commas.")

p = parser.parse_args()
regex = re.compile(r'-layer_(\d+)')

skyglow = ddict(float)
if p.params_filename is not None:
    contrib = ddict(partial(from_domain,p.params_filename))

for dirpath,dirnames,filenames in os.walk(p.exec_dir):
    if not os.path.isfile(os.path.join(dirpath,'illumina.in')):
        continue

    with open(os.path.join(dirpath,'illumina.in')) as f:
        basename = f.readlines()[1].split()[0]

    out_names = filter(lambda fname: fname.endswith(".out") and \
        fname.startswith(basename+'_'), filenames)
    if not out_names:
        out_names = [ basename + '.out' ]

    for oname in out_names:
        if oname == basename + ".out":
            params = dirpath.split('exec'+os.sep)[1].replace(os.sep,'-')
        else:
            params = oname.rstrip('.out').strip(basename+'_')

        try:
            for pname,pvals in p.param:
                if pname not in params:
                    print "ERROR: Parameter '%s' not found." % pname
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

        val = float(lines[-4])
        skyglow[regex.sub('',params)] += val
        if p.params_filename is not None and val != 0:
            n_layer = int(regex.search(params).groups()[0])
            key = regex.sub('',params)
            pix_size = ( contrib[key].pixel_size(n_layer) / 1000. ) ** 2 # in km^2
            if oname == basename + ".out":
                pcl_name = filter(lambda s: "pcl.bin" in s, filenames)[0]
            else:
                pcl_name = '_'.join([ basename,'pcl',params+".bin" ])
            pcl_path = os.path.join(dirpath,pcl_name)
            pcl_data = load_bin(pcl_path)
            pcl_data *= val / pix_size / pcl_data.sum()
            contrib[key][n_layer] = pcl_data

for key,val in skyglow.iteritems():
    print key,val
    if p.params_filename is not None:
        contrib[key].save(key)
