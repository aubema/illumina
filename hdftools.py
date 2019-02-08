#!/usr/bin/env python
#
# hdf5 tools for Illumina
#
# Author : Alexandre Simoneau
# unless noted otherwise
#
# January 2019

import MultiScaleData as MSD
import numpy as _np
import matplotlib.pyplot as _plt
from yaml import load as _yaml_load

def OpenCached(filename,cached={}):
    if filename in cached:
        return cached[filename]
    ds = MSD.Open(filename)
    cached[filename] = ds
    return ds

def plot(ds,n_layer=None,**options):
    _plt.gca().set_aspect(1)
    size = ds[0].shape[0]
    ind = _np.arange(-size/2,size/2+1)+0.5
    X,Y = _np.meshgrid(ind,ind)

    for i,layer in reversed(list(enumerate(ds[:n_layer]))):
        psize = ds.pixel_size(i)/1000.
        _plt.pcolor(X*psize,Y*psize,layer[::-1],**options)

def from_domain(params,data=None):
    if isinstance(params,str):
        with open(params) as f:
            params = _yaml_load(f)
    attrs = { k:v for k,v in params.iteritems() \
        if k not in ['extents','observers'] }
    attrs['obs_lat'] = [ d['latitude'] for d in params['observers'] ]
    attrs['obs_lon'] = [ d['longitude'] for d in params['observers'] ]
    attrs['layers'] = [ {k:v for k,v in d.iteritems() if k != 'layer'} \
        for d in params['extents'] ]
    return MSD.MultiScaleData(attrs,data)
