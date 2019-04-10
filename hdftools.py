#!/usr/bin/env python2
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
import matplotlib.colors as _colors
import yaml as _yaml

def OpenCached(filename,cached={}):
    if filename in cached:
        return cached[filename]
    ds = MSD.Open(filename)
    cached[filename] = ds
    return ds

def plot(ds,n_layer=None,log=False,**options):
    _plt.gca().set_aspect(1)
    n = ds[0].shape[0]
    buff = ds._attrs['buffer']
    N = n/2 - buff
    ind = _np.arange(-N-1,N+1)+0.5
    X,Y = _np.meshgrid(ind,ind)

    if not options.has_key('vmin'):
        options['vmin'] = ds.min()
    if not options.has_key('vmax'):
        options['vmax'] = ds.max()

    if log:
        options['norm'] = _colors.LogNorm(
            vmin=options['vmin'],
            vmax=options['vmax']
        )

    for i,layer in reversed(list(enumerate(ds[:n_layer]))):
        psize = ds.pixel_size(i)/1000.
        _plt.pcolor(
            X*psize,
            Y*psize,
            layer[buff:-buff,buff:-buff][::-1],
            **options
        )

def from_domain(params,data=None):
    if isinstance(params,str):
        with open(params) as f:
            params = _yaml.load(f,Loader=yaml.BaseLoader)
    attrs = { k:v for k,v in params.iteritems() \
        if k not in ['extents','observers'] }
    attrs['obs_lat'] = [ d['latitude'] for d in params['observers'] ]
    attrs['obs_lon'] = [ d['longitude'] for d in params['observers'] ]
    attrs['layers'] = [ {k:v for k,v in d.iteritems() if k != 'layer'} \
        for d in params['extents'] ]
    return MSD.MultiScaleData(attrs,data)
