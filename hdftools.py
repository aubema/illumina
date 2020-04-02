#!/usr/bin/env python3
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

def plot(ds,n_layer=None,log=False,area=False,**options):
    _plt.gca().set_aspect(1)

    for i,layer in reversed(list(enumerate(ds[:n_layer]))):
        layer = layer.copy()
        n = layer.shape[0]
        buff = ds._attrs['layers'][i]['buffer']

        psize = ds.pixel_size(i)/1000.
        if area:
            layer /= psize**2

        N = psize * (n/2 - buff)

        if log:
            norm = _np.array([ (ds.pixel_size(i)/1000)**2 for i in range(len(ds)) ])
        else:
            norm = 1.
        if 'vmin' not in options:
            options['vmin'] = min(_np.array([_np.min(l) for l in ds])/norm)
        if 'vmax' not in options:
            options['vmax'] = max(_np.array([_np.max(l) for l in ds])/norm)

        if log:
            options['norm'] = _colors.LogNorm(
                vmin=options['vmin'],
                vmax=options['vmax']
            )

        _plt.imshow(
            layer[buff:n-buff,buff:n-buff],
            extent=(-N,N,-N,N),
            **options
        )

    N = ds[-1].shape[0]/2 - ds._attrs['layers'][-1]['buffer']
    N *= ds.pixel_size(len(ds)-1)/1000
    _plt.xlim(-N,N)
    _plt.ylim(-N,N)

def from_domain(params,data=None):
    if isinstance(params,str):
        with open(params) as f:
            params = _yaml.safe_load(f)
    attrs = { k:v for k,v in list(params.items()) \
        if k not in ['extents','observers'] }
    attrs['obs_lat'] = [ d['latitude'] for d in params['observers'] ]
    attrs['obs_lon'] = [ d['longitude'] for d in params['observers'] ]
    attrs['obs_x'] = [ d['x'] for d in params['observers'] ]
    attrs['obs_y'] = [ d['y'] for d in params['observers'] ]
    attrs['layers'] = [ {k:v for k,v in list(d.items()) if k != 'layer'} \
        for d in params['extents'] ]
    return MSD.MultiScaleData(attrs,data)
