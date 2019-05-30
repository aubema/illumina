#!/usr/bin/env python2
#
# Multi-scale data handling
#
# Author : Alexandre Simoneau
# unless noted otherwise
#
# January 2019

import numpy as _np
import pyproj as _pyproj
from h5py import File as _HDFile
from numbers import Integral as _Integral
from copy import deepcopy as _clone
from itertools import izip as _izip

class MultiScaleData(_np.ndarray):
    def __new__(cls, params, data=None):
        if data is None:
            n_layers = params['nb_layers']
            nb_pixels = 2*(params['nb_pixels'] + params['buffer'])
            data = [ _np.zeros((
                nb_pixels + p['observer_size_y'],
                nb_pixels + p['observer_size_x']
            )) for p in params['layers'] ]
        obj = _np.asarray(data).view(cls)
        obj._attrs = _clone(params)

        return obj

    def __array_finalize__(self, obj):
        if obj is None: return
        self._attrs = getattr(obj, '_attrs', None)

    def _project(self,coords):
        lat,lon = coords
        wgs84 = _pyproj.Proj(init="epsg:4326")
        proj  = _pyproj.Proj(init=self._attrs['srs'])
        x,y   = _pyproj.transform(wgs84,proj,lon,lat)
        return x,y

    def _get_layer(self,coords):
        x,y = self._project(coords)
        for i,layer in enumerate(self._attrs['layers']):
            attrs = layer
            if attrs['xmin'] < x < attrs['xmax'] and \
               attrs['ymin'] < y < attrs['ymax']:
                return i
        raise IndexError('Coordinate out of range')

    def _get_col_row(self,coords,layer,asfloat=False,proj=False):
        if not proj:
            x,y = self._project(coords)
        else:
            x,y = coords
        attrs = self._attrs['layers'][layer]
        col = (x-attrs['xmin'])/attrs['pixel_size'] - 0.5
        row = (attrs['ymax']-y)/attrs['pixel_size'] - 0.5
        if not asfloat:
            try:
                col,row = int(round(col)),int(round(row))
            except TypeError:
                col = _np.round(col).astype(int)
                row = _np.round(row).astype(int)
        return col,row

    def _view_latlon(self,coords):
        n_layer = self._get_layer(coords)
        col,row = self._get_col_row(coords,n_layer)
        return data[row:row+1,col:col+1]

    def scale_factor(self):
        return (self._attrs['nb_pixels']+0.5) / (self._attrs['nb_core']+0.5)

    def pixel_size(self,index):
        if isinstance(index, _Integral):
            n_layer = index
        else:
            n_layer = self._get_layer(index)
        return self._attrs['layers'][n_layer]['pixel_size']

    def get_obs_pos(self,proj=False):
        if proj:
            return self._attrs['obs_x'],self._attrs['obs_y']
        else:
            return self._attrs['obs_lat'], self._attrs['obs_lon']

    def set_circle(self,center,radii,value):
        """Set a circle to a constant value on all layers.

        center: (lat,lon) tuple
        radii: Radius in meters
        value: Value to set the circle to"""
        for i in xrange(len(self)):
            ny,nx = self[i].shape
            X0, Y0 = self._get_col_row(center,i,asfloat=True)
            R = float(radii) / self.pixel_size(i)
            Y, X = _np.ogrid[:ny,:nx]
            d2 = (X-X0)**2 + (Y-Y0)**2
            self[i][d2 <= R**2] = value

    def set_overlap(self,value=0):
        nb_core = self._attrs['nb_core']
        for i in xrange(1,len(self)):
            ny,nx = self[i].shape
            obs_x = self._attrs['layers'][i]['observer_size_x']
            obs_y = self._attrs['layers'][i]['observer_size_y']
            self[i][
                (ny-obs_y)/2 - nb_core : (ny+obs_y)/2 + nb_core,
                (nx-obs_x)/2 - nb_core : (nx+obs_x)/2 + nb_core
            ] = value

    def set_buffer(self,value=0):
        buff = self._attrs['buffer']
        for i in xrange(len(self)):
            ny,nx = self[i].shape
            self[i][:buff] = value
            self[i][ny-buff:] = value
            self[i][:,:buff] = value
            self[i][:,nx-buff:] = value

    def split_observers(self):
        for coords in _izip( *self.get_obs_pos(proj=True) ):
            yield self.extract_observer(coords,proj=True)

    def extract_observer(self,coords,proj=False):
        n = self._attrs['nb_pixels']
        b = self._attrs['buffer']
        lat,lon = coords

        new = MultiScaleData(self._attrs)
        for l in xrange(len(new)):
            xc,yc = new._get_col_row(coords,l,proj=proj)
            new[l] = self[l][yc-n-b:yc+n+b+1,xc-n-b:xc+n+b+1].copy()
            new._attrs['layers'][l]['observer_size_x'] = 1
            new._attrs['layers'][l]['observer_size_y'] = 1
        new._attrs['obs_x'] = [xc]
        new._attrs['obs_y'] = [yc]
        return new

    def save(self,filename):
        if '.' not in filename or \
            'hdf' not in filename.rsplit('.',1)[-1]:
            filename += ".hdf5"
        with _HDFile(filename,'w') as File:
            for key,val in self._attrs.iteritems():
                if key != "layers":
                    File.attrs[key] = val
            for i in xrange(len(self)):
                ds = File.create_dataset("layer_%d" % i, data=self[i])
                for key,val in self._attrs['layers'][i].iteritems():
                    ds.attrs[key] = val

def Open(filename):
    """Open a multiscale HDF5 data fileself.

    Returns a MultiScaleData object."""
    ds = _HDFile(filename)
    data = _np.array([ ds[ly][:] for ly in ds ])
    params = dict(ds.attrs)
    params['layers'] = [ dict(ds[ly].attrs) for ly in ds ]
    return MultiScaleData(params, data)
