#!/usr/bin/env python
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

class MultiScaleData(_np.ndarray):
    def __new__(cls, data):
        obj = _np.array([ data[ly][:] for ly in data ]).view(cls)
        obj._attrs = dict(data.attrs)
        obj._attrs['layers'] = [ dict(data[ly].attrs) for ly in data ]

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

    def _get_col_row(self,coords,layer,asfloat=False):
        x,y = self._project(coords)
        attrs = self._attrs['layers'][layer]
        col = (x-attrs['xmin'])/attrs['pixel_size'] - 0.5
        row = (attrs['ymax']-y)/attrs['pixel_size'] - 0.5
        if not asfloat:
            col,row = int(round(col)),int(round(row))
        return col,row

    def _view_latlon(self,coords):
        n_layer = self._get_layer(coords)
        col,row = self._get_col_row(coords,n_layer)
        return data[col:col+1,row:row+1]

    def pixel_size(self,index):
        if isinstance(index, _Integral):
            n_layer = index
        else:
            n_layer = self._get_layer(index)
        return self._attrs['layers'][n_layer]['pixel_size']

    def get_obs_pos(self):
        return self._attrs['obs_lat'], self._attrs['obs_lon']

    def set_circle(self,center,radii,value):
        """Set a circle to a constant value on all layers.

        center: (lat,lon) tuple
        radii: Radius in meters
        value: Value to set the circle to"""
        H = self._attrs['nb_pixels']
        for i in range(len(self)):
            X0, Y0 = self._get_col_row(center,i,asfloat=True)
            R = float(radii) / self.pixel_size(i)
            Y, X = _np.ogrid[:H,:H]
            d2 = (X-X0)**2 + (Y-Y0)**2
            self[i][d2 <= R**2] = value

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

    def copy(self):
        return _clone(self)

def Open(filename):
    """Open a multiscale HDF5 data fileself.

    Returns a MultiScaleData object."""
    return MultiScaleData(_HDFile(filename))
