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

class MultiScaleData:
    def __init__(self, data):
        self._attrs = dict(data.attrs)
        self._levels = list()
        for i in range(data.attrs['nb_layers'],0,-1):
            self._levels.append((data['level_%i'%i][:],
                                 dict(data['level_%i'%i].attrs)))

    def __len__(self):
        return self._attrs['nb_layers']

    def _project(self,coords):
        lat,lon = coords
        wgs84 = _pyproj.Proj(init="epsg:4326")
        proj  = _pyproj.Proj(init=self._attrs['srs'])
        x,y   = _pyproj.transform(wgs84,proj,lon,lat)
        return x,y

    def _get_level(self,coords):
        x,y = self._project(coords)
        for i,level in enumerate(self._levels):
            data,attrs = level
            if attrs['xmin'] < x < attrs['xmax'] and \
               attrs['ymin'] < y < attrs['ymax']:
                return i
        raise IndexError('Coordinate out of range')

    def _get_col_row(self,coords,level,dtype=int):
        x,y = self._project(coords)
        data,attrs = self._levels[level]
        col = dtype((x-attrs['xmin'])/attrs['pixel_size'] - 0.5)
        row = dtype((y-attrs['ymin'])/attrs['pixel_size'] - 0.5)
        return col,row

    def _view_latlon(self,coords):
        n_level = self._get_level(coords)
        col,row = self._get_col_row(coords,n_level)
        return data[col:col+1,row:row+1]

    def __getitem__(self, index):
        if isinstance(index, _Integral):
            return self._levels[index][0][:]
        else:
            return self._view_latlon(index)[0,0]

    def __setitem__(self, index, value):
        if isinstance(index, _Integral):
            self._levels[index][0][:] = value
        else:
            self._view_latlon(index)[:] = value

    def pixel_size(self,index):
        if isinstance(index, _Integral):
            n_level = index
        else:
            n_level = self._get_level(index)
        return self._levels[n_level][1]['pixel_size']

    def get_obs_pos(self):
        return self._attrs['obs_lat'], self._attrs['obs_lon']

    def set_circle(self,center,radii,value):
        H = self._attrs['scale_factor']**2
        for i in range(len(self)):
            X0, Y0 = self._get_col_row(center,i,float)
            R = float(radii) / self.pixel_size(i)
            Y, X = _np.ogrid[:H,:H]
            d2 = (X-X0)**2 + (Y-Y0)**2
            self[i][d2 <= R**2] = value

    def fill(self,value):
        for i in range(len(self)):
            self[i] = value

    def clear(self):
        self.fill(0)

    def save(self,filename):
        if '.' not in filename or 'hdf' not in filename.rsplit('.',1)[-1]:
            filename += ".hdf5"
        with _HDFile(filename,'w') as File:
            for key,val in self._attrs.iteritems():
                File.attrs[key] = val
            for i in range(len(self)):
                ds = File.create_dataset("level_%d" % (len(self)-i), data=self[i])
                for key,val in self._levels[i][1].iteritems():
                    ds.attrs[key] = val

    def copy(self):
        return _clone(self)

def Open(filename):
    """Open a multiscale HDF5 data fileself.

    Returns a MultiScaleData object."""
    return MultiScaleData(_HDFile(filename))
