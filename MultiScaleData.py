#!/usr/bin/env python3
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
from fractions import Fraction as _frac

class MultiScaleData:
    def __init__(self, params, data=None):
        if data is None:
            n_layers = params['nb_layers']
            data = [ _np.zeros((
                2*(params['nb_pixels'] + p['buffer']) + p['observer_size_y'],
                2*(params['nb_pixels'] + p['buffer']) + p['observer_size_x']
            )) for p in params['layers'] ]
        self._data = data
        self._attrs = _clone(params)

    def __len__(self):
        return len(self._data)

    def __iter__(self):
        return self._data.__iter__()

    def __getitem__(self,idx):
        return self._data[idx]

    def __repr__(self):
        p  = f"nb_layers: {self._attrs['nb_layers']}, "
        p += f"nb_pixels: {self._attrs['nb_pixels']*2+1}, "
        p += f"scale_min: {self._attrs['scale_min']:d}, "
        sf = self.scale_factor()
        p += f"scale_factor: {sf} ({float(sf):.3f}), "
        p += f"srs: {self._attrs['srs']}"
        return f"MultiScaleData{{{p}}}"

    def __setitem__(self,idx,value):
        self._data[idx] = value

    def _project(self,coords):
        lat,lon = coords
        wgs84 = _pyproj.Proj("epsg:4326")
        proj = _pyproj.Proj(self._attrs['srs'])
        x,y   = _pyproj.transform(wgs84,proj,lon,lat,always_xy=True)
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
        return self[n_layer][row:row+1,col:col+1]

    def at(self,lat,lon):
        return self._view_latlon((lat,lon))[0,0]

    def scale_factor(self):
        n_pix = 2*self._attrs['nb_pixels']+1
        n_core = 2*self._attrs['nb_core']+1
        return _frac(n_pix,n_core)

    def pixel_size(self,index):
        if isinstance(index, _Integral):
            n_layer = index
        else:
            n_layer = self._get_layer(index)
        return self._attrs['layers'][n_layer]['pixel_size']

    def copy(self):
        return _clone(self)

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
        for i in range(len(self)):
            ny,nx = self[i].shape
            X0, Y0 = self._get_col_row(center,i,asfloat=True)
            R = float(radii) / self.pixel_size(i)
            Y, X = _np.ogrid[:ny,:nx]
            d2 = (X-X0)**2 + (Y-Y0)**2
            self[i][d2 <= R**2] = value

    def set_overlap(self,value=0):
        nb_core = self._attrs['nb_core']
        for i in range(1,len(self)):
            ny,nx = self[i].shape
            obs_x = self._attrs['layers'][i]['observer_size_x']
            obs_y = self._attrs['layers'][i]['observer_size_y']
            self[i][
                (ny-obs_y)//2 - nb_core : (ny+obs_y)//2 + nb_core,
                (nx-obs_x)//2 - nb_core : (nx+obs_x)//2 + nb_core
            ] = value

    def set_buffer(self,value=0):
        for i in range(len(self)):
            buff = self._attrs['layers'][i]['buffer']
            ny,nx = self[i].shape
            self[i][:buff] = value
            self[i][ny-buff:] = value
            self[i][:,:buff] = value
            self[i][:,nx-buff:] = value

    def split_observers(self):
        for obs_id in range(len(self.get_obs_pos()[0])):
            yield self.extract_observer(obs_id)

    def extract_observer(self,obs_id):
        n = self._attrs['nb_pixels']

        x = self._attrs['obs_x'][obs_id]
        y = self._attrs['obs_y'][obs_id]
        lat = self._attrs['obs_lat'][obs_id]
        lon = self._attrs['obs_lon'][obs_id]

        new = MultiScaleData(self._attrs)
        for l in range(len(new)):
            b = self._attrs['layers'][l]['buffer']
            xc,yc = new._get_col_row((lat,lon),l)
            new[l] = self[l][yc-n-b:yc+n+b+1,xc-n-b:xc+n+b+1].copy()

            pix_size = new._attrs['layers'][l]['pixel_size']
            ny,nx = self[l].shape
            new._attrs['layers'][l]['xmin'] += (xc-n-b)*pix_size
            new._attrs['layers'][l]['xmax'] -= (nx-(xc+n+b+1))*pix_size
            new._attrs['layers'][l]['ymin'] += (yc-n-b)*pix_size
            new._attrs['layers'][l]['ymax'] -= (ny-(yc+n+b+1))*pix_size
            new._attrs['layers'][l]['observer_size_x'] = 1
            new._attrs['layers'][l]['observer_size_y'] = 1

        new._attrs['obs_x'] = [x]
        new._attrs['obs_y'] = [y]
        new._attrs['obs_lat'] = [lat]
        new._attrs['obs_lon'] = [lon]
        return new

    def save(self,filename):
        if '.' not in filename or \
            'hdf' not in filename.rsplit('.',1)[-1]:
            filename += ".hdf5"
        with _HDFile(filename,'w') as File:
            for key,val in list(self._attrs.items()):
                if key != "layers":
                    if "obs" not in key:
                        File.attrs[key] = val
                    else:
                        ds = File.create_dataset(key.replace('_','/'), data=val)
            for i in range(len(self)):
                ds = File.create_dataset("layers/%d" % i, data=self[i])
                for key,val in list(self._attrs['layers'][i].items()):
                    ds.attrs[key] = val

def Open(filename):
    """Open a multiscale HDF5 data fileself.

    Returns a MultiScaleData object."""
    ds = _HDFile(filename,'r')
    data = [ ds['layers'][n][:] for n in sorted(ds['layers'],key=int) ]
    params = dict(ds.attrs)
    params['layers'] = [ dict(ds['layers'][n].attrs) for n in sorted(ds['layers'],key=int) ]
    params.update( ('obs_'+k, ds['obs'][k][:]) for k in ds['obs'] )
    try:
        params['srs'] = params['srs'].decode("utf-8")
    except AttributeError:
        pass
    return MultiScaleData(params, data)
