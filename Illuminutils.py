#!/usr/bin/env python

from glob import glob
import gdal, yaml, h5py

def warp(srcfiles, projection=None, extent=None):
    bounding_box = [extent["xmin"],extent["ymin"],extent["xmax"],extent["ymax"]]

    vrt = gdal.BuildVRT('',srcfiles)
    ds = gdal.Warp('', vrt, format="VRT", dstSRS=projection,
                   outputBounds=bounding_box, xRes=extent["pixel_size"],
                   yRes=extent["pixel_size"], resampleAlg="cubicspline")

    return ds.GetRasterBand(1).ReadAsArray()

def get_MYD09A1_band_name(fname, band_n):
    sub_ds = gdal.Open(fname).GetSubDatasets()
    return filter(lambda sds: ("b%02d" % band_n) in sds[0], sub_ds)[0][0]

class Illuminutils:
    def __init__(self, ini):
        with open(ini) as f:
            self.params = yaml.load(f)

    def save(self, data, dstname, scale_factor=1.):
        with h5py.File(dstname+".hdf5",'w') as f:
            f.attrs['nb_layers'] = self.params['nb_layers']
            f.attrs['nb_pixels'] = self.params['nb_pixels']
            f.attrs['scale_factor'] = self.params['scale_factor']
            f.attrs['scale_min'] = self.params['scale_min']
            f.attrs['srs'] = self.params['srs']
            f.attrs['obs_lat'] = self.params['observer']['latitude']
            f.attrs['obs_lon'] = self.params['observer']['longitude']
            f.attrs['obs_x'] = self.params['observer']['x']
            f.attrs['obs_y'] = self.params['observer']['y']

            for i,layer in enumerate(data):
                P = self.params['extents'][i]
                ds = f.create_dataset("level_%d" % P['level'], data=layer*scale_factor)
                for key in P:
                    if key == "level":
                        continue
                    ds.attrs[key] = P[key]

    def srtm(self, folder):
        files = glob(folder+"/*.hgt")
        data = list()
        for i in range(self.params["nb_layers"]):
            data.append(warp(files, self.params['srs'], self.params["extents"][i]))
        self.save(data,"srtm")

    def modis(self, folder, band_n):
        files = glob(folder+"/*.hdf")
        fname = "refl_b%02d" % band_n
        band_names = map(lambda f: get_MYD09A1_band_name(f, band_n), files)
        data = list()
        for i in range(self.params["nb_layers"]):
            data.append(warp(band_names, self.params['srs'], self.params["extents"][i]))
        self.save(data, fname, scale_factor=0.0001)

    def viirs(self, folder):
        files = glob(folder+"/*.tif")
        data = list()
        for i in range(self.params["nb_layers"]):
            data.append(warp(files, self.params['srs'], self.params["extents"][i]))
        self.save(data,"stable_lights")

if __name__ == "__main__":
    I = Illuminutils(glob("*.ini")[0])

    I.srtm("SRTM")
    I.viirs("VIIRS-DNB")
    for i in [1,2,3,4]:
        I.modis("MODIS",i)
