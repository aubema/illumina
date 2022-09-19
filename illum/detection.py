#!/usr/bin/env python3

import astropy.convolution as apconv
import numpy as np
import rasterio as rio
import scipy
import yaml

import illum.utils as u


def convolve(im, window):
    return apconv.convolve(im, apconv.Box2DKernel(width=window))


def threshold_detection(image, window):
    im = np.nan_to_num(image)
    local_emission = scipy.ndimage.median_filter(im, window)
    return im > local_emission


def rings(image, num_rings):
    mean = [image] + [
        convolve(np.nan_to_num(image), 2 * n + 3) for n in range(num_rings)
    ]
    rings = [image] + [
        ((2 * i + 1) ** 2 * mean[i] - (2 * i - 1) ** 2 * mean[i - 1]) / (8 * i)
        for i in range(1, len(mean))
    ]
    return (image > 0) & np.prod(
        [rings[i - 1] > rings[i] for i in range(1, len(rings))], axis=0
    )


def gaussian(image, sigma, window):
    laplacian = scipy.ndimage.gaussian_laplace(np.nan_to_num(image), sigma)
    minima = scipy.ndimage.minimum_filter(laplacian, size=(window, window))
    return (laplacian - minima) == 0


if __name__ == "__main__":
    with open("iss_params.in") as f:
        p = yaml.safe_load(f)

    rst = rio.open(f"{p['wd']}/Vrad.tiff")
    im = rst.read(1)

    im_threshold = threshold_detection(im, 3)
    im_rings = rings(im, 1)
    im_gaus = gaussian(im, 1, 3)
    im_lamps = im_threshold * im_rings * im_gaus
    u.save_geotiff(f"{p['wd']}/lamps.tiff", im_lamps, rst)
