import errno
import os
from glob import glob

import astropy.io.fits as pyfits
import numpy as np
import yaml
from astropy.convolution import Box2DKernel, convolve
from osgeo import gdal, osr
from scipy import stats


def fits2tiff(filename):
    hdulist = pyfits.open(filename)
    other = np.asarray(hdulist[0].data)
    earth = other.copy()
    earth = np.float32(earth)
    nrows, ncols = earth.shape[0], earth.shape[1]
    filename = filename[:-5] + ".tiff"
    dst_ds = gdal.GetDriverByName("GTiff").Create(
        filename, ncols, nrows, 1, gdal.GDT_Float32
    )
    sr = osr.SpatialReference()
    sr.ImportFromEPSG(4326)
    dst_ds.SetProjection(sr.ExportToWkt())
    dst_ds.SetGeoTransform(
        (
            hdulist[0].header["CRVAL1"],
            hdulist[0].header["CRPIX1"],
            0,
            hdulist[0].header["CRVAL2"],
            0,
            hdulist[0].header["CRPIX2"],
        )
    )
    dst_ds.GetRasterBand(1).WriteArray(earth)
    return filename


def open_tiff(filename, dtype=np.float32):
    if filename.endswith("fits"):
        filename = fits2tiff(filename)

    # Load file, and access the band and get a NumPy array
    if not os.path.isfile(filename):
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT), filename
        )
    src = gdal.Open(filename, gdal.GA_Update)
    band = src.GetRasterBand(1)
    ar = band.ReadAsArray()
    return src, ar


with open("iss_params.in") as f:
    p = yaml.safe_load(f)

exp_folds = glob("iss???e*/")
if not exp_folds:
    print("ERROR: Data folder not found. Please make sure it was not renamed.")
    exit()
if len(exp_folds) > 1:
    print(f"ERROR: Too many data folders found. {exp_folds}")
    exit()

exp_fold = exp_folds[0]
basename, ts = exp_fold.split("_")
print(f"Data `{basename}` ({ts[:4]}-{ts[4:6]}-{ts[6:8]}) found.")

# open intensity and technology (spectral class) images
src, image_intensity = open_tiff(
    f"{exp_fold}/popular/Corr_{basename}_ImpactVlG_GR.fits"
)
src, image_tech = open_tiff(
    f"{exp_fold}/popular/Corr_{basename}CompositeW.fits"
)

# Defining domain

src, domain = open_tiff(f"{exp_fold}/quality/{basename}RGBdisto2_rect.tiff")
domain[domain != 0] = 1

# 0. Find saturated pixels that have no value in the intensity image and
#    replace the nan by the maximum of the image
sat = os.sys.float_info.max
novalue = -1e30

# 1. Start of treatment : elimination of noise
# 1a. Eliminate negative data, mistake resulting of the pre-treatment
image_intensity[image_intensity < 0] = np.nan


# 1c. We create a temporary map with the pixels of less value since the noise
#     should be smaller or around the same value as the less valuable pixels.
#     We then extract statistical data of these small values

image_intensity -= np.nanmax(image_intensity) * p["threshold"]

# 1e. Eliminating negative pixels created by our treatment of the noise
#     in the image
image_intensity[image_intensity < 0] = np.nan

# 1f. Compare your image to the ones in ReadMe for reference as well as Google
#     maps to know which zones should emit light or not. We want to eliminate a
#     certain quantity of small value pixels in aeras that shouldn't emit light
#     At the same time, we want to limit the values eliminated in aeras where
#     light should be emitted

# 1.5. Finding unexplicable (so far) void pixels surrounded by high intensity
#      pixels and filling them with the mean of the pixels around
#      creating binary images in intensity, value=1, nan=0


def image_to_binary(image):
    im = image.copy()
    im[im >= 0] = 1
    im[np.isnan(image)] = 0
    return im


def binary_mode_classes(im_tech, window):
    tech_conv = [
        convolve(im_tech == i, Box2DKernel(width=window)) for i in range(1, 7)
    ]
    im_t = np.argmax(tech_conv, axis=0) + 1
    im_t[~np.any(tech_conv, axis=0)] = 0
    return im_t


def convolution_nb_void(image, im_tech, window, keep_value):
    im = image.copy()
    nb_nan_binary = convolve(image_to_binary(image), Box2DKernel(width=window))
    nb_nan_real = convolve(np.nan_to_num(image), Box2DKernel(width=window))
    mean = np.nanmean(image)
    nb_nan_binary_without0 = nb_nan_binary.copy()
    nb_nan_binary_without0[nb_nan_binary_without0 == 0] = 1
    mask = (
        (np.nan_to_num(image) == 0)
        & (nb_nan_binary > keep_value / window ** 2)
        & ((nb_nan_real / nb_nan_binary_without0) > mean)
    )
    tech = im_tech.copy()
    im_t = binary_mode_classes(im_tech, window)
    im[mask] = nb_nan_real[mask] / nb_nan_binary[mask]
    tech[mask] = im_t[mask]

    return im, tech, np.sum(mask)


image_tech[image_tech == 0] = np.nan
n_changed = 1
while n_changed:
    image_intensity, image_tech, n_changed = convolution_nb_void(
        image_intensity, image_tech, window=3, keep_value=3
    )
    print("Void pixel filled:", n_changed)


# 2. Elimination of remaining values in dark aeras
#    If your image is already free of valued pixels in dark aeras following
#    step 1, comment this section and skip to step 3

# 2a. We define our convolution fonction
def convolution_nb_nan(image, window, keep_value):
    im = image.copy()
    nb_nan = convolve(image_to_binary(image), Box2DKernel(width=window))
    mask = (np.nan_to_num(im) != 0) & ((window ** 2 * nb_nan) < keep_value)
    im[mask] = np.nan
    return im, np.sum(mask)  # np.nan_to_num(im)


# 2b. We create a binary image where value pixels are equal to 1 and NaN pixels
#     are equal to 0
n_changed = 1
while n_changed:
    image_intensity, n_changed = convolution_nb_nan(
        image_intensity, window=5, keep_value=3
    )
    print("Isolated pixel removed:", n_changed)


# 3. Concordance between intensity and technology images


def int_tech_comparison(intensity, im_tech):
    tech = im_tech.copy()
    im_t = binary_mode_classes(im_tech, window=3)
    tech[((np.nan_to_num(intensity) == 0) & (np.nan_to_num(im_tech) != 0))] = 0
    mask = (np.nan_to_num(intensity) != 0) & (np.nan_to_num(im_tech) == 0)
    tech[mask] = im_t[mask]
    return tech


image_tech = int_tech_comparison(image_intensity, image_tech)
image_tech[image_tech == 0] = np.nan

# 7. Saving results


def save_geotiff(filename, data):
    nband = 1
    nrow, ncol = data.shape
    driver = gdal.GetDriverByName("GTiff")
    dst_dataset = driver.Create(
        filename + ".tiff", ncol, nrow, nband, gdal.GDT_Float32
    )
    # sets same geotransform as input
    dst_dataset.SetGeoTransform(src.GetGeoTransform())
    # sets same projection as input
    dst_dataset.SetProjection(src.GetProjection())
    dst_dataset.GetRasterBand(1).WriteArray(data.astype(float))
    dst_dataset = None


if not os.path.isdir(p["wd"]):
    os.makedirs(p["wd"])

save_geotiff(p["wd"] + "/Vrad", image_intensity)
save_geotiff(p["wd"] + "/tech", image_tech)
save_geotiff(p["wd"] + "/domain", domain)
