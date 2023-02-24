import numpy as np
import pandas as pd
import rasterio as rio
import skimage
from scipy import signal


def gaussian_kernel(n, std, normalised=False):
    """
    Generates a n x n matrix with a centered gaussian
    of standard deviation std centered on it. If normalised,
    its volume equals 1."""
    gaussian1D = signal.gaussian(n, std)
    gaussian2D = np.outer(gaussian1D, gaussian1D)
    if normalised:
        gaussian2D /= 2 * np.pi * (std**2)
    return gaussian2D


rst = rio.open("Vrad.tiff")
img = np.nan_to_num(rst.read()[0])

psf = gaussian_kernel(13, 1, normalised=True)

deconv = skimage.restoration.richardson_lucy(
    img, psf, 30, clip=False, filter_epsilon=0.0001
)

# Alternate method to explicitly find the peaks

# peaks = skimage.feature.peak_local_max(deconv)

# mask = np.zeros_like(img, dtype=int)
# mask[tuple(peaks.T)] = np.arange(len(peaks))

# ws = skimage.segmentation.watershed(-deconv, mask)

ws = skimage.segmentation.watershed(-deconv)

df = pd.DataFrame(
    skimage.measure.regionprops_table(
        ws, deconv, properties=["centroid_weighted", "image_intensity"]
    )
)
df["image_intensity-sum"] = df["image_intensity"].map(np.sum)
df["longitude"], df["latitude"] = rio.transform.xy(
    rst.transform, df["centroid_weighted-0"], df["centroid_weighted-1"]
)

df.to_csv(
    "peaks.csv",
    index=False,
    columns=[
        "longitude",
        "latitude",
        "image_intensity-sum",
    ],
    header=["lon", "lat", "flux"],
)
