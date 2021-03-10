import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.utils.data import download_file
from scipy import stats
from astropy.convolution import convolve
from scipy.stats import norm
from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans, Box2DKernel
from osgeo import gdal

pix_size= 8 # taille du pixel en micron
focal_distance= 400 # distance en mm
ISS_altitude= 400000 # altutide de l'ISS en m
deviation= (25* (focal_distance/1000))/((pix_size/1000000)*ISS_altitude) # stdev to Blur final image with a 25m resolution
# on estime que la distance entre les lampadaires est typiquement de 50m dans une ville

file_name_intensity = 'Corr_iss035e017088_ImpactVlG_GR.fits' # nom de l'image d'intensité de la ville à étudier (ImpactVlG_GR)
file_name_technology = 'Corr_iss035e017088Composite.fits' # nom de l'image des technologies d'éclairement de la ville à étudier (Composite)



file_intensity = fits.open(file_name_intensity) 
file_tech = fits.open(file_name_technology)
# file_intensity.info()



image_intensity =file_intensity[0].data.astype(np.float32)
# print(type(image_intensity))
# print(image_intensity.shape)
image_tech = file_tech[0].data.astype(np.float32)
plt.figure()
plt.imshow(image_tech, cmap="rainbow")
plt.colorbar()
plt.title('')
plt.show()



file_intensity.close() #											pas pour tech?
ncol=np.size(image_tech,1) #										pour la sauvegarde à la fin
nrow=np.size(image_tech,0)





## 1. Start of treatment : elimination of noise

## 1a. We first eliminate negative data, as these are a mistake resulting of the pre-treatment
image_intensity [image_intensity<0] = np.nan

## 1b. Statistical data
mean=np.nanmean(image_intensity)
median=np.nanmedian(image_intensity)
standard_deviation=np.nanstd(image_intensity)
mode= stats.mode (image_intensity[~np.isnan(image_intensity)],axis=None)[0][0]
print("mode:", mode, "mean:", mean, "median:", median, "standard_deviation:", standard_deviation)


## 1c. We create a temporary map with the pixels of less value since the noise should be smaller or around the same value as the less valuable pixels
##     We then extract statistical data of these small values
small_values=image_intensity[image_intensity<standard_deviation]
mean2=np.nanmean(small_values)
median2=np.nanmedian(small_values)
standard_deviation2=np.nanstd(small_values)
mode2= stats.mode (small_values[~np.isnan(small_values)],axis=None)[0][0]
print("mean2:", mean2 , "median2:",  median2, "mode2:",  mode2, "standard_deviation2:", standard_deviation2)


## 1d. Using the statistical data on hand, we estimate the value of the noise. 
Noise= mode + standard_deviation2 ## this value can change. Modifie accordingly in order to achieve the best results (see step 1.f)
image_intensity-= Noise 

## 1e. Eliminating negative pixels created by our treatment of the noise in the image
image_intensity [image_intensity<0] = np.nan

## 1f. Compare your image to the ones in ReadMe for reference as well as Google maps to know which zones should emit light or not
##     We want to eliminate a certain quantity of small value pixels in aeras that shouldn't emit light
##     At the same time, we want to limit the values eliminated in aeras where light should be emitted
plt.figure()
plt.imshow(image_intensity, cmap="rainbow")
plt.colorbar()
plt.title('Intensity without noise')
# plt.show()


## 2. Elimination of remaining values in dark aeras
##    If your image is already free of valued pixels in dark aeras following step 1, comment this section and skip to step 3

## We need copies of the image we currently have. The number of copies needed is equal to the number of convolution exectuted
image_copy1=image_intensity.copy()
image_copy2=image_intensity.copy()


## 2a. We define our convolution fonction
def Convolution_nb_nan (image, width):
	image=image.copy()
	kernel=Box2DKernel(width=width)
	nb_nan=convolve(image, kernel)
	return nb_nan


## 2b. We set the values needed to do a first convolution. 
window=3 ## Size of the window, example 3X3, 5X5, 7X7, etc. 
keep_Value=4 ## Minimum of surrounding value pixels accepted in the window size in order to keep the convolued pixel
Threshold=(keep_Value)/window**2 

## 2c. We create a binary image where value pixels are equal to 1 and NaN pixels are equal to 0
image_intensity_NaN = np.isnan(image_intensity)
image_intensity [image_intensity>=0] = 1
image_intensity[image_intensity_NaN] = 0
image_binary = image_intensity

## 2d. First convolution
image_binary=Convolution_nb_nan(image_binary, window)
image_copy1[image_binary<Threshold]=np.nan
image_intensity=image_copy1 

## 2e. Second convolution if necessary. If not, skip to step 3

## 2f. We set the new values needed for the second convolution
window=5
keep_Value=12
Threshold=(keep_Value)/window**2

## 2g. We create a new binary image from the one we currently have (resulting from the first convolution)
image_intensity_NaN = np.isnan(image_intensity)
image_intensity [image_intensity>=0] = 1
image_intensity[image_intensity_NaN] = 0
image_binary = image_intensity

## 2h. Second convolution
image_binary=Convolution_nb_nan(image_binary, window)
image_copy2[image_binary<Threshold]=np.nan
image_intensity=image_copy2

## 2i. Compare your image to the examples in the ReadMe
plt.figure()
plt.imshow(image_intensity, cmap="rainbow")
plt.colorbar()
plt.title('Intensity with clean dark aeras')
# plt.show()



## 3. Concordance...

#3a. We make a copy of our image to create a new binary map without changing our current image
image_intensity_temporary=image_intensity.copy()
image_intensity_temporary_NaN = np.isnan(image_intensity_temporary)
image_intensity_temporary[image_intensity_temporary>=0] = 1
image_intensity_temporary[image_intensity_temporary_NaN] = 0
image_binary=np.nan_to_num(image_copy2)

# 3b. eliminating points of value=0
image_tech [image_tech==0] = np.nan

# 3c. By multiplying our technology map with the binary of intensity, we eliminate data where we have none in the intensity image
image_tech = image_tech * image_binaire2
image_tech[image_tech==0] = np.nan

plt.figure()
plt.imshow(image_tech2, cmap='rainbow')
plt.title("Technologies")
# plt.show()
image_tech_non_floue=image_tech2.copy()

# Copie de l'image des technologie initiale pour manipulations ultérieures
# image_tech_traite = image_tech2.copy() 
# plt.figure()
# plt.imshow(image_tech_final, cmap='rainbow')
# plt.colorbar()
# plt.title('Technologies, final')

# # Création image des technologies binaire
# image_tech_NaN = np.isnan(image_tech2)
# image_tech2[image_tech_NaN] = 0
# image_tech2[image_tech2>0] = 1
# image_binaire_tech=image_tech2

# # Concordance entre les deux images 
# image_intensity_finale = image_binaire_tech * image_intensity_treated
# plt.figure()
# plt.imshow(image_intensity_finale, cmap='rainbow')
# plt.colorbar()
# plt.title('Intensité lumineuse, final')

# Fonction de floutage

def Convolution_sans_zero(image,stdev=deviation):
	image=image.copy()
	mask = image==0
	image[mask]=np.nan
	kernel = Gaussian2DKernel(x_stddev=stdev)
	Blur = convolve(image, kernel)
	Blur[mask]=np.nan
	return Blur



# Floutage de l'image de l'intensité
Blur_intensity=Convolution_sans_zero(image_intensity_traite)
plt.figure()
plt.imshow(Blur_intensity, cmap= 'rainbow')
plt.colorbar()
plt.title('Intensité flouée')

# Variables des valeurs de MSI associées aux catégories de technologies
MSI_1=0.539
MSI_2=0.446
MSI_3=0.274
MSI_4=0.043
MSI_5=0.118
MSI_6=0.017

# Création carte MSI par remplacement des valeurs de 1 à 6 par les valeurs de MSI associées
image_MSI = np.zeros_like(image_tech2)
image_MSI[image_tech2 == 1] = MSI_1
image_MSI[image_tech2 == 2] = MSI_2
image_MSI[image_tech2 == 3] = MSI_3
image_MSI[image_tech2 == 4] = MSI_4
image_MSI[image_tech2 == 5] = MSI_5
image_MSI[image_tech2 == 6] = MSI_6
plt.figure()
plt.imshow(image_MSI, cmap='rainbow')
plt.colorbar()
plt.title('Carte MSI')

# Floutage de la carte des MSI
# stdev dépend de la résolution de pixels
Blur_MSI = Convolution_sans_zero(image_MSI)
plt.figure()
plt.imshow(Blur_MSI, cmap = 'rainbow')
plt.colorbar()
plt.title('MSI flouée')

# Création de la carte de l'impact final (MSI et intensité combinés)
Impact_MSI=Blur_MSI*Blur_intensity
plt.figure()
plt.imshow(Impact_MSI, cmap="rainbow")
plt.colorbar()
plt.title('Impact MSI')
# plt.show()

plt.figure()
plt.imshow(Blur_intensity, cmap="rainbow")
plt.colorbar()
plt.title("blur intensity")
plt.show()


 


final_data = Blur_intensity
# get parameters
# geotransform = src_dataset.GetGeoTransform()
# spatialreference = src_dataset.GetProjection()
nband = 1
print(ncol,nrow)
print(type(ncol))
# create dataset for output
fmt = 'GTiff'
driver = gdal.GetDriverByName(fmt)
dst_dataset = driver.Create('C:/Users/nikki/OneDrive/Bureau/Blur_intensity.tiff', ncol, nrow, nband, gdal.GDT_Float32)
# dst_dataset.SetGeoTransform(geotransform)
# dst_dataset.SetProjection(spatialreference)
dst_dataset.GetRasterBand(1).WriteArray(final_data.astype(float))
dst_dataset = None

np.save('Blur_Vrad',Blur_intensity) 




final_data = Blur_MSI
# get parameters
# geotransform = src_dataset.GetGeoTransform()
# spatialreference = src_dataset.GetProjection()
nband = 1
print(ncol,nrow)
print(type(ncol))
# create dataset for output
fmt = 'GTiff'
driver = gdal.GetDriverByName(fmt)
dst_dataset = driver.Create('C:/Users/nikki/OneDrive/Bureau/Impact_MSI.tiff', ncol, nrow, nband, gdal.GDT_Float32)
# dst_dataset.SetGeoTransform(geotransform)
# dst_dataset.SetProjection(spatialreference)
dst_dataset.GetRasterBand(1).WriteArray(final_data.astype(float))
dst_dataset = None

np.save('Impact_MSI',Impact_MSI)








final_data = Impact_MSI
# get parameters
# geotransform = src_dataset.GetGeoTransform()
# spatialreference = src_dataset.GetProjection()
nband = 1
print(ncol,nrow)
print(type(ncol))
# create dataset for output
fmt = 'GTiff'
driver = gdal.GetDriverByName(fmt)
dst_dataset = driver.Create('C:/Users/nikki/OneDrive/Bureau/Impact_MSI.tiff', ncol, nrow, nband, gdal.GDT_Float32)
# dst_dataset.SetGeoTransform(geotransform)
# dst_dataset.SetProjection(spatialreference)
dst_dataset.GetRasterBand(1).WriteArray(final_data.astype(float))
dst_dataset = None

np.save('Impact_MSI',Impact_MSI)