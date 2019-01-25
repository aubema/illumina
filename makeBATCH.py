#!/usr/bin/env python
#
# batch processing
#
# Author : Alexandre Simoneau
# unless noted otherwise
#
# January 2019

import MultiScaleData as MSD
import sys, yaml, os, shutil
from pytools import save_bin
from glob import glob
from itertools import product as comb
import numpy as np

def getp(name,name_list,param_list,params):
    return param_list[name_list.index(name)] if name in name_list else params[name]

# Load parameters
if len(sys.argv) < 2:
    print "Error: Must provide input parameter file location."
    exit()

with open(sys.argv[1]+"/inputs_params.in") as f:
    params = yaml.load(f)

if len(sys.argv) > 2:
    params['batch_file_name'] = sys.argv[2]

for fname in glob(sys.argv[1]+"/%s*" % params['batch_file_name']):
    os.remove(fname)

# Add wavelength and multiscale
params['wavelength'] = np.loadtxt("wav.lst").tolist()
params['layer'] = range(len(MSD.Open("stable_lights.hdf5")))

with open("zon.lst") as f:
    zones = f.read().split()

# Clear and create execution folder
dir_name = "exec/"
shutil.rmtree(dir_name,True)
os.makedirs(dir_name)

count = 0
multival = filter(lambda k: isinstance(params[k],list),params)
multival = sorted(multival,key=lambda s:len(s),reverse=True) # Semi-arbitrary sort
param_space = [params[k] for k in multival]
print "Number of executions:", np.prod(map(len,param_space))
for p in comb(*param_space):
    if "azimuth_angles" in multival:
        el = getp("elevation_angles",multival,p,params)
        if el == 90 and params['azimuth_angles'].index(p[multival.index("azimuth_angles")]) != 0:
            continue

    fold_name = dir_name + '/'.join([ multival[i]+"_%g" % p[i] for i in xrange(len(multival))]) + '/'

    os.makedirs(fold_name)
    print fold_name

    wavelength = "%03d" % p[multival.index("wavelength")]

    # Linking files
    mie_file = "%s_RH%02d_0.%s0um.mie.out" % ( params['aerosol_profile'],
                                                params['relative_humidity'],
                                                wavelength )
    os.symlink(os.path.abspath(mie_file),fold_name+"aerosol.mie.out")

    for zon in zones:
        os.symlink(os.path.abspath("fctem_wl_%s_zon_%s.dat" % (wavelength,zon)), fold_name+params['exp_name']+"_fctem_"+zon+".dat")

    ppath = os.environ['PATH'].split(os.pathsep)
    illumpath = filter(lambda s: "illumina" in s and "bin" in s, ppath)[0]
    os.symlink(os.path.abspath(illumpath+"/illumina"),fold_name+"illumina")

    # Copying layer data
    layer = p[multival.index("layer")]

    ds = MSD.Open("srtm.hdf5")
    save_bin(fold_name+params['exp_name']+"_topogra.bin", ds[layer])

    for name in ["obstd","obsth","obstf","altlp"]:
        ds = MSD.Open("%s_%s.hdf5" % (params['exp_name'],name))
        save_bin(fold_name+"%s_%s.bin" % (params['exp_name'],name), ds[layer])

    for zon in zones:
        ds = MSD.Open("%s_%s_lumlp_%s.hdf5" % (params['exp_name'],wavelength,zon))
        save_bin(fold_name+"%s_lumlp_%s.bin" % (params['exp_name'],zon), ds[layer])

    ds = MSD.Open("modis_%s.hdf5" % wavelength)
    save_bin(fold_name+"%s_reflect.bin" % params['exp_name'], ds[layer])

    # Create illumina.in
    n_space = 30
    with open(fold_name+"illumina.in",'w') as f:
        f.write("%-*s ! Input file for ILLUMINA\n" % (n_space,''))
        f.write("%-*s ! Root file name\n" % (n_space, params['exp_name']))
        f.write("%-*s ! Cell size [m]\n" % \
                (n_space, "%(s)s %(s)s" % {'s':ds.pixel_size(layer)}))
        f.write("%-*s ! Aerosol optical cross section file\n" % \
                (n_space, "aerosol.mie.out"))
        f.write("%-*s !\n" % (n_space,''))
        f.write("%-*s ! Double scattering radius [m] ; Scattering step\n" % \
                (n_space, "%g %d" % (getp("scattering_radius",multival,p,params)*1000,
                                     getp("scattering_skip",multival,p,params))))
        f.write("%-*s !\n" % (n_space,''))
        f.write("%-*s ! Wavelength [nm]\n" % (n_space, wavelength))
        f.write("%-*s ! Ground level pressure [kPa]\n" % \
                (n_space, getp("air_pressure",multival,p,params)))
        f.write("%-*s ! 500nm aerosol optical depth ; Angstrom exponent\n" % \
                (n_space, "%g %g" % (getp("aerosol_optical_depth",multival,p,params),
                                     getp("angstrom_coefficient",multival,p,params))))
        f.write("%-*s ! Number of source types\n" % (n_space, len(zones)))
        f.write("%-*s ! Contribution threshold\n" % \
                (n_space, getp("stop_limit",multival,p,params)))
        f.write("%-*s !\n" % (n_space,''))
        f.write(("%-*s ! Observer x position ; Observer y position ; "
                        "Observer elevation above ground [m] ; "
                        "Beginning cell along line of sight\n") % \
                (n_space, "%d %d %s %s" % (ds._attrs['nb_pixels']/2,
                                           ds._attrs['nb_pixels']/2,
                                           getp("observer_elevation",multival,p,params),
                                           1)))
        f.write("%-*s !\n" % (n_space,''))
        f.write("%-*s ! Elevation viewing angle ; Azimutal viewing angle\n" % \
                (n_space, "%d %d" % (getp("elevation_angles",multival,p,params),
                                     getp("azimuth_angles",multival,p,params))))
        f.write("%-*s !\n" % (n_space,''))
        f.write(("%-*s ! Slit width [m] ; Pixel size [m] ; "
                        "Focal length [m] ; Apperture diameter [m]\n") % \
                (n_space, "%g %g %g %g" % (1., 1., 1., 1.1283791671)))
        f.write("%-*s ! to get W/str/m**2, SW*PS*pi*AD**2/4/FL**2 = 1\n" % \
                (n_space,''))
        f.write("%-*s ! 1. 1. 1. 1.1283791671 is doing exactly that\n" % \
                (n_space,''))
        f.write(("%-*s ! Cloud model: 0=clear, 1=Thin Cirrus/Cirrostratus, "
                        "2=Thick Cirrus/Cirrostratus, 3=Altostratus/Altocumulus, "
                        "4=Cumulus/Cumulonimbus, 5=Stratocumulus\n") % \
                (n_space, getp("cloud_model",multival,p,params)))
        f.write("%-*s ! Minimal distance to nearest light source [m]\n" % \
                (n_space, getp("nearest_source_distance",multival,p,params)))

    # Write execute script
    with open(fold_name+"execute",'w') as f:
        f.write("#!/bin/sh\n")
        f.write("#SBATCH --job-name=Illumina\n")
        f.write("#SBATCH --time=%d:00:00\n" % params["estimated_computing_time"])
        f.write("#SBATCH --mem-per-cpu=1920\n")
        f.write("cd %s\n" % os.path.abspath(fold_name))
        f.write("umask 0011\n")
        f.write("./illumina\n")
    os.chmod(fold_name+"execute",0o777)

    # Append execution to batch list
    with open(sys.argv[1]+'/'+params['batch_file_name']+"_%d" % ((count/300)+1),'a') as f:
        f.write("cd %s\n" % os.path.abspath(fold_name))
        f.write("sbatch ./execute\n")
        f.write("sleep 0.05\n")

    count += 1
