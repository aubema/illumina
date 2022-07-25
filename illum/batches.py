#!/usr/bin/env python3
#
# batch processing
#
# Author : Alexandre Simoneau
# unless noted otherwise
#
# March 2021

import os
import shutil
from collections import ChainMap, OrderedDict
from functools import partial
from glob import glob
from itertools import product

import click
import illum
import numpy as np
import yaml
from illum import MultiScaleData as MSD
from illum.pytools import save_bin
from progressbar import progressbar

progress = partial(progressbar, redirect_stdout=True)


def input_line(val, comment, n_space=30):
    value_str = " ".join(str(v) for v in val)
    comment_str = " ; ".join(comment)
    return "%-*s ! %s" % (n_space, value_str, comment_str)


def MSDOpen(filename, cached={}):
    if filename in cached:
        return cached[filename]
    ds = MSD.Open(filename)
    cached[filename] = ds
    return ds


@click.command(name="batches")
@click.argument("input_path", type=click.Path(exists=True), default=".")
@click.argument("batch_name", required=False)
@click.option(
    "-c",
    "--compact",
    is_flag=True,
    help="If given, will chain similar executions. Reduces the overall number "
    "of runs at the cost of longuer individual executions.",
)
@click.option(
    "-N",
    "--batch_size",
    type=int,
    default=300,
    show_default=True,
    help="Number of runs per produced batch file.",
)
def CLI_batches(input_path, compact, batch_size, batch_name=None):
    """Makes the execution batches.

    INPUT_PATH is the path to the folder containing the inputs.

    BATCH_NAME is an optional name for the produced batch files.
    It overwrites the one defined in 'inputs_params.in' is given.
    """
    batches(input_path, compact, batch_size, batch_name)


def batches(input_path=".", compact=False, batch_size=300, batch_name=None):
    os.chdir(input_path)

    with open("inputs_params.in") as f:
        params = yaml.safe_load(f)

    if batch_name is not None:
        params["batch_file_name"] = batch_name

    for fname in glob("%s*" % params["batch_file_name"]):
        os.remove(fname)

    exp_name = params["exp_name"]

    ds = MSD.Open(glob("*.hdf5")[0])

    # Add wavelength
    params["wavelength"] = np.loadtxt("wav.lst", ndmin=1).tolist()
    bandwidth = (params["lambda_max"] - params["lambda_min"]) / params[
        "nb_bins"
    ]

    wls = params["wavelength"]
    refls = np.loadtxt("refl.lst", ndmin=1).tolist()

    for pname in ["layer", "observer_coordinates"]:
        if len(params[pname]) == 1:
            params[pname] = params[pname][0]

    with open("lamps.lst") as f:
        lamps = f.read().split()

    if os.path.isfile("brng.lst"):
        brng = np.loadtxt("brng.lst", ndmin=1)

    # Clear and create execution folder
    dir_name = "exec" + os.sep
    shutil.rmtree(dir_name, True)
    os.makedirs(dir_name)

    count = 0
    multival = [k for k in params if isinstance(params[k], list)]
    multival = sorted(multival, key=len, reverse=True)  # Semi-arbitrary sort
    param_space = [params[k] for k in multival]

    N = np.product([len(p) for p in param_space])
    for param_vals in progress(product(*param_space), max_value=N):
        local_params = OrderedDict(zip(multival, param_vals))
        P = ChainMap(local_params, params)
        if (
            "azimuth_angle" in multival
            and P["elevation_angle"] == 90
            and params["azimuth_angle"].index(P["azimuth_angle"]) != 0
        ):
            continue

        if os.path.isfile("brng.lst"):
            obs_index = (
                0
                if "observer_coordinates" not in multival
                else params["observer_coordinates"].index(
                    P["observer_coordinates"]
                )
            )
            bearing = brng[obs_index]
        else:
            bearing = 0

        coords = "%6f_%6f" % P["observer_coordinates"]
        if "observer_coordinates" in multival:
            P["observer_coordinates"] = coords

        if compact:
            fold_name = (
                dir_name
                + os.sep.join(
                    "%s_%s" % (k, v)
                    for k, v in local_params.items()
                    if k in ["observer_coordinates", "wavelength", "layer"]
                )
                + os.sep
            )
        else:
            fold_name = (
                dir_name
                + os.sep.join(
                    "%s_%s" % (k, v) for k, v in local_params.items()
                )
                + os.sep
            )

        unique_ID = "-".join("%s_%s" % item for item in local_params.items())
        wavelength = "%g" % P["wavelength"]
        layer = P["layer"]
        reflectance = refls[wls.index(P["wavelength"])]

        if not os.path.isdir(fold_name):
            os.makedirs(fold_name)
            # Linking files
            mie_file = "%s_%s.txt" % (params["aerosol_profile"], wavelength)
            os.symlink(
                os.path.relpath(mie_file, fold_name), fold_name + "aerosol.txt"
            )
            layer_file = "%s_%s.txt" % (params["layer_type"], wavelength)
            os.symlink(
                os.path.relpath(layer_file, fold_name), fold_name + "layer.txt"
            )

            os.symlink(
                os.path.relpath("MolecularAbs.txt", fold_name),
                fold_name + "MolecularAbs.txt",
            )

            for i, lamp in enumerate(lamps, 1):
                os.symlink(
                    os.path.relpath(
                        "fctem_wl_%s_lamp_%s.dat" % (wavelength, lamp),
                        fold_name,
                    ),
                    fold_name + exp_name + "_fctem_%03d.dat" % i,
                )

            illumpath = os.path.dirname(illum.__path__[0])
            os.symlink(
                os.path.abspath(illumpath + "/bin/illumina"),
                fold_name + "illumina",
            )

            # Copying layer data
            obs_fold = os.path.join("obs_data", coords, str(layer))

            os.symlink(
                os.path.relpath(os.path.join(obs_fold, "srtm.bin"), fold_name),
                fold_name + exp_name + "_topogra.bin",
            )

            os.symlink(
                os.path.relpath(
                    os.path.join(obs_fold, "origin.bin"), fold_name
                ),
                fold_name + "origin.bin",
            )

            for name in ["obstd", "obsth", "obstf", "altlp"]:
                os.symlink(
                    os.path.relpath(
                        os.path.join(obs_fold, "%s_%s.bin" % (exp_name, name)),
                        fold_name,
                    ),
                    fold_name + "%s_%s.bin" % (exp_name, name),
                )

            for i, lamp in enumerate(lamps, 1):
                os.symlink(
                    os.path.relpath(
                        os.path.join(
                            obs_fold,
                            "%s_%s_lumlp_%s.bin"
                            % (exp_name, wavelength, lamp),
                        ),
                        fold_name,
                    ),
                    fold_name + "%s_lumlp_%03d.bin" % (exp_name, i),
                )

        # Create illum-health.in
        input_data = (
            (("", "Input file for ILLUMINA"),),
            ((exp_name, "Root file name"),),
            (
                (ds.pixel_size(layer), "Cell size along X [m]"),
                (ds.pixel_size(layer), "Cell size along Y [m]"),
            ),
            ((5, "Number of light sources types"),),
            ((512), "Domain size"),  # TODO: Put real value in
            (
                (P["aerosol_optical_depth"], "Aerosol optical depth at 500nm"),
                (P["angstrom_coefficient"], "Angstrom exponent"),
                (P["aerosol_height"], "Aerosol scale height [m]"),
            ),
            ((wavelength, "Wavelength [nm]"), (bandwidth, "Bandwidth [nm]")),
            (
                (reflectance, "Street reflectance"),  # TODO: Put real value in
                (
                    reflectance,
                    "Building reflectance",
                ),  # TODO: Put real value in
                (reflectance, "Ground reflectance"),
            ),  # TODO: Put real value in
            ((P["air_pressure"], "Ground level pressure [kPa]"),),
            (
                (
                    P["observer_elevation"],
                    "Observer elevation above ground [m]",
                ),
            ),
        )

        with open(fold_name + "illum-health.in", "w") as f:
            lines = (input_line(*zip(*line_data)) for line_data in input_data)
            f.write("\n".join(lines))

        # Write execute script
        if not os.path.isfile(fold_name + "execute"):
            with open(fold_name + "execute", "w") as f:
                f.write("#!/bin/sh\n")
                f.write("#SBATCH --job-name=Illumina\n")
                f.write(
                    "#SBATCH --time=%d:00:00\n"
                    % params["estimated_computing_time"]
                )
                # f.write("#SBATCH --mem=2G\n") # TODO: Evaluate ram requirements dynamically
                f.write("cd %s\n" % os.path.abspath(fold_name))
                f.write("umask 0011\n")
            os.chmod(fold_name + "execute", 0o777)

            # Append execution to batch list
            with open(
                f"{params['batch_file_name']}_{(count//batch_size)+1}", "a"
            ) as f:
                f.write("cd %s\n" % os.path.abspath(fold_name))
                f.write("sbatch ./execute\n")
                f.write("sleep 0.05\n")

            count += 1

        # Add current parameters execution to execution script
        with open(fold_name + "execute", "a") as f:
            f.write("cp %s.in illum-health.in\n" % unique_ID)
            f.write("./illum-health\n")
            f.write("mv %s.out %s_%s.out\n" % (exp_name, exp_name, unique_ID))
            f.write(
                "mv %s_pcl.bin %s_pcl_%s.bin\n"
                % (exp_name, exp_name, unique_ID)
            )

    print("Final count:", count)

    print("Done.")
