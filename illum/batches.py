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
import numpy as np
import yaml
from progressbar import progressbar

import illum

progress = partial(progressbar, redirect_stdout=True)


def input_line(val, comment, n_space=30):
    value_str = " ".join(str(v) for v in val)
    comment_str = " ; ".join(comment)
    return "%-*s ! %s" % (n_space, value_str, comment_str)


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

    ds = illum.utils.load_geotiff(glob("*.tif"))

    # Add wavelength
    spectral_bands = np.loadtxt("wav.lst", ndmin=2)
    params["wavelength"] = spectral_bands[:, 0].tolist()

    wls = params["wavelength"]
    refls = np.loadtxt("refl.lst", ndmin=2).tolist()

    with open("lamps.lst") as f:
        lamps = f.read().split()

    # Clear and create execution folder
    dir_name = "exec"
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

        filter = ["wavelength"] if compact else local_params.keys()
        fold_name = os.path.join(
            dir_name,
            ("%s_%s" % (k, v) for k, v in local_params.items() if k in filter),
            "",
        )

        unique_ID = "-".join("%s_%s" % item for item in local_params.items())
        wavelength = "%g" % P["wavelength"]
        layer = P["layer"]
        reflectance = refls[wls.index(P["wavelength"])]
        bandwidth = spectral_bands[wls.index(P["wavelength"]), 1]

        if not os.path.isdir(fold_name):
            os.makedirs(fold_name)

            binfiles = (
                "topogra",
                "obsth",
                "altlp",
                "azimu",
                "lmpty",
                "gndty",
                "lumlp",
            )

            for name in binfiles:
                name += ".bin"
                os.symlink(
                    os.path.relpath(name, fold_name),
                    os.path.join(fold_name, name),
                )

            for fname in glob("fctem_*.dat"):
                os.symlink(
                    os.path.relpath(fname, fold_name),
                    os.path.join(fold_name, "_".join(fname.split("_")[::2])),
                )

            os.symlink(
                os.path.abspath(os.path.join(illum.path, "bin", "illumina")),
                os.path.join(fold_name, "illum-health"),
            )

        # Create illum-health.in
        input_data = (
            (("", "Input file for ILLUMINA"),),
            ((exp_name, "Root file name"),),
            (
                # TODO: Put real sizes in
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
                (reflectance[0], "Street"),
                (reflectance[1], "Building"),
                (reflectance[2], "Ground reflectance"),
            ),
            ((P["air_pressure"], "Ground level pressure [kPa]"),),
            (
                (
                    P["observer_elevation"],
                    "Observer elevation above ground [m]",
                ),
            ),
        )

        with open(os.path.join(fold_name, "in_" + unique_ID), "w") as f:
            lines = (input_line(*zip(*line_data)) for line_data in input_data)
            f.write("\n".join(lines))

        # Write execute script
        exec_file = os.path.join(fold_name, "execute")
        if not os.path.isfile(exec_file):
            with open(exec_file, "w") as f:
                f.write("#!/bin/sh\n")
                f.write("#SBATCH --job-name=Illumina\n")
                f.write(
                    "#SBATCH --time=%d:00:00\n"
                    % params["estimated_computing_time"]
                )
                # TODO: Evaluate ram requirements dynamically
                # f.write("#SBATCH --mem=2G\n")
                f.write("cd %s\n" % os.path.abspath(fold_name))
                f.write("umask 0011\n")
            os.chmod(exec_file, 0o777)

            # Append execution to batch list
            with open(
                f"{params['batch_file_name']}_{(count//batch_size)+1}", "a"
            ) as f:
                f.write("cd %s\n" % os.path.abspath(fold_name))
                f.write("sbatch ./execute\n")
                f.write("sleep 0.05\n")

            count += 1

        # TODO: Manage outputs properly
        # Add current parameters execution to execution script
        with open(exec_file, "a") as f:
            f.write("cp in_%s illum-health.in\n" % unique_ID)
            f.write("./illum-health\n")
            f.write("mv %s.out %s_%s.out\n" % (exp_name, exp_name, unique_ID))
            f.write(
                "mv %s_pcl.bin %s_pcl_%s.bin\n"
                % (exp_name, exp_name, unique_ID)
            )

    print("Final count:", count)
    print("Done.")
