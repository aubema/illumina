#!/usr/bin/env python3

from importlib import import_module

import click

import illum


@click.group()
@click.version_option(illum.__version__, prog_name="Illumina model")
def illum():
    r"""Illumina model.

    See 'illum COMMAND --help' to read about specific subcommand.
    The user's guide is available at:
    http://dome.obsand.org:2080/wiki/pmwiki.php/Prof/IlluminaGuidev22

    Below is detailled the typical use case:

    1. init -- Prep the experiment folder.

    2. domain -- Define the simulation domain once the desired values
    are set in `domain_params.in`.

    3. warp -- Process the satellite imagery.

    4a. inputs -- Prepare the execution folder once the desired values
    are set in `inputs_params.in`.

    4b. alternate -- (Optionnal) Prepare conversion scenarios to be executed.

    At this point, the `Inputs` folder can be moved to a compute cluster.

    5. batches -- Prepare the model to run in parallel.

    At this point, the produced batches files can be executed.

    6a. failed -- (Optionnal) Identify failed execution that need to be rerun.

    6b. extract -- Extract the results from the executions.

    At this point, the analysis of the data is left to the user.
    We recommend Python and provide some functions to assist in that regard
    in the `MultiScaleData` module.


    One may want to convert the custom HDF5 format used by Illumina to a more
    standard format for use with GIS programs.
    This can be acheived with `convert`.
    """
    pass  # Entry point


functions = (
    (".alternate", "CLI_alternate"),
    (".batches", "CLI_batches"),
    (".convert", "CLI_convert"),
    (".domain", "CLI_domain"),
    (".failed", "CLI_failed"),
    (".extract", "CLI_extract"),
    (".init", "CLI_init"),
    (".inputs", "CLI_inputs"),
    (".warp", "CLI_warp"),
)

for module_name, method in functions:
    module = import_module(module_name, package="illum")
    illum.add_command(getattr(module, method))
