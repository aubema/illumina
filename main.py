#!/usr/bin/env python3

import click
from importlib import import_module

@click.group()
def illum():
    """Illumina model.

    See 'illum COMMAND --help' to read about specific subcommand.
    """
    pass # Entry point

functions = (
    ("defineDomain","domain"),
    ("Illuminutils","warp"),
    ("make_inputs","inputs"),
    ("makeBATCH","batches"),
    ("find-failed-runs","failed"),
    ("extract-output-data","extract")
)

for module_name,method in functions:
    module = import_module(module_name)
    illum.add_command(getattr(module,method))
