#!/usr/bin/env python3

import fnmatch
import os
from glob import glob

import click


def recursive_glob(rootdir=".", pattern="*"):
    for root, dirnames, filenames in os.walk(rootdir):
        for filename in fnmatch.filter(filenames, pattern):
            yield os.path.join(root, filename)


@click.command(name="failed")
@click.option(
    "-e",
    "--executable",
    is_flag=True,
    help="If given, returns the executable code to rerun failed executions.",
)
def CLI_failed(executable):
    "Find failed ILLUMINA executions."
    failed(executable)


def failed(executable=False):
    if executable:

        def failed(dirname):
            print("cd " + os.path.abspath(dirname))
            print("sbatch ./execute")
            print("sleep 0.05")

    else:

        def failed(dirname):
            print(dirname)

    for dirname in recursive_glob(pattern="illumina"):
        dirname = os.path.dirname(dirname)
        if not os.path.isfile(os.path.join(dirname, "illumina.in")):
            failed(dirname)
        else:
            with open(os.path.join(dirname, "illumina.in")) as f:
                basename = f.readlines()[1].split()[0]

            outnames = glob(os.path.join(dirname, basename + "*.out"))
            nb_in = len(glob(os.path.join(dirname, "*.in")))

            if nb_in <= 2:
                if len(outnames) == 0:
                    failed(dirname)
                else:
                    with open(outnames[0]) as f:
                        lines = f.readlines()
                    if len(lines) < 2 or "Diffuse radiance" not in lines[-2]:
                        failed(dirname)
            elif (
                os.path.isfile(os.path.join(dirname, basename + ".out"))
                or len(outnames) + 1 != nb_in
            ):
                failed(dirname)
