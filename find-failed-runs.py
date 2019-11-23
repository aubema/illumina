#!/usr/bin/env python2

import os, re, errno
import fnmatch
from glob import glob
import argparse

parser = argparse.ArgumentParser(description="Find failed ILLUMINA executions.")
parser.add_argument("-e", dest="batch", action="store_true",
    help="If given, will return the executable code to rerun the failed executions.")

def recursive_glob(rootdir='.', pattern='*'):
    for root, dirnames, filenames in os.walk(rootdir):
        for filename in fnmatch.filter(filenames,pattern):
            yield os.path.join(root,filename)

p = parser.parse_args()

if p.batch:
    def failed(dirname):
        print "cd " + os.path.abspath(dirname)
        print "sbatch ./execute"
        print "sleep 0.05"
else:
    def failed(dirname):
        print dirname

for dirname in recursive_glob(pattern='illumina'):
    dirname = os.path.dirname(dirname)
    if not os.path.isfile(os.path.join(dirname,"illumina.in")):
        failed(dirname)
    else:
        with open(os.path.join(dirname,"illumina.in")) as f:
            basename = f.readlines()[1].split()[0]

        outnames = glob(os.path.join(dirname,basename+"*.out"))
        nb_in = len(glob(os.path.join(dirname,"*.in")))

        if nb_in <= 2:
            if len(outnames) == 0:
                failed(dirname)
            else:
                with open(outnames[0]) as f:
                    lines = f.readlines()
                    if len(lines) < 2 or "Sky radiance" not in lines[-2]:
                        failed(dirname)
        elif os.path.isfile(os.path.join(dirname,basename+".out")) \
            or len(outnames)+1 != nb_in:
                failed(dirname)
