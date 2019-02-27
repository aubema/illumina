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

for fname in recursive_glob(pattern='illumina'):
    fname = os.path.dirname(fname)
    filenames = glob(os.path.join(fname,"*.out"))
    filenames = [ f for f in filenames if not f.endswith(".mie.out") ]

    try:
        with open(filenames[0]) as f:
            if "Sky radiance" not in f.readlines()[-5]:
                raise IOError(errno.EIO,"Aborted run",fname)
    except (IndexError,IOError) as err:
        if p.batch:
            print "cd " + fname
            print "sbatch ./execute"
            print "sleep 0.05"
        else:
            print fname
