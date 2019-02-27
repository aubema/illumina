#!/usr/bin/env python2

# Get the size of a ILLUMINA binary
from pytools import load_bin
import sys

shape = load_bin(sys.argv[1]).shape
print shape[1], shape[0]
