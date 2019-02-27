#!/usr/bin/env python2

import numpy as np
from pytools import load_pgm, save_pgm

cutoff = 256
infilename = "pgms/stable_lights.pgm"
outfilename = "pgms/stable_lights_lowcut.pgm"

head,p,data = load_pgm(infilename)
data[data<(cutoff*gain+offset)] = 0.
save_pgm(outfilename,head,p,data)
