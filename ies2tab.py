#!/usr/bin/python2
# cedi est un commentaire
import numpy as np
import sys

nV,nH = 0,0
aV,aH = [],[]
data = []
dump = []

filename = raw_input("input filename : ")

with open(filename) as f:
    for line in f.readlines():
        vals = line.split()

        try:
            vals = [ int(float(v)) for v in vals ]
        except ValueError:
            continue

        if nV == 0:
            nV,nH = vals[3],vals[4]
        elif len(dump) == 0:
            dump.extend(vals)
        elif len(aV) < nV:
            aV.extend(vals)
        elif len(aH) < nH:
            aH.extend(vals)
        else:
            data.extend(vals)

aV = np.array(aV)
aH = np.array(aH)
data = np.array(data).reshape((nH,nV))

outname = raw_input("out filename : ")
with open(outname,'w') as f:
	f.write("%d %d\n" % (nV,nH))
	f.write(' '.join(map(str,aV))+'\n')
	f.write(' '.join(map(str,aH))+'\n')
	f.write('\n'.join([' '.join(map(str,d)) for d in data])+'\n')
