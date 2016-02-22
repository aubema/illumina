#!/usr/bin/env python

import numpy as np
import subprocess as sub
import os, re, pyfits
import matplotlib.pyplot as plt
from glob import glob
from pytools import load_pgm, save_fits
from collections import defaultdict as ddict

def parse_data(s):
    s = s.split()
    if len(s)==1:
        val = 0.
        s = s[0]
    else:
        s,val = s
    sp = re.split('--|-',s)
    wl = int(sp.pop(3)[2:])
    sp = [sp.pop(1)]+sp
    return '-'.join(sp),[wl,float(val),s]

dir = raw_input("Result folder path : ")
integ = raw_input("Integration_limits file path : ")+"/integration_limits.dat"
lims = np.genfromtxt(integ,skiprows=1)
bw =  lims[1:]-lims[:-1]
wl = (lims[1:]+lims[:-1])/2

files = glob(dir+"/*/data.txt")
files = map(lambda s:(s.split('/')[-2].split('_')[-1],s), files)

tree = ddict(lambda: ddict(list))
ftree = ddict(lambda: ddict(list))

for cat,fname in files:
    with open(fname) as f:
        data = map(parse_data,f.readlines())
    for s,d in data:
        tree[s][cat].append(d)

files = dict(files)
for cat in files:
    files[cat] = os.path.dirname(files[cat])

for a in tree:
    for b in tree[a]:
        l = sorted(tree[a][b])
        ftree[a][b] = map(list.pop, l)
        tree[a][b] = np.asarray(l)[:,1]/bw

scen = sorted(files)

try:
    plot_ratio = raw_input("Generate ratio graph [Y/n] : ")[0] in "Yy"
except IndexError:
    plot_ratio = True

if plot_ratio:
    print "Experiments : " + str(scen)[1:-1]
    basename = raw_input("Experiment to compare to : ")
    if basename not in scen:
        print "Error : Experiment '%s' not found." % basename
        plot_ratio = False

inter = plt.isinteractive()
plt.ioff()

spcts = dict()
for d in tree:
    temp = [wl]
    for s in sorted(tree[d].keys()):
        temp.append(tree[d][s])
    spcts[d] = np.asarray(temp).T
    np.savetxt(d+".dat",spcts[d],header="wavelength\t"+'\t'.join(sorted(tree[d].keys())))
    if plot_ratio:
        plt.figure()
        for s in sorted(tree[d].keys()):
            plt.plot(wl,tree[d][s]/tree[d][basename],'-',label=s+" : %.02f" % np.mean(tree[d][s]/tree[d][basename]))
        plt.legend(loc=0)
        plt.xlabel("Wavelenght (nm)")
        plt.ylabel("LP ratio after/before")
        plt.title("LP ratio to "+s+" scenario\n"+d.replace('-',' '))
        plt.xlim((np.min(wl),np.max(wl)))
        plt.savefig(d+".png")
        plt.close()

if inter:
    plt.ion()

try:
    mk_cube = raw_input("Generate fits data cubes [Y/n] : ")[0] in "Yy"
except IndexError:
    mk_cube = True

if mk_cube:
    for sc in ftree:
        for cat in ftree[sc]:
            name = cat+'-'+sc.split('-',1)[1]+".fits"
            print "Generating '%s'" % name
            cube = np.asarray(map(lambda s: load_pgm(files[cat]+"/PCL-%s.pgm"%s)[2], ftree[sc][cat]))[:,::-1,:]
            save_fits([(1,1),(1,1),(wl[0],bw[0])],cube,name)

print "\nDone."
    
