import numpy as np, matplotlib.pyplot as plt, re, os
from subprocess import call

# Load data
# IMPORTANT : x axis of all similar data must be the same

print "Loading data..."

# Angular distribution (normalised to 1)
u0,angles = np.loadtxt("Lights/0pcUPLIGHT.dat").T
u5 = np.loadtxt("Lights/5pcUPLIGHT.dat")[:,0]
u10 = np.loadtxt("Lights/10pcUPLIGHT.dat")[:,0]
u15 = np.loadtxt("Lights/15pcUPLIGHT.dat")[:,0]
u20 = np.loadtxt("Lights/20pcUPLIGHT.dat")[:,0]

upL = { 'x':angles, 0:u0/sum(u0), 5:u5/sum(u5), 10:u10/sum(u10), 15:u15/sum(u15), 20:u20/sum(u20) }

# Spectral distribution (normalised with scotopric vision to 1)
wavelenght, scotopic = np.loadtxt("Lights/scotopic.dat", skiprows=1).T
hps = np.loadtxt("Lights/HPS_Helios.dat", skiprows=1)[:,1]
lps = np.loadtxt("Lights/LPS_Dline.dat", skiprows=1)[:,1]
mh  = np.loadtxt("Lights/MH.dat", skiprows=1)[:,1]
flu = np.loadtxt("Lights/Fluorescent_4kK.dat", skiprows=1)[:,1]

lamps = { 'x':wavelenght, 'H':hps/sum(scotopic*hps), 'L':lps/sum(scotopic*lps), 'M':mh/sum(scotopic*mh), 'F':flu/sum(scotopic*flu) }

# lamps distribution
with open("inventaire.txt") as file:
	zonData = file.readlines()
zonData = map(lambda s: s.split('\t',1)[1].strip().split(), zonData[1:])

p = re.compile(ur"(\d+)(\D)(\d+)")
zonData = [ map(lambda s: re.match(p,s).groups(),i) for i in zonData ]
zonData = [ map(lambda s: [int(s[0]),s[1],int(s[2])],i) for i in zonData ]

print "Calculating the generalized lamps..."

# Calculate zones lamps
zones = np.array([ sum([ (l[0]/100.)*lamps[l[1]]*upL[l[2]][:,np.newaxis] for l in lampData]) for lampData in zonData])

for i in xrange(len(zones)):
	np.savetxt("Lamps/zone%i_lamp.dat"%(i+1),zones[i])

call(["gnuplot","zones_plot.gp"])

print "Splitting in a few wavelenghts..."
n = int(raw_input("    Number of wavelenghts to use : "))

# Create the desired lamp files
x = np.array(map(np.mean,np.array_split(lamps['x'],n,-1),[-1]*n))
y = np.array(map(np.sum,np.array_split(zones,n,-1),[-1]*n)).transpose(1,2,0)

print "Creating files..."

for l in xrange(n):
	os.makedirs("Intrants/%0.01fnm"%x[l])
	for z in xrange(len(zones)):
		np.savetxt( "Intrants/%0.01fnm/lamp_zon%02d.dat"%(x[l],z+1), np.concatenate([ y[z,:,l],upL['x'] ]).reshape((2,-1)).T )

print "Done."


