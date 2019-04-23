#!/usr/bin/env python2
#
# Preprocessing of zones inventory for Illumina
#
# Author : Alexandre Simoneau
#
# April 2019

print "Building inputs from zones inventory."

# lamps distribution
inv_name = params['zones_inventory']
zonData = pt.parse_inventory(inv_name,7)

print "Calculating the generalized lamps."

# Calculate zones lamps
zones = pt.make_zones(angles, lop, wav, spct, zonData )

# Create the desired lamp files
y = np.array(map(
	np.mean,
	np.array_split(
		zones[:,:,bool_array],
		n_bins,
		-1),
	[-1]*n_bins )).transpose(1,2,0)

print "Creating files."

for l in xrange(n_bins):
	for z in xrange(len(zones)):
		np.savetxt(
			dir_name+"fctem_wl_%03d_zon_%03d.dat" % (x[l],z+1),
			np.concatenate([ y[z,:,l],angles ]).reshape((2,-1)).T )

try:
	os.symlink(
		os.path.abspath("stable_lights.hdf5"),
		dir_name+"stable_lights.hdf5" )
except OSError as e:
	if e[0] != 17:
		raise

print "Making zone properties files."

circles = hdf.from_domain("domain.ini") # Same geolocalisation

zonfile = np.loadtxt(inv_name,usecols=range(7),ndmin=2)

# zone number
for i,dat in enumerate(zonfile,1):
	circles.set_circle((dat[0],dat[1]),dat[2]*1000,i)
circles.save(dir_name+out_name+"_zone")

for n,name in izip(xrange(3,7),['obsth','obstd','obstf','altlp']):
    for i,dat in enumerate(zonfile,1):
    	circles.set_circle((dat[0],dat[1]),dat[2]*1000,dat[n])
    circles.save(dir_name+out_name+"_"+name)

with open(dir_name+"zon.lst",'w') as zfile:
	zfile.write('\n'.join( map(
		lambda n: "%03d"%n,
		xrange(1,len(zones)+1) ) ) + '\n' )

print "Inverting lamp intensity."

viirs_dat = MSD.Open(dir_name+"stable_lights.hdf5") * 1e-5 #nW/cm^2/sr -> W/m^2/sr

# Water mask
for i in xrange(len(viirs_dat)):
	viirs_dat[i][refl[i](700) < 0.01] = 0

circles = MSD.Open(dir_name+out_name+"_zone.hdf5")
zon_mask = np.empty(len(circles),dtype=object)
for i in xrange(len(zon_mask)):
	zon_mask[i] = np.arange(1,len(zones)+1)[:,None,None] == circles[i]

a = np.deg2rad(angles)
mids = np.concatenate([[a[0]],np.mean([a[1:],a[:-1]],0),[a[-1]]])
sinx = 2*np.pi*(np.cos(mids[:-1])-np.cos(mids[1:]))

# Pixel size in m^2
S = np.array([
	viirs_dat.pixel_size(i)**2 \
	for i in xrange(len(viirs_dat)) ])

# phie = DNB * S / int( R ( rho/pi Gdown + Gup ) ) dlambda
gdown = (zones * sinx[:,None])[:,angles>90].sum(1)
gup = (zones * sinx[:,None])[:,angles<70].sum(1) / sinx[angles<70].sum()

# Has to be done this way or risk MemoryError
def integral():
	for i,wl in enumerate(wav):
		arr = np.empty(len(refl),dtype=object)
		for j,r in enumerate(refl):
			arr[j] = zon_mask[j] * viirs[i] * \
				(r(wl)/np.pi * gdown[:,i,None,None] + gup[:,i,None,None])
		yield arr

phie = sum(integral()) * (wav[1]-wav[0])
# layer, zone, x, y

for i,p in enumerate(phie):
	phie[i] = pt.safe_divide(viirs_dat[i] * S[i], p)

wl_bin = np.array_split(wav[bool_array],n_bins,-1)
fctem_bin = np.array_split(zones[:,:,bool_array],n_bins,-1)
# bin, zone, ang, wl

for n in xrange(n_bins):
	ratio = (fctem_bin[n].mean(-1) * sinx).sum(1)
	for z,r in enumerate(ratio):
		new = hdf.from_domain("domain.ini")
		for i in xrange(len(phie)):
			new[i] = (phie[i][z] * r)
		new.save(dir_name+"%s_%03d_lumlp_%03d" % (out_name,x[n],z+1))
