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

sources = np.unique([ lamp[2] for zd in zonData for lamp in zd ])

for l in xrange(n_bins):
	for s in sources:
		np.savetxt(
			dir_name+"fctem_wl_%03d_lamp_%s.dat" % (x[l],s),
			np.concatenate([ lop[s],angles ]).reshape((2,-1)).T )

with open(dir_name+"lamps.lst",'w') as zfile:
	zfile.write('\n'.join( sources ) + '\n' )

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

print "Inverting lamp intensity."

# Calculate zones lamps
zones = pt.make_zones(angles, lop, wav, spct, zonData, sources )

viirs_dat = MSD.Open("stable_lights.hdf5") * 1e-5 #nW/cm^2/sr -> W/m^2/sr

water_mask = MSD.Open("water_mask.hdf5")
for i,wm in enumerate(water_mask):
	viirs_dat[i][wm==0] = 0.

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
Gdown = np.tensordot(zones[:,:,angles>90],sinx[angles>90],axes=([2],[0]))
Gup = np.tensordot(zones[:,:,angles<70],sinx[angles<70],axes=([2],[0])) \
	/ sinx[angles<70].sum()
integral = np.sum(viirs * (Gdown*refl/np.pi + Gup),(1,2)) * (wav[1]-wav[0])

phie = [
	pt.safe_divide(
		viirs_dat[i] * S[i],
		np.sum(zon_mask[i] * integral[:,None,None], 0)
	) \
	for i in xrange(len(S))
]

ratio = [
	fctem_bin.mean(-1) \
	for fctem_bin in np.array_split(
		(zones*sinx[:,None]).sum(2)[:,:,bool_array],
		n_bins,-1
	)
]

for n in xrange(n_bins):
	r = [
		np.sum(
			zon_mask[layer][:,None] * \
			ratio[n][:,:,None,None],
			0
		) for layer in xrange(len(phie))
	]
	for i,s in enumerate(sources):
		new = hdf.from_domain("domain.ini")
		for layer in xrange(len(new)):
			new[layer] = phie[layer] * r[layer][i]
		new.save(dir_name+"%s_%03d_lumlp_%s" % (out_name,x[n],s))
