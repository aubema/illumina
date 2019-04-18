#!/usr/bin/env python2
#
# Preprocessing of lamp inventory for Illumina
#
# Author : Alexandre Simoneau
#
# April 2019

print "Building inputs from discrete inventory."

# lamps distribution
inv_name = params['lamps_inventory']
lampsData = np.loadtxt(inv_name, usecols=range(7))
photometry = np.loadtxt(inv_name, usecols=[-2,-1], dtype=str)
domain = hdf.from_domain("domain.ini")

sources = np.unique(photometry[:,1])

print "Classifying points."

points = ddict(list)
for layer in xrange(len(domain)):
    ysize,xsize = domain[layer].shape
    for i,lamp in enumerate(lampsData):
        col,row = domain._get_col_row(lamp[:2],layer)
        if 0 <= col < xsize and 0 <= row < ysize:
            points[layer,col,row].append(i)

print "Calculating the generalized lamps."

for l in xrange(n_bins):
	for s in sources:
		np.savetxt(
			dir_name+"fctem_wl_%03d_lamp_%s.dat" % (x[l],s),
			np.concatenate([ lop[s],angles ]).reshape((2,-1)).T )

with open(dir_name+"lamps.lst",'w') as zfile:
	zfile.write('\n'.join( sources ) + '\n' )

geometry = dict()
for geo in ['obsth','obstd','obstf','altlp']:
    geometry[geo] = hdf.from_domain("domain.ini")

lumlp = dict()
for s in sources:
    for wl in x:
        lumlp[s,wl] = hdf.from_domain("domain.ini")

for key,ind in points.iteritems():
    layer,col,row = key
    lumens = lampsData[:,2][ind]

    for n,geo in izip(xrange(3,7),['obsth','obstd','obstf','altlp']):
        geometry[geo][layer][row,col] = np.average(
            lampsData[:,n][ind],
            weights=lumens
        )

    local_sources = np.unique(photometry[ind][:,1])
    for s in local_sources:
        mask = photometry[:,1][ind] == s
        fctem = np.array([
            spct[type] for type in photometry[:,0][ind][mask]
        ])
        fctem = np.sum(fctem*lumens[mask,None],0)

        y = [ np.mean(x) for x in \
            np.array_split(
        		fctem[bool_array],
        		n_bins,
        		-1
            )
        ]

        for wl in x:
            lumlp[s,wl][layer][row,col] = y

print "Saving data."

for geo,ds in geometry.iteritems():
    ds.save(dir_name+out_name+"_"+geo)

for key,ds in lumlp.iteritems():
    s,wl = key
    ds.save(dir_name+"%s_%03d_lumlp_%s" % (out_name,wl,s))
