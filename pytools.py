#!/usr/bin/env python

import numpy as _np
import pyfits as _fits

def safe_divide(a,b):
	with _np.errstate(divide='ignore', invalid='ignore'):
		c = _np.true_divide(a,b)
		c[c == _np.inf] = 0
		c = _np.nan_to_num(c)
	return c

def LOP_norm(angles,x):
	sinx = 2*_np.pi*_np.sin(_np.deg2rad(angles))
	dtheta = angles[1]-angles[0]
	return safe_divide( x, _np.sum(x*sinx)*dtheta )

def SPD_norm(wavelenght, norm_spct, x):
	dlambda = wavelenght[1]-wavelenght[0]
	return safe_divide( x, 683.002*_np.sum(norm_spct*x)*dlambda )

def zon_norm(angles, wavelenght, zone):
	sinx = 2*_np.pi*_np.sin(_np.deg2rad(angles))
	dlambda = wavelenght[1]-wavelenght[0]
	dtheta = angles[1]-angles[0]
	return safe_divide( zone, _np.sum(zone.T*sinx)*dtheta*dlambda )

def parse_inventory(filename,n):
	"""Parse an inventory type file."""
	def lamp_norm(lampsData):
		trans = map(list, zip(*lampsData))
		norm = sum(trans[0])
		if norm != 0.:
			trans[0] = [n/norm for n in trans[0]]
		return map(list, zip(*trans))
		
	with open(filename) as inv_file:
		zonData = strip_comments(inv_file.readlines())
	zonData = map(lambda s: s.split()[n:], zonData)
	zonData = [ map(lambda s: s.split('_'),i) for i in zonData ]
	zonData = [ map(lambda s: [float(s[0]),s[1],s[2]],i) for i in zonData ]
	zonData = map(lamp_norm, zonData)
	
	return zonData
	
def make_zones(theta, lop, wl, spct, ivtr):
	return _np.asarray([ zon_norm(theta, wl, sum(l[0]*spct[l[1]]* \
		lop[l[2]][:,_np.newaxis] for l in lampData)) for lampData in ivtr ])
	
def load_pgm(filename):
	with open(filename) as file:
		data = file.read().split('\n')[:-1]

	head = ( s for s in data if s[0]=='#' )
	head = { s.split()[1]:s.split()[2] for s in head }

	gain = float(head.setdefault('gain',0))
	offset = float(head.setdefault('offset',0))

	header = data[:len(head)+2]
	data = _np.loadtxt(filename,skiprows=len(head)+2,ndmin=2)
	data = data*gain+offset
	
	p = map(int,header[-1].split())

	return head,p,data.reshape(p[1::-1])

def save_pgm(filename,head,p,data):
	new_gain = _np.max(data)/float(p[2])
	
	if new_gain != 0:
		data = (data / new_gain)
	data = data.astype(int)

	head['gain'] = str(new_gain)
	head['offset'] = "0.0"

	headstring = "P2\n" + '\n'.join( map(lambda s: '# '+s[0]+' '+s[1], head.iteritems() ) ) + \
			'\n' + ' '.join(map(str,p))
	_np.savetxt(filename,data,fmt="%d",header=headstring,comments='')

def save_fits(axis,data,filename):
	"""Save a data cube to a fits file.
	
	  axis : a list of 2-tuple containing the base value and the increment for each axis
	"""
	hdu = _fits.PrimaryHDU()
	hdu.data = data
	for i in xrange(len(axis)):
		hdu.header['CRPIX%d'%(i+1)] = 1
		hdu.header['CRVAL%d'%(i+1)] = axis[i][0]
		hdu.header['CDELT%d'%(i+1)] = axis[i][1]
	hdu.writeto(filename,clobber=True)

def strip_comments(item, token='#'):
    """Generator. Strips comments and whitespace from input lines.

    This generator strips comments, leading/trailing whitespace, and
    blank lines from its input.

    Arguments:
        item (obj):  Object to strip comments from.
        token (str, optional):  Comment delimiter.  Defaults to ``#``.

    Yields:
        str:  Next non-blank line from ``item`` with comments and
            leading/trailing whitespace removed.

    Credits: Doug R., StackOverflow
    """

    for line in item:
        s = line.split(token, 1)[0].strip()
        if s != '':
            yield s

