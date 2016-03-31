#!/usr/bin/env python
#
# Library of usefull function related to Illumina and PGM and FITS handling
#
# Author : Alexandre Simoneau
# unless noted otherwise
#
# March 2016

import numpy as _np
import pyfits as _fits
from scipy.interpolate import interp1d as _I

def safe_divide(a,b):
	"""Safely divide two arrays, with 0 as a result of a division by 0."""
	with _np.errstate(divide='ignore', invalid='ignore'):
		c = _np.true_divide(a,b)
		c[c == _np.inf] = 0
		c = _np.nan_to_num(c)
	return c

def LOP_norm(angles,x):
	"""Normalises 'x' as a function of theta over the full sphere.
	Uses the two first elements of 'angles' as the integration step."""
	sinx = 2*_np.pi*_np.sin(_np.deg2rad(angles))
	dtheta = angles[1]-angles[0]
	return safe_divide( x, _np.sum(x*sinx)*dtheta )

def SPD_norm(wav, norm_spct, x, factor=683.002):
	"""Normalises a spectrum 'x' with a normalisation spectrum to 'factor'.
	Both arrays must be of the same lenght.
	Uses the two first elements of 'wav' as the integration step."""
	dlambda = wav[1]-wav[0]
	return safe_divide( x, factor*_np.sum(norm_spct*x)*dlambda )

def spct_norm(wav, x):
	"""Normalises 'x' using the two first elements of 'wav' as the integration step."""
	dlambda = wav[1]-wav[0]
	return safe_divide( x, _np.sum(x)*dlambda )

def zon_norm(angles, wavelenght, zone):
	"""Normalises an Illumina zone LOP"""
	sinx = 2*_np.pi*_np.sin(_np.deg2rad(angles))
	dlambda = wavelenght[1]-wavelenght[0]
	dtheta = angles[1]-angles[0]
	return safe_divide( zone, _np.sum(zone.T*sinx)*dtheta*dlambda )

def parse_inventory(filename, n=0):
	"""Parse an inventory type file.
	Skips the first 'n' columns."""
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
	"""Returns an array of normalized zones.
	
	  theta : angles used to define the Light Output Pattern
	  lop : Light Output Pattern dictionnary
	  wl : wavelength used to define the spectrum
	  spct : Lamp spectrum dictionnary
	  ivtr : Parsed inventory"""
	return _np.asarray([ zon_norm(theta, wl, sum(l[0]*spct[l[1]]* \
		lop[l[2]][:,_np.newaxis] for l in lampData)) for lampData in ivtr ])
	
def load_pgm(filename):
	"""Opens a PGM file.

	Returns the header dictionnary, the image parameters and the data as a 3-ple.
	The image parameters are a list of the height, the width and the resolution."""
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
	"""Saves a PGM file.

	  head : A dictionnary of headers
	  p : The image parameters as returned by 'load_pgm'
	  data : Data to be saved"""
	new_gain = _np.max(data)/float(p[2])
	
	if new_gain != 0:
		data = (data / new_gain)
	data = data.astype(int)

	head['gain'] = str(new_gain)
	head['offset'] = "0.0"

	headstring = "P2\n" + '\n'.join( map(lambda s: '# '+s[0]+' '+s[1], head.iteritems() ) ) + \
			'\n' + ' '.join(map(str,p))
	_np.savetxt(filename,data,fmt="%d",header=headstring,comments='')

def load_fits(filename):
	"""Loads a FITS file.

	Returns a 2-ple containing
	  - A list of arrays defining each axis in order (x,y,z,...)
	  - The caintained data in transposed order (...,z,y,x)
	"""
	hdu = _fits.open(filename)[0]
	
	ax = [ _np.linspace( hdu.header['CRVAL%d'%(i+1)],
		hdu.header['CRVAL%d'%(i+1)]+hdu.header['CDELT%d'%(i+1)]*(hdu.header['NAXIS%d'%(i+1)]-1),
		hdu.header['NAXIS%d'%(i+1)] ) for i in xrange(hdu.header['NAXIS']) ]

	return ax,hdu.data.T[:,::-1].T
	

def save_fits(axis,data,filename):
	"""Save an array to a fits file. Must be at least 2D.
	
	  axis : a list of 2-tuple containing the base value and the increment for each axis.
	  data : data array. The dimensions must be ordered (...,z,y,x)
	  filename : name of the file to create
	"""
	hdu = _fits.PrimaryHDU()
	hdu.data = data.T[:,::-1].T
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

def load_lop(angles, filename, interp="cubic"):
	"""Load an LOP file interpolated to 'angles' and normalised.
	
	  interp : Interpolation kind

	See 'scipy.interpolate.interp1d for interpolation kinds."""
	data = _np.loadtxt(filename).T
	y = data[0] if _np.all(data[1] == angles) else \
		_I(data[1], data[0], kind=interp, bounds_error=False, fill_value=0.)
	return LOP_norm(angles,y)

def load_spct(wav, norm_spct, filename, interp="cubic", factor=683.002):
	"""Load a spectrum file interpolated to 'wav' and normalised.

	  interp : Interpolation kind
	  factor : Normalisation factor

	See 'scipy.interpolate.interp1d for interpolation kinds."""
	data = _np.loadtxt(filename,skiprows=1).T
	y = data[1] if _np.all(data[0] == wav) else \
		_I(data[0], data[1], kind=interp, bounds_error=False, fill_value=0.)(wav)
	return SPD_norm(wav, norm_spct, y, factor)

