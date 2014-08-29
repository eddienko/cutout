
import os, sys, time, json, md5
import argparse
import subprocess
import glob
import math
from multiprocessing import Pool, Process
import StringIO

import subprocess
from subprocess import Popen

# To make tar giles
import tarfile

try:
    import pyfits
except ImportError:
    print 'ERROR: pyfits modules not installed'
    sys.exit(1)

import numpy
import numpy as np

# Matplotlib
try:
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
    from matplotlib.figure import Figure
    from matplotlib import rcParams
    from matplotlib import cm as colormap
except ImportError:
    print 'ERROR: matplotlib module not installed'
    sys.exit(1)

# Kapteyn
try:
    from kapteyn import maputils, wcsgrat, wcs
    from kapteyn.mplutil import VariableColormap
except ImportError:
    print 'ERROR: kapteyn module not installed'
    sys.exit(1)

# Python Imaging Library
try:
    from PIL import Image as PILImage
    from PIL import ImageOps as PILImageOps
except ImportError:
    print 'ERROR: python imaging module (PIL) not installed'
    sys.exit(1)

import astats
import img_scale

CUTOUT='/Users/eglez/Development/Django/casu/bin/cutout'

hdr="""BITPIX  =                  -32 / number of bits per data pixel
NAXIS   =                    2 / number of data axes
NAXIS1  =               %(naxis)6d / length of data axis 1
NAXIS2  =               %(naxis)6d / length of data axis 2
CRVAL1  =        %(racen)f / Reference value on 1st axis in primary WCS
CRVAL2  =        %(deccen)f / Reference value on 2nd axis in primary WCS
CRPIX1  =               %(crpix)f / Reference pixel on 1st axis in primary WCS
CRPIX2  =               %(crpix)f / Reference pixel on 2nd axis in primary WCS
CD1_1   =           -%(cdelt)e / Transformation matrix for primary WCS
CD2_2   =           %(cdelt)e / Transformation matrix for primary WCS
CTYPE1  = 'RA---TAN'           / &
CTYPE2  = 'DEC--TAN'           / &
RADESYS = 'ICRS    '           / Reference frame for RA/DEC values
"""



def make_hdr(fileName, hduNum, cenra, cendec, psize, md5hash):
    """Given the input FITS file name and extension, coordinates of the center of the stamp and its size and output
    name it creates a FITS header with the WCS describing the cutout.

    :param fileName:
    :param hduNum:
    :param cenra:
    :param cendec:
    :param psize:
    :param md5hash:
    :return:
    """

    # Open the FITS file and read the header
    fh = pyfits.open(fileName)
    h = fh[hduNum].header

    # Extract the relevant keywords from the header
    keys=['CTYPE1','CTYPE2','CRVAL1','CRVAL2','CRPIX1','CRPIX2','CD1_1','CD2_2', 'CD2_1','CD1_2','PV2_1','PV2_2',
          'PV2_3','PV2_4','PV2_5','NAXIS','NAXIS1','NAXIS2', 'CDELT1', 'CDELT2', 'CROTA1', 'CROTA2']
    header = {}
    header['NAXIS'] = 2
    for k in keys:
        if k in h:
            header[k]=h[k]

    fh.close()

    # Construct WCS structure and create header of stamp
    w = wcs.Projection(header)
    a, d = w.crval
    x1, y1 = w.topixel((a,d))
    x2, y2 = w.topixel((a,d+1/60.))
    scl = 1./60./math.sqrt((x2-x1)**2+(y2-y1)**2)
    h = {}
    h['cdelt'] = scl
    h['naxis'] = int(psize/scl/3600.)
    h['crpix'] = int(psize/scl/3600.)/2.0
    h['racen'] = cenra
    h['deccen'] = cendec

    # Write header to file
    open('%s.hdr' % md5hash, 'w').write(hdr % h)


def make_fits(fileName, md5hash):
    """Run the external application to make a postage stamp

    :param fileName:
    :param md5hash:
    :return:
    """

    args = [CUTOUT, fileName, "%s.hdr" % md5hash, "!%s.fit" % md5hash]
    p = subprocess.call(args)


def make_png(md5hash, cmap='gray', sigma=1.0, crop=False, options='', annotation=''):
    """Create a PNG stamp from a FITS file.

    :param md5hash:
    :param cmap:
    :param sigma:
    :param crop:
    :param options:
    :param annotation:
    :return:
    """
    outFile = md5hash

    # Read
    fitsim = maputils.FITSimage(outFile+'.fit')
    fig = Figure(figsize=(5,5), facecolor="#ffffff")
    canvas = FigureCanvas(fig)
    frame = fig.add_subplot(1,1,1)
    
    # Work out scaling
    dat = fitsim.dat.ravel()
    dat = dat.compress(dat>0.0)
    if len(dat)>10:
        z1, z2 = astats.zscale(dat, nsample=min(1000,len(dat)))
    else:
        z1 = z2 = 0.0

    # Options
    cmapinverse=False
    if ('H' in options):
        cmap='hot'
    if ('W' in options):
        cmap='winter'
    if ('I' in options):
        cmap=cmap+'_r'


    annim = fitsim.Annotatedimage(frame, cmap=cmap.replace('_r',''), clipmin=z1, clipmax=z2*sigma)
    if cmap.find('_r')>0:
        annim.cmap.set_inverse(True)
    annim.Image()
    
    # Work out grid
    a,d = fitsim.convproj.toworld(((0,0),(0,fitsim.hdr['NAXIS2'])))
    b = numpy.array([5,10,15,20,25,30,35,40,45,50,55,60])
    j = b.searchsorted(int(math.sqrt((a[0]-a[1])**2+(d[0]-d[1])**2)*3600.)/6)
    deltay = b[j]
    
    grat = annim.Graticule( 
                    starty=fitsim.hdr['CRVAL2'], deltay=deltay/3600.,
                    startx=fitsim.hdr['CRVAL1'], deltax=deltay/3600./math.cos(math.radians(d[0])))
    grat.setp_lineswcs0(visible=False)
    grat.setp_lineswcs1(visible=False)
    grat.setp_plotaxis(0, mode=0, label='')
    grat.setp_plotaxis(1, mode=0, label='')
    grat.setp_plotaxis(2, mode=0, label='')
    grat.setp_plotaxis(3, mode=0, label='')
    grat.setp_tick(visible=False)
    grat.setp_tickmark(color='#99FF99', markersize=4, markeredgewidth=2)
    
    # Plot center cross
    xcen, ycen = fitsim.hdr['NAXIS1']/2.+0.5,fitsim.hdr['NAXIS2']/2.+0.5
    frame.plot([xcen,xcen],[ycen/6.,ycen-ycen/6.],linewidth=1,color='#99FF99')
    frame.plot([xcen,xcen],[ycen+ycen/6.,ycen*2-ycen/6.],linewidth=1,color='#99FF99')
    frame.plot([xcen/6.,xcen-xcen/6.],[ycen,ycen],linewidth=1,color='#99FF99')
    frame.plot([xcen+ycen/6.,xcen*2-ycen/6.],[ycen,ycen],linewidth=1,color='#99FF99')
    
    annim.plot()

    if (annotation):
        frame.fill([0.001,0.999,0.999,0.001], [0.999,0.999,0.90-0.08, 0.90-0.08], transform = frame.transAxes, edgecolor='none', facecolor='#000000', alpha=0.4)
        i = 0
        for item in annotation.split('\\n'):
            frame.text(0.04,0.92-i*0.06, item, transform = frame.transAxes, color='#FFFFFF')
            i = i + 1

    canvas.draw()
    if crop:
        size = canvas.get_renderer().get_canvas_width_height()
        buf = canvas.tostring_rgb()
        im = PILImage.fromstring('RGB', map(int, size), buf, 'raw', 'RGB', 0, 1)
        im2 = im.convert('L'); im2 = PILImageOps.invert(im2)
        im=im.crop(im2.getbbox())
        im.save(outFile+'.png', 'PNG')
    else:
        s = StringIO.StringIO()
        canvas.print_png(s)
        s.flush()
        s.seek(0)
        open(outFile+'.png','w').write(s.read())
    time.sleep(0.25)
        
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='MapMaker')
    parser.add_argument('inputFile', type=str,
                       help='Input FITS file')
    parser.add_argument('ra', type=float,
                       help='R.A.')
    parser.add_argument('dec', type=float,
                       help='Dec')
    parser.add_argument('psize', type=float,
                       help='Stamp Size')
    parser.add_argument('--hdu', dest='hdu', type=int,
                       default=0,
                       help='Extension number')
    parser.add_argument('--options', dest='options', type=str,
                       default='',
                       help='Options')
    args = parser.parse_args()

    md5hash = md5.new(args.__str__()).hexdigest()

    make_hdr(args.inputFile, args.hdu, args.ra, args.dec, args.psize, md5hash)
    make_fits(args.inputFile, md5hash)
    print args
    make_png(md5hash, options=args.options)

