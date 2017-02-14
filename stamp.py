#!/usr/bin/env python -tt

from astropy.nddata.utils import Cutout2D
from astropy.nddata.utils import NoOverlapError
from astropy import units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import sys
from astropy.table import Table
import matplotlib
matplotlib.use('Agg')
import os
import aplpy
import matplotlib.pyplot as plt
import glob
from astropy.stats import sigma_clipped_stats
from photutils import make_source_mask

def cutout(imgs, ra, dec):

  coords = SkyCoord(ra = ra, dec = dec, frame = 'fk5')

  for j, img in enumerate(imgs):

    imgflag = 0
    try:
      imgdata = img[0].data[0,0,:,:]
    except IndexError:
      imgdata = img[0].data
      imgflag = 1

    if imgflag == 0:
      wcs = WCS(img[0].header).dropaxis(3).dropaxis(2)
    elif imgflag == 1:
      wcs = WCS(img[0].header)

    bmaj = None
    filter = None

    if 'BMAJ' in img[0].header:
      bmaj = img[0].header['BMAJ']
      bmin = img[0].header['BMIN']
      bpa = img[0].header['BPA']

    if 'FILTER' in img[0].header:
      filter = img[0].header['FILTER']
    
    if 'FILTER2' in img[0].header:
      filter = img[0].header['FILTER2']
      if filter[0] is not 'F':
        filter = img[0].header['FILTER1']

    for i, coord in enumerate(coords):

      overlapflag = 0
      try:
        cutout = Cutout2D(imgdata, coord, 10.*u.arcsec, wcs = wcs)
      except NoOverlapError:
        overlapflag = 1
        
      if overlapflag == 0:
      
        header = cutout.wcs.to_header()

        if bmaj:
          header['BMAJ'] = bmaj
          header['BMIN'] = bmin
          header['BPA'] = bpa
        if filter:
          header['FILTER'] = filter
        
        hdu = fits.PrimaryHDU(header = header, data = cutout.data)
        mkdir_err = ''
        if not os.path.exists('output/'):
          try:
            os.mkdir('output/')
          except OSError as mkdir_err:
            print mkdir_err
            sys.exit(1)
        hdu.writeto('output/'+str(i)+'img'+str(j)+'.fits')
    
  outputeps(len(coords))
    
def outputeps(num_srcs):

  for src in xrange(num_srcs):
  
    files = glob.glob('output/'+str(src)+'img*.fits')
  
    fig = plt.figure(figsize = (7.5*len(files), 7.5))

    rms = None

    if len(files) > 1:
      marg = .06
      if len(files) >= 4:
        rgbflag = 1
      elif len(files) < 4:
        rgbflag = 0
    elif len(files) == 1:
      marg = .2
      rgbflag = 0

    for imgind, file in enumerate(files):
    
      img = fits.open(file)
      f = aplpy.FITSFigure(img[0], figure = fig, 
      subplot = [marg+(1.-2.*marg)/(len(files)+rgbflag)*imgind, .1,
      (1.-2.*marg)/(len(files)+rgbflag), .8])

      if 'BMAJ' in img[0].header:
        f.add_beam()
        f.beam.set_major(img[0].header['BMAJ'])
        f.beam.set_minor(img[0].header['BMIN'])
        f.beam.set_angle(img[0].header['BPA'])
        f.beam.show(corner = 'top left', color = 'white', pad = 4)

      if 'FILTER' in img[0].header:
        f.add_label(0.2, 0.9, img[0].header['FILTER'], relative = True, 
        color = 'white')

      if imgind is not 0:
        img_contour = fits.open(files[0])
        if not rms:
          mask = make_source_mask(img_contour[0].data, snr = 2.,
          npixels = 5., dilate_size = 11.)
          mean, median, rms = sigma_clipped_stats(img_contour[0].data,
          sigma = 3., mask = mask)
        f.show_contour(img_contour[0], levels = (10.*rms, 20.*rms,
        30.*rms), colors = 'red')
        f.hide_yaxis_label()
        f.hide_ytick_labels()
        f.hide_xaxis_label()
        f.hide_xtick_labels()
 
      f.show_grayscale(interpolation = 'none')

    if rgbflag:
      aplpy.make_rgb_image(files[-3:], 'output/'+str(src)+'rgb.eps')
      f = aplpy.FITSFigure(files[-1], figure = fig,
      subplot = [marg+(1.-2.*marg)/(len(files)+rgbflag)*len(files), .1,
      (1.-2.*marg)/(len(files)+rgbflag), .8])
      f.hide_yaxis_label()
      f.hide_ytick_labels()
      f.hide_xaxis_label()
      f.hide_xtick_labels()
      f.show_rgb('output/'+str(src)+'rgb.eps')

    fig.canvas.draw()
    fig.savefig('output/'+str(src)+'.eps')

def main():

  args = sys.argv[1:]

  if not args:
    print "Usage: --radio_img file --cat catalog [--imgs file1 file2 ...]"
    sys.exit(1)
    
  if args[0] == '--radio_img':
    radioimgfname = args[1]
    del args[0:2]
  else:
    print 'Incorrect command line usage, radio img required'
    sys.exit(1)
    
  if args[0] == '--cat':
    catfname = args[1]
    del args[0:2]
  else:
    print 'Incorrect command line usage, catalog required'
    sys.exit(1)

  imgs = [fits.open(radioimgfname)]

  if args:
    imgfnames = args[1:]
    for imgfname in imgfnames:
      imgs.append(fits.open(imgfname))
      
  cat = Table.read(catfname, format = 'ascii')

# radio img will always be first

  ra = cat['col2']*u.degree
  dec = cat['col3']*u.degree
  ra = ra[0:20]
  dec = dec[0:20]
  
  cutout(imgs, ra, dec)


if __name__ == '__main__':
  main()