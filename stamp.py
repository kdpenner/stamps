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
from aplpy import image_util
from numpy import count_nonzero
from astropy.io.fits.hdu.image import PrimaryHDU

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
        
      if count_nonzero(cutout.data) == 0:
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

    img.close()
    
  outputeps(len(coords))

def sort_src_waves(src_wave):
  return src_wave[0]

def outputeps(num_srcs):

  filter_waves = {}
  filter_waves['F435W'] = 0.435
  filter_waves['F606W'] = 0.606
  filter_waves['F814W'] = 0.814
  filter_waves['F105W'] = 1.05
  filter_waves['F125W'] = 1.25
  filter_waves['F140W'] = 1.40
  filter_waves['F160W'] = 1.60

  for src in xrange(num_srcs):
  
    files = glob.glob('output/'+str(src)+'img*.fits')
  
    fig = plt.figure(figsize = (7.5*len(files), 7.5))

    rms = None
    
    counter = 0

    if len(files) > 1:

      marg = .06

      imgs = []
      src_waves = {}

      for fileind, file in enumerate(files):
        imgs.extend(fits.open(file))
        if fileind is not 0:
          if imgs[-1].header['FILTER'] in src_waves:
            src_waves[imgs[-1].header['FILTER']].extend( \
            [imgs[-1], file])
          else:
            src_waves[imgs[-1].header['FILTER']] = \
            [filter_waves[imgs[-1].header['FILTER']], imgs[-1], file]
      sorted_src_waves = sorted(src_waves.values(),
      key = sort_src_waves)
    
      sorted_src_waves.insert(0, [0, imgs[0], files[0]])

      if len(files) >= 4:
        rgbflag = 1
      elif len(files) < 4:
        rgbflag = 0

    elif len(files) == 1:

      marg = .2

      sorted_src_waves = [[0, fits.open(files[0])[0], files[0]]]

      rgbflag = 0
      
    for imgind, wave in enumerate(sorted_src_waves):
      
      find_imgs = [each for each in wave if type(each) is PrimaryHDU]
      
      for img in find_imgs:
        f = aplpy.FITSFigure(img, figure = fig, 
        subplot = [marg+(1.-2.*marg)/(len(files)+rgbflag)*counter, .1,
        (1.-2.*marg)/(len(files)+rgbflag), .8])
      
        vmin = None
        vmax = None
        pmin = .25
        pmax = 97.

        if 'BMAJ' in img.header:
          f.add_beam()
          f.beam.set_major(img.header['BMAJ'])
          f.beam.set_minor(img.header['BMIN'])
          f.beam.set_angle(img.header['BPA'])
          f.beam.show(corner = 'top left', color = 'white', pad = 4)

        if 'FILTER' in img.header:
          f.add_label(0.2, 0.9, img.header['FILTER'], relative = True, 
          color = 'white')

        if imgind is not 0:
          img_contour = sorted_src_waves[0][1]
          if not rms:
            mask = make_source_mask(img_contour.data, snr = 2.,
            npixels = 5., dilate_size = 11.)
            mean, median, rms = sigma_clipped_stats(img_contour.data,
            sigma = 3., mask = mask)
          f.show_contour(img_contour, levels = (3.*rms, 5.*rms,
          10.*rms), colors = 'red')
          f.hide_yaxis_label()
          f.hide_ytick_labels()
          f.hide_xaxis_label()
          f.hide_xtick_labels()

        f.tick_labels.set_xformat('ddd.ddd')
        f.tick_labels.set_yformat('ddd.ddd')
        f.show_grayscale(interpolation = 'none', vmin = vmin,
        vmax = vmax, pmin = pmin, pmax = pmax)
        counter += 1

    if rgbflag:
      just_imgs = [blah[2] for blah in sorted_src_waves]
      aplpy.make_rgb_image(just_imgs[-3:], 'output/'+str(src)+'rgb.eps',
      pmin_r = pmin, pmin_g = pmin, pmin_b = pmin,
      pmax_r = pmax, pmax_g = pmax, pmax_b = pmax)
      f = aplpy.FITSFigure(just_imgs[-1], figure = fig,
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
  
  cutout(imgs, ra, dec)


if __name__ == '__main__':
  main()