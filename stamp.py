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
from numpy import count_nonzero
from astropy.io.fits.hdu.image import PrimaryHDU
import pdb

def cutout(imgs, ra, dec, targname):

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

    if set(['CRVAL3', 'CUNIT3']).issubset(set(img[0].header)):
      blah = img[0].header['CRVAL3']*u.Unit(img[0].header['CUNIT3'])
      blah = blah.to(u.GHz)
      filter = "{0:0.2f}".format(blah)+" ({0:0.2f})".format(blah.to(u.cm,
      equivalencies = u.spectral()))

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
          
        header['TARGNAME'] = targname[i]
        
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
  
  mkdir_err = ''
  if not os.path.exists('output/no_hst_counterpart/'):
    try:
      os.mkdir('output/no_hst_counterpart/')
    except OSError as mkdir_err:
      print mkdir_err
      sys.exit(1)

  for src in xrange(num_srcs):
  
    files = glob.glob('output/'+str(src)+'img*.fits')

    rms = None
    
    counter = 0

    if len(files) > 1:

      imgs = []
      src_waves = {}

      for fileind, file in enumerate(files):
        imgs.extend(fits.open(file))
        if fileind is not 0:
          if 'FILTER' in imgs[-1].header:
            if imgs[-1].header['FILTER'] in src_waves:
              src_waves[imgs[-1].header['FILTER']].extend( \
              [imgs[-1], file])
            else:
              try:
                src_waves[imgs[-1].header['FILTER']] = \
                [filter_waves[imgs[-1].header['FILTER']], imgs[-1], file]
              except KeyError:
                src_waves['NOTAVAIL'] = [0, imgs[-1], file]
          else:
            src_waves['NOTAVAIL'] = [0, imgs[-1], file]
      sorted_src_waves = sorted(src_waves.values(),
      key = sort_src_waves)
    
      sorted_src_waves.insert(0, [0, imgs[0], files[0]])

      key_array = src_waves.keys()
      if 'NOTAVAIL' in key_array:
        key_array.remove('NOTAVAIL')

      if len(key_array) >= 3:
        rgbflag = 1
        append = ''
      elif len(key_array) < 3 and len(key_array) > 0:
        rgbflag = 0
        append = ''
      elif len(key_array) == 0:
        rgbflag = 0
        append = 'no_hst_counterpart/'

      marg_y = .15
      marg_in = 7.5/(1.-2.*marg_y)*marg_y
      width = 7.5*(len(files)+rgbflag)+2.*marg_in
      marg_x = marg_in/width

    elif len(files) == 1:

      marg_x = .15
      marg_y = .15

      sorted_src_waves = [[0, fits.open(files[0])[0], files[0]]]

      rgbflag = 0
      
      append = 'no_hst_counterpart/'
    
    fig = plt.figure(figsize = (7.5*(len(files)+rgbflag)/(1.-2.*marg_x), 7.5/(1.-2.*marg_y)))

    for imgind, wave in enumerate(sorted_src_waves):
      
      find_imgs = fits.HDUList([each for each in wave if type(each) is PrimaryHDU])
      
      for img in find_imgs:

        f = aplpy.FITSFigure(img, figure = fig, 
        subplot = [marg_x+(1.-2.*marg_x)/(len(files)+rgbflag)*counter, marg_y,
        (1.-2.*marg_x)/(len(files)+rgbflag), 1.-2.*marg_y])

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
          
        titleadd = ''
        if 'FILTER' in img.header:
          titleadd += img.header['FILTER']
          
        if imgind is not 0:
          img_contour = sorted_src_waves[0][1]
          f.recenter(ra_cen, dec_cen, radius = (5.*u.arcsec).to(u.deg).value)
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
          f.set_title(titleadd)
        else:
          f.set_title(img.header['TARGNAME']+'\n'+titleadd)
          indcen = len(img.data)/2
          wcs = WCS(img.header)
          ra_cen, dec_cen = wcs.all_pix2world(indcen, indcen, 0)
          f.recenter(ra_cen, dec_cen, radius = (5.*u.arcsec).to(u.deg).value)
          

        f.tick_labels.set_xformat('ddd.ddd')
        f.tick_labels.set_yformat('ddd.ddd')
        f.show_grayscale(interpolation = 'none', vmin = vmin,
        vmax = vmax, pmin = pmin, pmax = pmax)
        counter += 1
        del img.data
        del img
      
      find_imgs.close()


    if rgbflag:
      just_imgs = [blah[2] for blah in sorted_src_waves if blah[0] is not 0]
      just_imgs_rgb = just_imgs[-3:]
      just_imgs_rgb.reverse()
      aplpy.make_rgb_image(just_imgs_rgb,
      'output/'+append+str(src)+'rgb.eps',
      pmin_r = pmin, pmin_g = pmin, pmin_b = pmin,
      pmax_r = pmax, pmax_g = pmax, pmax_b = pmax)
      f = aplpy.FITSFigure(just_imgs[-1], figure = fig,
      subplot = [marg_x+(1.-2.*marg_x)/(len(files)+rgbflag)*len(files), marg_y,
      (1.-2.*marg_x)/(len(files)+rgbflag), 1.-2.*marg_y])
      f.recenter(ra_cen, dec_cen, radius = (5.*u.arcsec).to(u.deg).value)
      f.hide_yaxis_label()
      f.hide_ytick_labels()
      f.hide_xaxis_label()
      f.hide_xtick_labels()
      f.show_rgb('output/'+append+str(src)+'rgb.eps')

    fig.canvas.draw()
    fig.savefig('output/'+append+str(src)+'.png')
    plt.close(fig)

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

  ra = cat['col2'][0:100]*u.degree
  dec = cat['col3'][0:100]*u.degree
  targname = cat['col1'][0:100]
  
  cutout(imgs, ra, dec, targname)
  outputeps(len(ra))


if __name__ == '__main__':
  main()