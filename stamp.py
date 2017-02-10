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
from numpy import std

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

    for imgind, file in enumerate(files):

      img = fits.open(file)
      f = aplpy.FITSFigure(img[0], figure = fig, 
      subplot = [.04+.92/len(files)*imgind, .1, .92/len(files), .8])
      if 'BMAJ' in img[0].header:
        f.add_beam()
        f.beam.set_major(img[0].header['BMAJ'])
        f.beam.set_minor(img[0].header['BMIN'])
        f.beam.set_angle(img[0].header['BPA'])
        f.beam.show(corner = 'top left', color = 'white', pad = 4)
        rms = std(img[0].data)
        plt.contour(img[0].data, levels = (3.*rms, 4.*rms, 5.*rms), 
        colors = 'white')
      if 'FILTER' in img[0].header:
        f.add_label(0.2, 0.9, img[0].header['FILTER'], relative = True, 
        color = 'white')
      f.show_grayscale(interpolation = 'none')
  
    fig.canvas.draw()
    fig.savefig('output/'+str(src)+'.eps')

def main():

  args = sys.argv[1:]

  imgfnames = args[0:8]
  catfname = args[8]
  
  cat = Table.read(catfname, format = 'ascii')
  imgs = []
  for imgfname in imgfnames:
    imgs.append(fits.open(imgfname))

  ra = cat['col2']*u.degree
  dec = cat['col3']*u.degree
  
  cutout(imgs, ra, dec)


if __name__ == '__main__':
  main()