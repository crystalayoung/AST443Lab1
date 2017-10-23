#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 07:53:12 2017

@author: crystalyoung
"""
############################################################################
#   Takes in all the science images and uses the master dark and
#   master flat to correct the science images.
#   Be sure to check pathways and file names so Python knows where
#   your files are and what to name the new files.
############################################################################

### I don't think we need numpy for this
import numpy as np

### For operations on fits
from astropy.io import fits

### For opening multiple fits
import glob

### Arrays for science images
image_list = []
image_data = []

### Open the science images
### (change the directory before the file name below if needed:)
for filename in glob.glob('wasp_93_transit_10s.*.FIT'):
    im=fits.open(filename)
    image_list.append(im)

### Get the data from the images
for i in image_list:
    data = i[0].data
    image_data.append(data)
    
### Open the master dark image -- check pathways
imagedark = fits.open('master_dark.fits')
imagedarkdata = imagedark[0].data
### Open the master flat image
imageflat = fits.open('master_flat.fits')
imageflatdata = imageflat[0].data

### write out the corrected fits
for i in range(len(image_data)):
    newhdu = fits.PrimaryHDU((image_data[i]-imagedarkdata)/imageflatdata)
    newhdu.writeto("new_wasp_93_transit_10s" + str(i) + ".fits", overwrite=True)
    
### Close all open fits
imagedark.close()
imageflat.close()
im.close()