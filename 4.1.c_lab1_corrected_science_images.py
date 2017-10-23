#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 07:53:12 2017

@author: crystalyoung
"""

import numpy as np
from astropy.io import fits
import glob

image_list = []
image_data = []

# Open the science images (change the directory before the file name below if needed:)
for filename in glob.glob('wasp_93_transit_10s.*.FIT'):
    im=fits.open(filename)
    image_list.append(im)

# Open the master dark image -- need to change below!!
##### imagedark = fits.open('master_dark.fits')
##### imagedarkdata = imagedark[0].data
# Open the master flat image
##### imageflat = fits.open('master_flat.fits')
##### imageflatdata = imageflat[0].data

# Get the data from the images
for i in image_list:
    data = i[0].data
    image_data.append(data)

for i in range(len(image_data)):
    newhdu = fits.PrimaryHDU(image_data[i]-imagedarkdata-imageflatdata)
    newhdu.writeto("new_wasp_93_transit_10s" + str(i) + ".fits", overwrite=True)