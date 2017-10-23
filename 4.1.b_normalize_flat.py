#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  8 13:36:17 2017

@author: chris
"""

############################################################################
# Creates a normalized master flat field from the median combined master flat
# (Median combine done with iraf)
############################################################################

### for array operations
import numpy as np

### for operations on FITS images
from astropy.io import fits

### statistics functions
from scipy import stats

### Open median combined master flat
hdulist = fits.open('ast443_lab1_fits/lab1_science_images/master_flat.fits')

### get the data as an array
imagedata = hdulist[0].data

### change to a 1d list, for calculating mode
countvalues = imagedata.flatten()

### Find the mode
mode = stats.mode(countvalues)[0][0]

### Divide by the mode
imagedata = np.divide(imagedata,mode)

### write to a new file
newhdu = fits.PrimaryHDU(imagedata)
newhdu.writeto('ast443_lab1_fits/lab1_science_images/master_flat_normalized.fits', clobber=True)

### you open it, you close it
hdulist.close()
