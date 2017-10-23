#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 23:19:24 2017

@author: chris
"""
##############################################################################
#   Identifies hot/dead pixels in a fit and creates a bad pixel map.
#   The new fit is put in the lab 1 science images folder(you can change this) 
#   and named 'badmap.fits'
##############################################################################


### for array operations
import numpy as np
### for sigma clipping
from scipy.stats import sigmaclip
### for operations on FITS images
from astropy.io import fits

### a and b are 1d lists. Identifies hot pixels
def hot(a, b):
    # return an array where a 1 represents an element which was not removed 
    # and a 0 is an element which was removed
    i = np.where(np.in1d(a, b))[0]
    hot = np.zeros_like(a)
    hot[i] = 1
    return hot

### f and g are 1d lists. Identifies dead pixels
def dead(f,g):
    # return an array where a 1 represents an element which was not removed 
    # and a 0 is an element which was removed
    i = np.where(np.in1d(f, g))[0]
    index = np.zeros_like(f)
    index[i] = 1
    return dead

### i and j are 1d lists. Puts together the hot and dead pix.
def hot_dead(i,j):
    # return an array where a 1 represents an element which was not removed 
    # and a 0 is an element which was removed
    a = np.where(j==0)[0]
    b = np.where(i==0)[0]
    hot_dead = np.ones_like(i)
    hot_dead[a] = 0
    hot_dead[b] = 0
    return hot_dead

### open the master dark
hdulist = fits.open('ast443_lab1_fits/lab1_science_images/master_dark.fits')
### open the master flat
hdulist2 = fits.open('ast443_lab1_fits/lab1_science_images/master_flat_normalized.fits')

### get the image data as an array
imagedata = hdulist[0].data
imagedata2 = hdulist2[0].data

### convert the 2d array into a 1d list
countvalues = imagedata.flatten()
countvalues2 = imagedata2.flatten()

### sigma clip to set hot pix threshold, fact=sigma
fact = 5.0
c, low, upp = sigmaclip(countvalues, fact, fact)

### get values for pixels that aren't dead
countvalues3 = countvalues2[countvalues2 > 0.8]

### define the 1d lists for hot(a,b)
a=countvalues[:]
b=c[:]

### Call hot
i = hot(a,b)

### define the 1d lists for dead(f,g)
f=countvalues[:]
g=countvalues3[:]

### Call dead
j = dead(f,g)

### Call hot_dead
k = hot_dead(i,j)

### Get the flattened data back to an array
k = k.reshape(imagedata[0,:].size,-1)

### write into a new fits file
newhdu = fits.PrimaryHDU(k)
newhdu.writeto('ast443_lab1_fits/lab1_science_images/badmap.fits', clobber=True)

### you open it, you close it
hdulist.close() 
hdulist2.close() 
newhdu.close()
