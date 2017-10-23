#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 15:30:35 2017

@author: crystalyoung
"""

#subtract the master bias from all of the dark frames in order to get the dark current

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

differences_5 = []
differences_10 = []
median5 = []
median10 = []

bias = fits.open('masterbias.fits')
image1 = fits.open('10deg_lab0.00000154.DARK_10s.FIT')
image2 = fits.open('10deg_lab0.00000155.DARK_50s.FIT')
image3 = fits.open('10deg_lab0.00000156.DARK_90s.FIT')
image4 = fits.open('10deg_lab0.00000157.DARK_130s.FIT')
image5 = fits.open('10deg_lab0.00000158.DARK_170s.FIT')
image6 = fits.open('10deg_lab0.00000159.DARK_220s.FIT')
image7 = fits.open('10deg_lab0.00000160.DARK_270s.FIT')
image8 = fits.open('10deg_lab0.00000161.DARK_300s.FIT')
image9 = fits.open('5deg_lab0.00000126.DARK_10s.FIT')
image10 = fits.open('5deg_lab0.00000127.DARK_50s.FIT')
image11 = fits.open('5deg_lab0.00000129.DARK_90s.FIT')
image12 = fits.open('5deg_lab0.00000130.DARK_170s.FIT')
image13 = fits.open('5deg_lab0.00000131.DARK_220s.FIT')
image14 = fits.open('5deg_lab0.00000132.DARK_270s.FIT')
image15 = fits.open('5deg_lab0.00000133.DARK_300s.FIT')

biasdata = bias[0].data
imagedata1 = image1[0].data
imagedata2 = image2[0].data
imagedata3 = image3[0].data
imagedata4 = image4[0].data
imagedata5 = image5[0].data
imagedata6 = image6[0].data
imagedata7 = image7[0].data
imagedata8 = image8[0].data
imagedata9 = image9[0].data
imagedata10 = image10[0].data
imagedata11 = image11[0].data
imagedata12 = image12[0].data
imagedata13 = image13[0].data
imagedata14 = image14[0].data
imagedata15 = image15[0].data

darkdata = [imagedata1, imagedata2, imagedata3, imagedata4, imagedata5, imagedata6, imagedata7, imagedata8, imagedata9, imagedata10, imagedata11, imagedata12, imagedata13, imagedata14, imagedata15]
darkdata5deg = [imagedata9, imagedata10, imagedata11, imagedata12, imagedata13, imagedata14, imagedata15]
darkdata10deg = [imagedata1, imagedata2, imagedata3, imagedata4, imagedata5, imagedata6, imagedata7, imagedata8]

for i in range(len(darkdata5deg)):
    newhdu = fits.PrimaryHDU(darkdata5deg[i]-biasdata)
    newhdu.writeto("5_deg_new_darks" + str(i) + ".fits", overwrite=True)
    
for i in range(len(darkdata10deg)):
    newhdu = fits.PrimaryHDU(darkdata10deg[i]-biasdata)
    newhdu.writeto("10_deg_new_darks" + str(i) + ".fits", overwrite=True)

imagedata26 = (fits.open('5_deg_new_darks0.fits'))[0].data
imagedata27 = (fits.open('5_deg_new_darks1.fits'))[0].data
imagedata28 = (fits.open('5_deg_new_darks2.fits'))[0].data
imagedata29 = (fits.open('5_deg_new_darks3.fits'))[0].data
imagedata30 = (fits.open('5_deg_new_darks4.fits'))[0].data
imagedata31 = (fits.open('5_deg_new_darks5.fits'))[0].data
imagedata32 = (fits.open('5_deg_new_darks6.fits'))[0].data
imagedata33 = (fits.open('10_deg_new_darks0.fits'))[0].data
imagedata34 = (fits.open('10_deg_new_darks1.fits'))[0].data
imagedata35 = (fits.open('10_deg_new_darks2.fits'))[0].data
imagedata36 = (fits.open('10_deg_new_darks3.fits'))[0].data
imagedata37 = (fits.open('10_deg_new_darks4.fits'))[0].data
imagedata38 = (fits.open('10_deg_new_darks5.fits'))[0].data
imagedata39 = (fits.open('10_deg_new_darks6.fits'))[0].data
imagedata40 = (fits.open('10_deg_new_darks7.fits'))[0].data

new_images = [imagedata33, imagedata34, imagedata35, imagedata36, imagedata37, imagedata38, imagedata39, imagedata26, imagedata27, imagedata28, imagedata29, imagedata30, imagedata31, imagedata32]
new_images_5 = [imagedata26, imagedata27, imagedata28, imagedata29, imagedata30, imagedata31, imagedata32]
new_images_10 = [imagedata33, imagedata34, imagedata35, imagedata36, imagedata37, imagedata38, imagedata39, imagedata40]

for i in range(len(new_images_5)):
    median_5 = np.median(new_images_5[i])
    median5.append(median_5)
    avg = np.mean(new_images_5[i])
    avg = np.mean(avg)
    avg = np.mean(avg)
    diff_avg_median_5 = abs(avg - median_5)
    differences_5.append((diff_avg_median_5))

for i in range(len(new_images_10)):
    median_10 = np.median(new_images_10[i])
    median10.append(median_10)
    avg = np.mean(new_images_10[i])
    avg = np.mean(avg)
    avg = np.mean(avg)
    diff_avg_median_10 = abs(avg - median_10)
    differences_10.append((diff_avg_median_10))


exposuretime5 = [10, 50, 90, 170, 220, 270, 300]
exposuretime10 = [10, 50, 90, 130, 170, 220, 270, 300]

# Error in the measurement = differences_5 and differences_10

# Plot 5s
plt.plot(exposuretime5, median5, 'ro', color='red')
plt.title('5 second and 10 second darks')
plt.errorbar(exposuretime5, median5, xerr=None, yerr=differences_5, fmt='o')

# Plot 10s
plt.plot(exposuretime10, median10, 'ro', color='blue')
plt.errorbar(exposuretime10, median10, xerr=None, yerr=differences_10, fmt='o')
plt.show