#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 21 23:39:01 2017

@author: chris (crystal and Hu contributed)
"""

##############################################################################
# This was used for ast443 lab1:exoplanet transit. Given the ra and dec of 10
# reference stars, along with the target star, this program can adjust for
# differences in seeing and give the transit depth. From the tranist depth,
# the ratio between the the radius of the star and the radius of the planet
# can be calculated. Rp/Rs = sqrt(depth)
##############################################################################

### for array operations
import numpy as np
### load multiple files using in order sorted(glob.glob(*.fit))
import glob
### used to get the file names ordered by number
import re
### for utc to jd conversion
from astropy.time import Time
### for header info
import pyfits
### for plotting
import matplotlib.pyplot as plt
### used for creating tables and saving them
from astropy.table import Table, Column
from astropy.io import ascii

###Used to view an entire array without ...
np.set_printoptions(threshold=np.nan)

### Get ra and dec for 10 reference stars
### ra and dec obtained from .cat files
ra0, ra1, ra2, ra3, ra4, ra5, ra6, ra7, ra8, ra9 = [ 9.6202277,  9.4864205,  9.8800711,  9.3715302,  9.2912587,
        9.2740893,  9.5327732,  9.3045929,  9.5796475,  9.8891627]

dec0, dec1, dec2, dec3, dec4, dec5, dec6, dec7, dec8, dec9 = [ 51.176707 ,  51.1883238,  51.2919868,  51.3788731,  51.4390164,
        51.4464428,  51.5089537,  51.5157426,  51.5126596,  51.5392065]

### ra and dec of the target star
target_ra = 9.45833  
target_dec =   51.2889

### used to sort in files by image number
def keyFunc(afilename):
    nondigits = re.compile("\D")
    return int(nondigits.sub("",afilename))

### Used to find where to values are equal
def find_equal(x_label,y_label):
    equal = 0
    for i in x_label:
        for j in y_label:
            if i == j:
                equal = i
    return equal

### Get the normalized flux and flux_err(flux_std)
def flux_norm(ra,dec):
    flux = []
    flux_list = []
    flux_std = []
    flux_std_list = []
    ### Used sorted glob to take files in order
    for filename in sorted(glob.glob('ast443_lab1_fits/lab1_science_images/cat_files/*.new.cat'), key=keyFunc):
        file_cat = np.loadtxt(filename)
        ### find the index
        x_label = np.asarray(np.asarray(np.where(abs(file_cat[:,3]-ra)<0.001))[0])
        y_label = np.asarray(np.asarray(np.where(abs(file_cat[:,4]-dec)<0.001))[0])
        index = find_equal(x_label,y_label)
        ### find flux from index and put into a list
        flux = (file_cat[index,5])
        flux_std = (file_cat[index,6])
        flux_list.append(flux)
        flux_std_list.append(flux_std)
    ### normalize the data to the mean
    flux_norm = np.divide(flux_list,np.mean(flux_list))
    err_norm = np.divide(flux_std_list,np.mean(flux_list))
    return flux_norm , err_norm

### Get the flux and flux_err of the target star
def target_flux(ra,dec):
    flux = []
    flux_list = []
    flux_std = []
    flux_std_list = []
    ### Used sorted glob to take files in order
    for filename in sorted(glob.glob('ast443_lab1_fits/lab1_science_images/cat_files/*.new.cat'), key=keyFunc):
        file_cat = np.loadtxt(filename)
        ### find the index
        x_label = np.asarray(np.asarray(np.where(abs(file_cat[:,3]-ra)<0.001))[0])
        y_label = np.asarray(np.asarray(np.where(abs(file_cat[:,4]-dec)<0.001))[0])
        index = find_equal(x_label,y_label)
        ### find flux from index
        flux = (file_cat[index,5])
        flux_std = (file_cat[index,6])
        flux_list.append(flux)
        flux_std_list.append(flux_std)
    flux_array = np.asarray(flux_list)
    err_array = np.asarray(flux_std_list)
    return flux_array , err_array

 
### normalized flux(and errror) for all 10 reference stars
a0, e0 = flux_norm(ra0,dec0)
a1, e1 = flux_norm(ra1,dec1)
a2, e2 = flux_norm(ra2,dec2)
a3, e3 = flux_norm(ra3,dec3)
a4, e4 = flux_norm(ra4,dec4)
a5, e5 = flux_norm(ra5,dec5)
a6, e6 = flux_norm(ra6,dec6)
a7, e7 = flux_norm(ra7,dec7)
a8, e8 = flux_norm(ra8,dec8)
a9, e9 = flux_norm(ra9,dec9)
### flux and error for target star (not normalized)
target_flux, target_flux_err = target_flux(target_ra,target_dec)

### square the error lists. we will use this later to calculate the --
### -- weighted mean and its error
b0 = [x**2 for x in e0]
b1 = [x**2 for x in e1]
b2 = [x**2 for x in e2]
b3 = [x**2 for x in e3]
b4 = [x**2 for x in e4]
b5 = [x**2 for x in e5]
b6 = [x**2 for x in e6]
b7 = [x**2 for x in e7]
b8 = [x**2 for x in e8]
b9 = [x**2 for x in e9]

### convert the lists into arrays
a0 = np.asarray(a0)
a1 = np.asarray(a1)
a2 = np.asarray(a2)
a3 = np.asarray(a3)
a4 = np.asarray(a4)
a5 = np.asarray(a5)
a6 = np.asarray(a6)
a7 = np.asarray(a7)
a8 = np.asarray(a8)
a9 = np.asarray(a9)
b0 = np.asarray(b0)
b1 = np.asarray(b1)
b2 = np.asarray(b2)
b3 = np.asarray(b3)
b4 = np.asarray(b4)
b5 = np.asarray(b5)
b6 = np.asarray(b6)
b7 = np.asarray(b7)
b8 = np.asarray(b8)
b9 = np.asarray(b9)

### calculate the weighted mean
mun = (a0/b0)+(a1/b1)+(a2/b2)+(a3/b3)+(a4/b4)+(a5/b5)+(a6/b6)+(a7/b7)+(a8/b8)+(a9/b9)
mud = (1/b0)+(1/b1)+(1/b2)+(1/b3)+(1/b4)+(1/b5)+(1/b6)+(1/b7)+(1/b8)+(1/b9)
mu = mun/mud
### calculate the error on the weighted mean
sig = 1/mud
sigma = np.sqrt(sig)

### Calc the the weather corrected science images with its error
corrected_flux = target_flux/mu
mu2 = np.square(mu)
err1 = target_flux_err/mu
err2 = sigma*target_flux/mu2
err12 = np.square(err1)
err22 = np.square(err2)
err0 = np.add(err12, err22)
err = np.sqrt(err0)

### Make a list of the non-transit images and find the average
### This is the baseline flux
x = corrected_flux[0:39]
y = corrected_flux[530:734]
z = np.append(x,y)
baseline = np.mean(z)

### Calculate the error in the baseline flux
#a = err[0:40]
#b = err[530:734]
#c = np.append(a,b)
#c2 = np.square(c)
#s = np.sum(c2)
#sr = np.sqrt(s)
#baseline_err = sr/len(c)

### Normalize all the science images(and error) to the baseline flux
flux_norm = corrected_flux/baseline
err_norm = err/baseline

### Put the data into bins where the flux equals the avergae flux for that bin
flux_norm1 = np.mean(flux_norm[0:30])
flux_norm2 = np.mean(flux_norm[30:60])
flux_norm3 = np.mean(flux_norm[60:90])
flux_norm4 = np.mean(flux_norm[90:120])
flux_norm5 = np.mean(flux_norm[120:150])
flux_norm6 = np.mean(flux_norm[150:180])
flux_norm7 = np.mean(flux_norm[180:210])
flux_norm8 = np.mean(flux_norm[210:230])
flux_norm9 = np.mean(flux_norm[230:250])
flux_norm10 = np.mean(flux_norm[250:270])
flux_norm11 = np.mean(flux_norm[270:290])
flux_norm12 = np.mean(flux_norm[290:310])
flux_norm13 = np.mean(flux_norm[310:330])
flux_norm14 = np.mean(flux_norm[330:350])
flux_norm15 = np.mean(flux_norm[350:370])
flux_norm16 = np.mean(flux_norm[370:390])
flux_norm17 = np.mean(flux_norm[390:410])
flux_norm18 = np.mean(flux_norm[410:430])
flux_norm19 = np.mean(flux_norm[430:450])
flux_norm20 = np.mean(flux_norm[450:470])
flux_norm21 = np.mean(flux_norm[470:490])
flux_norm22 = np.mean(flux_norm[490:510])
flux_norm23 = np.mean(flux_norm[510:530])
flux_norm24 = np.mean(flux_norm[530:550])
flux_norm25 = np.mean(flux_norm[550:570])
flux_norm26 = np.mean(flux_norm[570:590])
flux_norm27 = np.mean(flux_norm[590:610])
flux_norm28 = np.mean(flux_norm[610:630])
flux_norm29 = np.mean(flux_norm[630:650])
flux_norm30 = np.mean(flux_norm[650:670])
flux_norm31 = np.mean(flux_norm[670:690])
flux_norm32 = np.mean(flux_norm[690:710])
flux_norm33 = np.mean(flux_norm[710:735])
### Put each bin into an array
flux_norm_groups = [flux_norm1, flux_norm2, flux_norm3, flux_norm4, flux_norm5, flux_norm6, flux_norm7, flux_norm8, flux_norm9, flux_norm10, flux_norm11, flux_norm12, flux_norm13, flux_norm14, flux_norm15, flux_norm16, flux_norm17, flux_norm18, flux_norm19, flux_norm20, flux_norm21, flux_norm22, flux_norm23, flux_norm24, flux_norm25, flux_norm26, flux_norm27, flux_norm28, flux_norm29, flux_norm30, flux_norm31, flux_norm32, flux_norm33]

### Do the same for the error
err_norm1 = np.mean(err_norm[0:30])
err_norm2 = np.mean(err_norm[30:60])
err_norm3 = np.mean(err_norm[60:90])
err_norm4 = np.mean(err_norm[90:120])
err_norm5 = np.mean(err_norm[120:150])
err_norm6 = np.mean(err_norm[150:180])
err_norm7 = np.mean(err_norm[180:210])
err_norm8 = np.mean(err_norm[210:230])
err_norm9 = np.mean(err_norm[230:250])
err_norm10 = np.mean(err_norm[250:270])
err_norm11 = np.mean(err_norm[270:290])
err_norm12 = np.mean(err_norm[290:310])
err_norm13 = np.mean(err_norm[310:330])
err_norm14 = np.mean(err_norm[330:350])
err_norm15 = np.mean(err_norm[350:370])
err_norm16 = np.mean(err_norm[370:390])
err_norm17 = np.mean(err_norm[390:410])
err_norm18 = np.mean(err_norm[410:430])
err_norm19 = np.mean(err_norm[430:450])
err_norm20 = np.mean(err_norm[450:470])
err_norm21 = np.mean(err_norm[470:490])
err_norm22 = np.mean(err_norm[490:510])
err_norm23 = np.mean(err_norm[510:530])
err_norm24 = np.mean(err_norm[530:550])
err_norm25 = np.mean(err_norm[550:570])
err_norm26 = np.mean(err_norm[570:590])
err_norm27 = np.mean(err_norm[590:610])
err_norm28 = np.mean(err_norm[610:630])
err_norm29 = np.mean(err_norm[630:650])
err_norm30 = np.mean(err_norm[650:670])
err_norm31 = np.mean(err_norm[670:690])
err_norm32 = np.mean(err_norm[690:710])
err_norm33 = np.mean(err_norm[710:735])
### Put each bin into an array
err_norm_groups = [err_norm1, err_norm2, err_norm3, err_norm4, err_norm5, err_norm6, err_norm7, err_norm8, err_norm9, err_norm10, err_norm11, err_norm12, err_norm13, err_norm14, err_norm15, err_norm16, err_norm17, err_norm18, err_norm19, err_norm20, err_norm21, err_norm22, err_norm23, err_norm24, err_norm25, err_norm26, err_norm27, err_norm28, err_norm29, err_norm30, err_norm31, err_norm32, err_norm33]


###########################################################################
# Gets UTC from header of a fit(s) and converts to JD
########################################################################### 

### set up a 1d lists to store the times in
OBStime = []

### loop over all files (in order)
for fitsName in sorted(glob.glob('ast443_lab1_fits/lab1_science_images/corrected_science_images/*.fits'),key=keyFunc):
    ### get header info
    hdulist = pyfits.open(fitsName)
    header = hdulist[0].header
    ### pull out UTC
    datetime = header['DATE-OBS']
    ### tell astropy the format of the date and time
    t = Time(datetime, format='isot', scale='utc')
    ### convert to JD
    julian = t.jd
    ### put the looped data into a 1d lists
    OBStime.append(julian)
    ### close the hdu after each loop
    hdulist.close()

### Put the data into bins where the time equals the avg time for that bin.
### This will put the data point in the middle of the bin
OBS1 = np.mean(OBStime[0:30])
OBS2 = np.mean(OBStime[30:60])
OBS3 = np.mean(OBStime[60:90])
OBS4 = np.mean(OBStime[90:120])
OBS5 = np.mean(OBStime[120:150])
OBS6 = np.mean(OBStime[150:180])
OBS7 = np.mean(OBStime[180:210])
OBS8 = np.mean(OBStime[210:230])
OBS9 = np.mean(OBStime[230:250])
OBS10 = np.mean(OBStime[250:270])
OBS11 = np.mean(OBStime[270:290])
OBS12 = np.mean(OBStime[290:310])
OBS13 = np.mean(OBStime[310:330])
OBS14 = np.mean(OBStime[330:350])
OBS15 = np.mean(OBStime[350:370])
OBS16 = np.mean(OBStime[370:390])
OBS17 = np.mean(OBStime[390:410])
OBS18 = np.mean(OBStime[410:430])
OBS19 = np.mean(OBStime[430:450])
OBS20 = np.mean(OBStime[450:470])
OBS21 = np.mean(OBStime[470:490])
OBS22 = np.mean(OBStime[490:510])
OBS23 = np.mean(OBStime[510:530])
OBS24 = np.mean(OBStime[530:550])
OBS25 = np.mean(OBStime[550:570])
OBS26 = np.mean(OBStime[570:590])
OBS27 = np.mean(OBStime[590:610])
OBS28 = np.mean(OBStime[610:630])
OBS29 = np.mean(OBStime[630:650])
OBS30 = np.mean(OBStime[650:670])
OBS31 = np.mean(OBStime[670:690])
OBS32 = np.mean(OBStime[690:710])
OBS33 = np.mean(OBStime[710:735])
### Put each bin into an array
OBStime_groups = [OBS1, OBS2, OBS3, OBS4, OBS5, OBS6, OBS7, OBS8, OBS9, OBS10, OBS11, OBS12, OBS13, OBS14, OBS15, OBS16, OBS17, OBS18, OBS19, OBS20, OBS21, OBS22, OBS23, OBS24, OBS25, OBS26, OBS27, OBS28, OBS29, OBS30, OBS31, OBS32, OBS33]

### Calculate the transit depth and the planet star ratio
depth = 1 - np.min(flux_norm_groups)
planet_star_ratio = np.sqrt(depth)

##############################################################################
# This is for plotting and making tables. Do not delete a plot, just comment
# it out. For the table, set clobber=True if you want to overwrite the file.
##############################################################################

### Table of JD, flux, flux_err
t = Table([OBStime, target_flux, target_flux_err], names=('JD','flux','flux_err'))
ascii.write(t, 'values2.txt')

### Plot all 10 normalized ref stars on the same plot
plt.errorbar(OBStime, a0, xerr=None, yerr=e0, fmt=',', color='red', ecolor='red')
plt.errorbar(OBStime, a1+1, xerr=None, yerr=e1, fmt=',', color='blue', ecolor='blue')
plt.errorbar(OBStime, a2+2, xerr=None, yerr=e2, fmt=',', color='green', ecolor='green')
plt.errorbar(OBStime, a3+3, xerr=None, yerr=e3, fmt=',', color='magenta', ecolor='magenta')
plt.errorbar(OBStime, a4+4, xerr=None, yerr=e4, fmt=',', color='cyan', ecolor='cyan')
plt.errorbar(OBStime, a5+5, xerr=None, yerr=e5, fmt=',', color='purple',ecolor='purple')
plt.errorbar(OBStime, a6+6, xerr=None, yerr=e6, fmt=',', color='red', ecolor='red')
plt.errorbar(OBStime, a7+7, xerr=None, yerr=e7, fmt=',', color='blue', ecolor='blue')
plt.errorbar(OBStime, a8+8, xerr=None, yerr=e8, fmt=',', color='magenta', ecolor='magenta')
plt.errorbar(OBStime, a9+9, xerr=None, yerr=e9, fmt=',', color='green', ecolor='green')
plt.title('Normalized Reference Star Flux vs Time')
plt.ylabel('Normalized Flux (Normalized to integer values)')
plt.xlabel('Julian Date')
y=[0,1,2,3,4,5,6,7,8,9,10]
plt.yticks(y)
plt.show()

## Plot flux vs time for target star
plt.errorbar(OBStime, corrected_flux, xerr=None, yerr=err, fmt=',', color='red', ecolor='black',elinewidth=0.5, capsize=5, capthick=0.5)
plt.title('Transit WASP-93 b (Corrected with avg of 10 reference stars)')
plt.xlabel('Julian Date')
plt.ylabel('Flux (Counts)')
plt.show()

### Plot normalized flux vs time for target star
plt.errorbar(OBStime, flux_norm, xerr=None, yerr=err_norm, fmt=',', color='red', ecolor='black',elinewidth=0.5, capsize=5, capthick=0.5)
plt.title('Normalized Transit WASP-93 b')
plt.xlabel('Julian Date')
plt.ylabel('Relative Flux')
plt.show()

### Plot binned, normalized flux vs time for target star 
plt.errorbar(OBStime_groups, flux_norm_groups, xerr=None, yerr=err_norm_groups, fmt=',', color='red', ecolor='black',elinewidth=0.5, capsize=5, capthick=0.5)
plt.title('Normalized Transit WASP-93 b')
plt.xlabel('Julian Date')
plt.ylabel('Relative Flux')
plt.show()

