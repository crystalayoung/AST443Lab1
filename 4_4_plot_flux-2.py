# This codes are used for plotting flux vs time

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import math as mt
import collections as col
from scipy import stats
from astropy.io import fits
from xlwt import Workbook

sys.path.append('C:/My Programs/Python27/PHY517/Lab1')

def convert(ra,dec):
    ra = str(ra)
    dec = str(dec)
    ra = np.asarray(np.float32(ra.split()))
    dec = np.asarray(np.float32(dec.split()))
    ra = round((ra[0]*15+ra[1]/4.0+ra[2]/240.0)*10**8)/10**8
    dec = round((dec[0]+dec[1]/60.0+dec[2]/3600.0)*10**8)/10**8
    return ra,dec

def import_data(cat_file):
    file_cat = np.loadtxt(cat_file)
    return file_cat

def create_list(path,extenstions):
    file_list = [fn for fn in os.listdir(path) if any(fn.endswith(ext) for ext in extenstions)]
    return file_list

def refer_star(path,file_list,ra,dec):
    ref_ra = []
    ref_dec = []
    if isinstance(ra,str):
        ra,dec = convert(ra,dec)
    prompt = 'Which image do you want to use?/Number:\n'
    i = input(prompt)
    cat_file = path+'/'+file_list[i]
    cat_value = import_data(cat_file)
    x_label = np.asarray(np.asarray(np.where(abs(cat_value[:,3]-ra)<0.001))[0])
    y_label = np.asarray(np.asarray(np.where(abs(cat_value[:,4]-dec)<0.001))[0])
    index = find_equal(x_label,y_label)
    flux = cat_value[index,5]
    prompt = 'How close do you want to use?/Number:\n'
    i = input(prompt)
    index = np.asarray(np.where(abs(cat_value[:,5]-flux)<i))
    for i in index:
        ref_ra.append(cat_value[i,3])
        ref_dec.append(cat_value[i,4])
    return ref_ra,ref_dec
    
def find_equal(x_label,y_label):
    equal = 0
    for i in x_label:
        for j in y_label:
            if i == j:
                equal = i
    return equal

def table_data(path,file_list,ra,dec):
    if isinstance(ra,str):
        ra,dec = convert(ra,dec)
    flux = []
    flux_std = []
    for cat_file in file_list:
        cat_file = path+'/'+cat_file
        cat_value = import_data(cat_file)
        x_label = np.asarray(np.asarray(np.where(abs(cat_value[:,3]-ra)<0.001))[0])
        y_label = np.asarray(np.asarray(np.where(abs(cat_value[:,4]-dec)<0.001))[0])
        index = find_equal(x_label,y_label)
        flux.append(cat_value[index,5])
        flux_std.append(cat_value[index,6])
    flux = np.asarray(flux)
    flux_std = np.asarray(flux_std)
    return flux,flux_std

def norm_data(path,file_list,ra,dec):
    flux,flux_std = table_data(path,file_list,ra,dec)
    flux_aver = np.average(flux)
    norm_flux = np.divide(flux,flux_aver)
    norm_flux_std = np.divide(flux_std,flux_aver)
    return norm_flux,norm_flux_std

def weighted_mean(path,file_list,n,ra,dec):
    flux,flux_std = norm_data(path,file_list,ra[0],dec[0])
    flux_s = np.zeros_like(flux)
    flux_std_s = np.zeros_like(flux_std)
    for i in range(n):
        flux_s = np.add(flux_s,np.divide(flux,flux_std**2))
        flux_std_s = np.add(flux_std_s,np.divide(1,flux_std**2))
        if i+1 <n:
            flux,flux_std = norm_data(path,file_list,ra[i+1],dec[i+1])
    mu = np.divide(flux_s,flux_std_s)
    mu_std = np.sqrt(np.divide(1,flux_std_s))
    return mu,mu_std

def Gre2Jul(year,month,day,hour,minute,second):
    a = (14-month)/12
    y = year+4800-a
    m = month+12*a-3
    JDN = day+(153*m+2)/5+365*y+y/4-y/100+y/400-32045
    JD = JDN+(hour-12)/24.0+minute/1440.0+second/86400.0
    return JD

def Time2Sec(hour,minute,second):
    Sec = hour*3600.0+minute*60.0+second
    return Sec

def table_time(path,file_list,form):
    Res_time = []
    if form == 'Jul':
        for fits_file in file_list:
            fits_file = path+'/'+fits_file[:-4]
            hdulist = fits.open(fits_file)
            date_time = hdulist[0].header['DATE-OBS']
            date_time = date_time.split('T')
            date = date_time[0].split('-')
            time = date_time[1].split(':')
            Jul = Gre2Jul(int(date[0]),int(date[1]),int(date[2]),int(time[0]),int(time[1]),float(time[2]))
            Res_time.append(Jul)
    else:
        for fits_file in file_list:
            fits_file = path+'/'+fits_file[:-4]
            hdulist = fits.open(fits_file)
            date_time = hdulist[0].header['TIME-OBS']
            time = date_time.split(':')
            Sec = Time2Sec(int(time[0]),int(time[1]),float(time[2]))
            Res_time.append(Sec)
    return np.asarray(Res_time)
        

def bad_index(flux,time):
    plot_sketch(flux,time)
    prompt = 'Please input low-clip:\n'
    low_clip = input(prompt)
    prompt = 'Please input high-clip:\n'
    high_clip = input(prompt)
    index1 = np.where(flux<low_clip)[0]
    index2 = np.where(flux>high_clip)[0]
    index = np.concatenate((index1,index2),axis=0)
    return index

def work_data(flux,flux_std,index):
    flux = np.delete(flux,index)
    flux_std = np.delete(flux_std,index)
    return flux,flux_std

def create_dic(flux,flux_std,time):
    star = {}
    for i in range(len(time)):
        star[time[i]] = [flux[i],flux_std[i]]
    return star

def change_order(star):
    od = col.OrderedDict(sorted(star.items()))
    return od

def order_data(star):
    od = col.OrderedDict(sorted(star.items()))
    time = od.keys()
    value = np.asarray(od.values())
    flux = value[:,0]
    flux_std = value[:,1]
    return flux,flux_std,time

def work_err(flux,flux_std,mu,mu_std):
    con1 = np.divide(flux_std,mu)
    con2 = np.divide(np.multiply(flux,mu_std),mu**2)
    return np.sqrt(np.add(con1**2,con2**2))

def plot_sketch(flux,time):
    plt.plot(time,flux)
    plt.show()

def plot_flux(flux,flux_std,time,form):
    plt.figure()
    plt.errorbar(time,flux,yerr=flux_std,xerr=None,fmt='o',color='k')
    plt.axvline(x=2458019.61082,label='Ingress Time',color = 'r')
    plt.axvline(x=2458019.66199,label='Midtime',color = 'y')
    plt.axvline(x=2458019.71360,label='Egress Time',color = 'b')
    plt.xlabel(form+' Date')
    plt.ylabel('Normalized Flux')
    plt.title('WASP-93b Transition')
    plt.legend()
    plt.show()

def baseline(star,ing_time,eg_time,form):
    flux = []
    if isinstance(ing_time,str):
        if form == 'Jul':
            ing_time = ing_time.split()
            eg_time = eg_time.split()
            ing_time = Gre2Jul(int(ing_time[0]),int(ing_time[1]),int(ing_time[2]),int(ing_time[3]),int(ing_time[4]),float(ing_time[5]))
            eg_time = Gre2Jul(int(eg_time[0]),int(eg_time[1]),int(eg_time[2]),int(eg_time[3]),int(eg_time[4]),float(eg_time[5]))
        elif form == 'Sec':
            ing_time = ing_time.split()
            eg_time = eg_time.split()
            ing_time = Time2Sec(int(ing_time[0]),int(ing_time[1]),int(eg_time[2]))
            eg_time = Time2Sec(int(eg_time[0]),int(eg_time[1]),float(eg_time[3]))
    key_times = star.keys()
    for time in key_times:
        if time<=ing_time or time>=eg_time:
            flux.append(star[time])
    flux = np.asarray(flux)[:,0]
    return np.mean(flux)

def bin_data(flux,time,diff_time):
    bin_time = []
    bin_data = []
    bin_err = []
    set_time = time[-1]-diff_time
    bin_time.append(time[-1])
    bin_in_da = []
    for i in reversed(range(len(time))):
        if time[i] > set_time:
            bin_in_da.append(flux[i])
        else:
            bin_in_da = np.asarray(bin_in_da)
            aver_da = np.mean(bin_in_da)
            l = float(len(bin_in_da))
            N = l*(l-1)
            aver_er = mt.sqrt(np.sum((bin_in_da-aver_da)**2)/N)
            bin_in_da = []
            bin_in_da.append(flux[i])
            set_time -= diff_time
            bin_time.append(time[i])
            bin_data.append(aver_da)
            bin_err.append(aver_er)
    bin_in_da = np.asarray(bin_in_da)
    aver_da = np.mean(bin_in_da)
    l = float(len(bin_in_da))
    N = l*(l-1)
    aver_er = mt.sqrt(np.sum((bin_in_da-aver_da)**2)/N)
    bin_data.append(aver_da)
    bin_err.append(aver_er)
    bin_data = np.asarray(bin_data)
    bin_err = np.asarray(bin_err)
    bin_time = np.asarray(bin_time)
    index = np.where(np.isnan(bin_err))[0]
    bin_data = np.delete(bin_data,index)
    bin_err = np.delete(bin_err,index)
    bin_time = np.delete(bin_time,index)
    return bin_data,bin_err,bin_time

def save_ref(path,file_list,ra,dec,n):
    wb = Workbook()
    sheet1 = wb.add_sheet('Sheet 1')
    sheet2 = wb.add_sheet('Sheet 2')
    sheet1.write(0,0,'RA')
    sheet1.write(1,0,'dec')
    sheet1.write(2,0,'flux')
    sheet2.write(0,0,'RA')
    sheet2.write(1,0,'dec')
    sheet2.write(2,0,'flux_err')
    for i in range(n):
        sheet1.write(0,i+1,ra[i])
        sheet1.write(1,i+1,dec[i])
        sheet2.write(0,i+1,ra[i])
        sheet2.write(1,i+1,dec[i])
        flux,flux_std = table_data(path,file_list,ra[i],dec[i])
        for j in range(len(flux)):
            sheet1.write(j+2,i+1,flux[j])
            sheet2.write(j+2,i+1,flux_std[j])
    wb.save('table of flux and error for 10 reference stars.xls')

def save_tab(flux,flux_std,mu,mu_std,nflux,nflux_std,time):
    wb = Workbook()
    sheet1 = wb.add_sheet('Sheet 1')
    sheet1.write(0,0,'Time')
    sheet1.write(0,1,'Flux')
    sheet1.write(0,2,'Flux_err')
    sheet1.write(0,3,'Mu')
    sheet1.write(0,4,'Mu_err')
    sheet1.write(0,5,'Norm_flux')
    sheet1.write(0,6,'Norm_flux_err')
    for i in range(len(time)):
        sheet1.write(i+1,0,time[i])
        sheet1.write(i+1,1,flux[i])
        sheet1.write(i+1,2,flux_std[i])
        sheet1.write(i+1,3,mu[i])
        sheet1.write(i+1,4,mu_std[i])
        sheet1.write(i+1,5,nflux[i])
        sheet1.write(i+1,6,nflux_std[i])
    wb.save('table of flux, mu, norm_flux and errors.xls')
        
