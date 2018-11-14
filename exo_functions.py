# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 10:53:42 2018

@author: jesper
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Get current time
import datetime
now = datetime.datetime.now()
timestamp = str(now.year) + '-' + str(now.month) + '-' + str(now.day) + '_' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second)

#%% Import function

def import_tess_fits(path,print_fig=False):
    hdulist = fits.open('data/'+ path)
    hdu = hdulist[1]
    
    hdu_data = hdu.data
    
    time = []
    sap_flux = []
    pdc_flux = []
    
    for i in range(len(hdu_data)):
        if np.isnan(hdu_data[i][3])==False:
            time.append(hdu_data[i][0])
            sap_flux.append(hdu_data[i][3])
            pdc_flux.append(hdu_data[i][7])
        
    time =      np.array(time)
    sap_flux =  np.array(sap_flux)
    pdc_flux =  np.array(pdc_flux)
    
    if print_fig==True:    
        plt.figure()
        plt.plot(time, sap_flux, ',')
        plt.title(path)
        plt.xlabel('time [days]')
        plt.ylabel('flux')
    
    return time, sap_flux

#%%
def normer_fluxes(times,fluxes,intervals,cutoff = 0.98,print_fig=False,save_fig=False):
    med_fluxes = []
    for i in range(len(fluxes)):
        med_flux = []
        for j in range(len(fluxes[i])):
            med_flux.append(np.median(fluxes[i][np.where(np.logical_and(times[i] > times[i][j]-intervals[i],
            times[i] < times[i][j]+intervals[i]))]))
        med_flux = np.array(med_flux)    
        med_fluxes.append(med_flux)
    
    med_fluxes=np.array(med_fluxes)
    
    norm_fluxes = []
    norm_times = []
    for i in range(len(fluxes)):
        norm_flux = fluxes[i]/med_fluxes[i]
        norm_time = times[i][np.where(norm_flux > cutoff)]
        norm_flux = norm_flux[np.where(norm_flux > cutoff)]
        norm_fluxes.append(norm_flux)
        norm_times.append(norm_time)
    norm_fluxes = np.array(norm_fluxes)
    norm_times = np.array(norm_times)
    
    if print_fig==True:        
        for i in range(len(norm_fluxes)):
            plt.figure()
            plt.subplot(2,1,1)
            plt.plot(times[i], fluxes[i], 'b,')
            plt.plot(times[i], med_fluxes[i], 'r-')
            plt.xlabel('Time [days]')
            plt.ylabel('Flux')
            plt.axis(xmin=min(times[i]),
                     xmax=max(times[i]),
                     ymin=np.median(fluxes[i])-300,
                     ymax=np.median(fluxes[i])+300
                     )
            plt.subplot(2,1,2)
            plt.plot(norm_times[i],norm_fluxes[i],'k,')
            plt.xlabel('Time [days]')
            plt.ylabel('Normalized Flux')
            plt.axis(xmin=min(norm_times[i]),
                     xmax=max(norm_times[i]),
                     ymin=0.98,
                     ymax=1.02
                     )            
            plt.tight_layout()
            plt.show()
            if save_fig==True:
                plt.savefig('figures/' + timestamp + '_norm_curve' + str(i) + '.pdf')

    return norm_times, norm_fluxes
    
#%% Define script that can remove data
    
def remove_bad_data(times, fluxes, bad_data):
    new_times = []
    new_fluxes = []
    for i in range(len(times)):
        bad_index = np.array([])
        for j in range(len(bad_data[:,0])):
            bad_index = np.append(bad_index, np.where(np.logical_and(times[i]>bad_data[j,0], times[i]<bad_data[j,1])))
        print(bad_index.shape)        
        new_time = np.delete(times[i], bad_index)
        new_times.append(new_time)
        new_flux = np.delete(fluxes[i], bad_index)
        new_fluxes.append(new_flux)
    new_times = np.array(new_times)
    new_fluxes = np.array(new_fluxes)
    
    return new_times, new_fluxes
    
#%% Define script that can interpolate data to even time step

def interpolate_tess(norm_times, norm_fluxes, print_fig=False, save_fig=False):

    even_times = []
    even_fluxes = []

    for i in range(len(norm_times)):
        even_time = np.linspace(norm_times[i][0], norm_times[i][-1], 22000)
        flux_func = interp1d(norm_times[i], norm_fluxes[i], kind='linear')
        even_flux = flux_func(even_time)
        even_times.append(even_time)
        even_fluxes.append(even_flux)

    if print_fig==True:
        for i in range(len(norm_times)):        
            plt.figure()
            plt.plot(even_times[i], even_fluxes[i], 'r,')
            plt.plot(norm_times[i], norm_fluxes[i], 'b,')
            plt.show()
            if save_fig==True:
                plt.savefig('figures/' + timestamp + '_inter_curve' + str(i) + '.pdf')
                
    return even_times, even_fluxes