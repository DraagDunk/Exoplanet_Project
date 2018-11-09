# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 10:53:42 2018

@author: jesper
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

# Get current time
import datetime
now = datetime.datetime.now()

#%% Import function

def import_tess_fits(path,print_fig=False):
    hdulist = fits.open(path)
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
def normer_fluxes(times,fluxes,intervals,print_fig=False,save_fig=False):
    med_fluxes = []
    for i in range(len(fluxes)):
        med_flux = []
        for j in range(len(fluxes[i])-2*intervals[i]):
            med_flux.append(np.median(fluxes[i][j:j+2*intervals[i]]))
        med_flux = np.array(med_flux)    
        med_fluxes.append(med_flux)
    
    med_fluxes=np.array(med_fluxes)
    
    norm_fluxes = []
    norm_times = []
    for i in range(len(fluxes)):
        norm_flux = fluxes[i][intervals[i]:-intervals[i]]/med_fluxes[i]
        norm_fluxes.append(norm_flux)
        norm_time = times[i][intervals[i]:-intervals[i]]
        norm_times.append(norm_time)
    norm_fluxes = np.array(norm_fluxes)
    norm_times = np.array(norm_times)
        
    if print_fig==True:        
        for i in range(len(norm_fluxes)):
            plt.figure()
            plt.subplot(2,1,1)
            plt.plot(norm_times[i], fluxes[i][intervals[i]:-intervals[i]], 'b,')
            plt.plot(norm_times[i], med_fluxes[i], 'r-')
            plt.xlabel('Time [days]')
            plt.ylabel('Flux')
            plt.axis(xmin=min(norm_times[i]),
                     xmax=max(norm_times[i]),
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
                plt.savefig('figures/' + 
                            str(now.year) + 
                            '-' + 
                            str(now.month) + 
                            '-' + 
                            str(now.day) +
                            '_' +
                            str(now.hour) +
                            ':' +
                            str(now.minute) +
                            ':' +
                            str(now.second) +
                            '_norm_curve' +
                            str(i) +
                            '.pdf')

    return norm_times, norm_fluxes