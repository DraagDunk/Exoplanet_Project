# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 10:53:42 2018

@author: jesper
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

#%% Import function

def import_tess_fits(path):
    hdulist = fits.open(path)
    hdu = hdulist[1]
    
    hdu_data = hdu.data
    
    time = []
    sap_flux = []
    pdc_flux = []
    
    for i in range(len(hdu_data)):
        time.append(hdu_data[i][0])
        sap_flux.append(hdu_data[i][3])
        pdc_flux.append(hdu_data[i][7])
        
    time =      np.array(time)
    sap_flux =  np.array(sap_flux)
    pdc_flux =  np.array(pdc_flux)
        
    plt.figure()
    plt.plot(time, sap_flux)
    plt.title(path)
    plt.xlabel('time [days]')
    plt.ylabel('flux')
    
    return time, sap_flux