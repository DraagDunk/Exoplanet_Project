# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 11:31:16 2018

@author: jesper
"""

import numpy as np
import matplotlib.pyplot as plt
import exo_functions as ex

plt.close('all')

paths = np.genfromtxt('data/fits.txt', dtype=str)

times =     []
fluxes =    [] 
                 
for i in range(len(paths)):
    time, flux = ex.import_tess_fits(paths[i], print_fig=False)
    times.append(time)
    fluxes.append(flux)

times = np.array(times)
fluxes = np.array(fluxes)

intervals = [500,250,150]

norm_times, norm_fluxes = ex.normer_fluxes(times,fluxes,intervals,print_fig=True, save_fig=False)