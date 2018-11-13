# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 11:31:16 2018

@author: jesper
"""

import numpy as np
import matplotlib.pyplot as plt
import exo_functions as ex
from scipy.interpolate import interp1d

plt.close('all')

#%% TODO

# Reparer medianfilter
# Autokorreler
# Find periode
# Fold til phaseplot

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

#even_times = np.linspace(min(norm_times[0]), max(norm_times[0]), 40000)
#flux_func = interp1d(norm_times[0], norm_fluxes[0], kind='square')
#even_fluxes = flux_func(even_times)

#plt.figure()
#plt.plot(norm_times, norm_fluxes, 'b,')
#plt.plot(even_times, even_fluxes, 'r-')
#plt.show()