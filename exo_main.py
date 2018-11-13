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

intervals = [1.42,0.35,0.2]

norm_times, norm_fluxes = ex.normer_fluxes(times,fluxes,intervals,cutoff = 0.99,print_fig=True, save_fig=False)

#%%
even_times = np.linspace(norm_times[0][0], norm_times[0][-1], 22000)
flux_func = interp1d(norm_times[0], norm_fluxes[0], kind='linear')
even_fluxes = flux_func(even_times)

plt.figure()
plt.plot(even_times, even_fluxes, 'r,')
plt.plot(norm_times[0], norm_fluxes[0], 'b,')
plt.show()