# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 11:31:16 2018

@author: jesper
"""

import numpy as np
import matplotlib.pyplot as plt
import exo_functions as ex

plt.close('all')

#%% TODO

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

norm_times, norm_fluxes = ex.normer_fluxes(times,fluxes,intervals,cutoff = 0.985,print_fig=False, save_fig=False)

bad_data = np.array([[1338.5, 1339.7],
                     [1347.1, 1349.4],
                     [1367.1, 1368.65]])
                     
norm_times, norm_fluxes = ex.remove_bad_data(norm_times, norm_fluxes, bad_data)

even_times, even_fluxes = ex.interpolate_tess(norm_times, norm_fluxes, print_fig=False, save_fig=False)

time_steps = even_times[:,1] - even_times[:,0]

correlation_x, correlation_y = ex.correlate_tess(even_fluxes, time_steps, print_fig=True, save_fig=False)

thresholds = np.array([0.003, 0.003, 0.003])

peaks = ex.find_peaks(correlation_x, correlation_y, thresholds, print_fig=True)