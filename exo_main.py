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

# Normer korrelationsspektrum

# EKSTRA
# Find mere nøjagtig periode (Lineær fit til korrelationspeaks)

paths = np.genfromtxt('data/fits.txt', dtype=str)

times =     []
fluxes =    [] 
                 
for i in range(len(paths)):
    time, flux = ex.import_tess_fits(paths[i], print_fig=False)
    times.append(time)
    fluxes.append(flux)

times = np.array(times)
fluxes = np.array(fluxes)

intervals = [1.42,0.35,0.7,0.2]

n_sigmas = [0.35, 0.3, 0.25, 0.5]

times, fluxes = ex.fine_mesh_filter_tess(times, fluxes, n_sigmas, print_fig=True, save_fig=True)
#%%
norm_times, norm_fluxes = ex.normer_fluxes(times,fluxes,intervals,cutoff = 0.985,print_fig=False, save_fig=False)

bad_data = np.array([[1338.5, 1339.7],
                     [1347.1, 1349.4],
                     [1367.1, 1368.65]])
                     
norm_fluxes = ex.remove_bad_data(norm_times, norm_fluxes, bad_data)

even_times, even_fluxes = ex.interpolate_tess(norm_times, norm_fluxes, print_fig=False, save_fig=False)

time_steps = even_times[:,1] - even_times[:,0]

correlation_x, correlation_y = ex.correlate_tess(even_fluxes, time_steps, print_fig=False, save_fig=False)

thresholds = np.array([0.003, 0.003, 0.003, 0.003])

centroids = ex.find_peaks(correlation_x, correlation_y, thresholds, print_fig=True, save_fig=True)

periods = centroids*time_steps
print('Periods:' + str(periods))

binned_times, binned_fluxes = ex.bin_fluxes_and_times_tess(norm_times, norm_fluxes, periods, [0, 0, 1, 0], 
                                                           print_fig=True, save_fig=True)