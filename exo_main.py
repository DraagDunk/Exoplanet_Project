# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 11:31:16 2018

@author: jesper
"""

import numpy as np
import matplotlib.pyplot as plt
import exo_functions as ex
#==============================================================================
# import pymc
# from PyAstronomy.modelSuite import forTrans as ft
# from PyAstronomy import funcFit as fuf
#==============================================================================

plt.close('all')

#%% TODO

# Normer korrelationsspektrum

# EKSTRA
# Find mere nøjagtig periode (Lineær fit til korrelationspeaks)

paths = np.genfromtxt('data/fits.txt', dtype=str)

TICs = []

for i in range(len(paths)):
    TICs.append(paths[i][24:40])
    while TICs[-1][0] == '0':
        TICs[-1] = TICs[-1][1:]

times =     []
fluxes =    [] 
                 
for i in range(len(paths)):
    time, flux = ex.import_tess_fits(paths[i], TICs[i], print_fig=True, save_fig=False)
    times.append(time)
    fluxes.append(flux)

times = np.array(times)
fluxes = np.array(fluxes)

intervals = [1.42,0.35,0.55,0.2]

n_sigmas = [0.35, 0.3, 0.25, 0.5]

# width SKAL være ulige!!!
times, fluxes = ex.fine_mesh_filter_tess(times, fluxes, n_sigmas, TICs, width=11, print_fig=True, save_fig=False)

norm_times, norm_fluxes = ex.normer_fluxes(times,fluxes,intervals,TICs,'lightcurve',cutoff = 0.985,print_fig=True, save_fig=False)
#%%
bad_data = np.array([[1338.5, 1339.7],
                     [1347.1, 1349.4],
                     [1367.1, 1368.65]])
                     
norm_fluxes = ex.remove_bad_data(norm_times, norm_fluxes, bad_data)

even_times, even_fluxes = ex.interpolate_tess(norm_times, norm_fluxes, TICs, print_fig=False, save_fig=False)

time_steps = even_times[:,1] - even_times[:,0]

correlation_x, correlation_y = ex.correlate_tess(even_fluxes, time_steps, TICs, print_fig=False, save_fig=False)

corr_intervals = 500*np.ones(4)
correlations_x = []
for i in range(4):
    correlations_x.append(correlation_x)
    
correlations_x = np.array(correlations_x)

norm_corr_x, norm_corr_y = ex.normer_fluxes(correlations_x, correlation_y, corr_intervals,TICs,'correlation',cutoff = -1,print_fig = False, save_fig = False)

thresholds = np.array([3, 1.1, 2.5, 2])

alt_peaks = [np.array([]), 
             np.array([(np.abs(norm_corr_x[1]-1673)).argmin(), (np.abs(norm_corr_x[1]-3341)).argmin(), (np.abs(norm_corr_x[1]-5008)).argmin(), (np.abs(norm_corr_x[1]-6676)).argmin()]), 
             np.array([]), 
             np.array([])
             ]

periods = ex.find_peaks(correlation_x, norm_corr_y, thresholds, alt_peaks, TICs, print_fig=False, save_fig=False)

periods = periods.T*time_steps
periods = periods.reshape(len(periods[0]))

print('Periods:' + str(periods))

binned_times, binned_fluxes = ex.bin_fluxes_and_times_tess(norm_times, norm_fluxes, periods, [0, 0, 1, 0], TICs, 
                                                           print_fig=False, save_fig=False)


#==============================================================================
# 
# #%%
# ma = ft.MandelAgolLC(orbit="keplerian", ld="quad")
# time = binned_times[0]
# flux = binned_fluxes[0]
# 
# period = periods[0]
# a_in_stellar_rad = 10
# inclination = 80
# radius_ratio = np.sqrt(1 - np.min(flux))
# lin_Limb = 0.4048
# quad_Limb = 0.2695
# T0 = np.mean(time)
# eccentricity = 0.10
# 
# 
# ma["per"] = period; ma["i"] = inclination; ma["a"] = a_in_stellar_rad;
# ma["p"] = radius_ratio; ma["linLimb"] = lin_Limb; ma["quadLimb"] = quad_Limb;
# ma["e"] = eccentricity; ma["tau"] = T0
# 
# data_til_err = flux[0:60]
# err = np.std(data_til_err)
# t_offset = time_steps[0]*len(time)/2
# 
# errors = err*np.ones(len(flux))
# ma.thaw(["i","a","p","e","tau"])
# 
# 
# X0 = {"i":ma["i"],"a":ma["a"],"p":ma["p"],"e":ma["e"],"tau":ma["tau"]}
# 
# Lims = {"i":[75.,90.],"a":[5.,100],"p":[0.,0.5],"e":[0.,0.7],"tau":[T0-t_offset,T0+t_offset]}
#     
# steps = {"i":0.5,"a":0.5,"p":0.01,"e":0.01,"tau":time_steps[0]}
# 
# ma.fitMCMC(time,flux,X0,Lims,steps,errors,iter = 100000,dbfile = "mcmc_test.tmp"\
#  ,burn = 0, quiet = True)
# 
# ma.parameterSummary()
# 
# MCMC_fit = ma.model
# 
# 
# #%%
# 
# plt.figure()
# plt.plot(time,flux,'.b')
# plt.plot(time,MCMC_fit,'-r')
# plt.title('MCMC fit of MA-model to transit')
# plt.xlabel('Time % period [days]')
# plt.ylabel('Normalized flux')
# plt.show()
# 
# db = pymc.database.pickle.load('mcmc_test.tmp')
# ta = fuf.TraceAnalysis("mcmc_test.tmp")
# 
# print("Available parameters: ", ta.availableParameters())
# 
# for p in ta.availableParameters():
#   hpd = ta.hpd(p, cred=0.95)
#   print("Parameter %5s, mean = % g, median = % g, std = % g, 95%% HPD = % g - % g" \
#         % (p, ta.mean(p), ta.median(p), ta.std(p), hpd[0], hpd[1]))
# 
# #%%
# plt.figure()
# plt.hist(db.trace("i", 0)[:])
# plt.show()
# 
# #%%
# plt.figure()
# ta.plotTrace("i")
# ta.show()
# 
# plt.figure()
# plt.hist(db.trace("a", 0)[:])
# plt.show()
# 
# plt.figure()
# ta.plotTrace("a")
# ta.show()
# 
# plt.figure()
# plt.hist(db.trace("p", 0)[:])
# plt.show()
# 
# plt.figure()
# ta.plotTrace("p")
# ta.show()
# 
# #%%
# 
# ta.correlationTable(coeff = "spearman")
# ta.correlationTable(coeff = "pearson")
# 
# #%%
# 
# 
# 
#==============================================================================
