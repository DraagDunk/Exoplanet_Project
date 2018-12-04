# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 11:31:16 2018

@author: jesper
"""

import numpy as np
import matplotlib.pyplot as plt
import exo_functions as ex

from PyAstronomy.modelSuite import forTrans as ft
from PyAstronomy import funcFit as fuf

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

times, fluxes = ex.fine_mesh_filter_tess(times, fluxes, n_sigmas, print_fig=False, save_fig=False)
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

centroids = ex.find_peaks(correlation_x, correlation_y, thresholds, print_fig=False, save_fig=False)

periods = centroids*time_steps
print('Periods:' + str(periods))

binned_times, binned_fluxes = ex.bin_fluxes_and_times_tess(norm_times, norm_fluxes, periods, [0, 0, 1, 0], 
                                                           print_fig=False, save_fig=False)


#%%
ma = ft.MandelAgolLC(orbit="keplerian", ld="quad")
data=np.genfromtxt("rap6.ascii", unpack=True)
time = binned_times[0]
flux = binned_fluxes[0]

period = periods[0]
a_in_stellar_rad = 10
inclination = 60
radius_ratio = np.sqrt(1 - np.min(flux))
lin_Limb = 0.4048
quad_Limb = 0.2695
T0 = np.mean(time)
eccentricity = 0.10


ma["per"] = period; ma["i"] = inclination; ma["a"] = a_in_stellar_rad;
ma["p"] = radius_ratio; ma["linLimb"] = lin_Limb; ma["quadLimb"] = quad_Limb;
ma["e"] = eccentricity; ma["T0"] = T0

data_til_err = flux[0:60]
err = np.std(data_til_err)


errors = err*np.ones(len(flux))
ma.thaw(["i","a","p","linLimb","quadLimb","e","T0"])

#==============================================================================
# 
# X0 = {"i_Occultquad":ma["i_Occultquad"],"a_Occultquad":ma["a_Occultquad"],\
#     "p_Occultquad":ma["p_Occultquad"]}
# 
# Lims = {"i_Occultquad":[80.,90.],"a_Occultquad":[0.,15],"p_Occultquad":[0.,0.5]}
#     
# steps = {"i_Occultquad":0.5,"a_Occultquad":0.5,"p_Occultquad":0.01}
# 
# ma.fitMCMC(time,flux,X0,Lims,steps,errors,iter = 100000,dbfile = "mcmctest3.tmp"\
#  ,burn = 0, quiet = True)
# 
# ma.parameterSummary()
# #%%
# MCMC_fit = ma.model
# 
# data_til_err = flux[0:60]
# err = np.std(data_til_err) 
# 
# n = len(flux)
# chi_square = sum(((flux-MCMC_fit)**2)/err**2)
# k = 4
# 
# plt.plot(time,flux,'.b')
# plt.plot(time,MCMC_fit,'-r')
# plt.show()
# 
# BIC = chi_square + k*np.log(n) 
# print(BIC) 
# #%%
# db = pymc.database.pickle.load('mcmctest3.tmp')
# ta = fuf.TraceAnalysis("mcmctest3.tmp")
# 
# plt.hist(db.trace("i_Occultquad", 0)[:])
# plt.show()
# ta.plotTrace("i_Occultquad")
# ta.show()
# 
# plt.hist(db.trace("a_Occultquad", 0)[:])
# plt.show()
# ta.plotTrace("a_Occultquad")
# ta.show()
# 
# plt.hist(db.trace("p_Occultquad", 0)[:])
# plt.show()
# ta.plotTrace("p_Occultquad")
# ta.show()
# 
# #%%
# 
# ta.correlationTable(coeff = "spearman")
# ta.correlationTable(coeff = "pearson")
# 
# #%%
# print("Available parameters: ", ta.availableParameters())
# 
# for p in ta.availableParameters():
#   hpd = ta.hpd(p, cred=0.95)
#   print("Parameter %5s, mean = % g, median = % g, std = % g, 95%% HPD = % g - % g" \
#         % (p, ta.mean(p), ta.median(p), ta.std(p), hpd[0], hpd[1]))
# 
#==============================================================================
