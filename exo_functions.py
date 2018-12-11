# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 10:53:42 2018

@author: jesper
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit

# Get current time
import datetime
now = datetime.datetime.now()
timestamp = str(now.year) + '-' + str(now.month) + '-' + str(now.day) + '_' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second)

#%% Function that imports data from .fits file, and returns time and SAP flux

# path: name of .fits file
def import_tess_fits(path,print_fig=False):
    hdulist = fits.open('data/'+ path)
    hdu = hdulist[1]
    
    hdu_data = hdu.data
    
    time = []
    sap_flux = []
    pdc_flux = []
    
    # Remove NaN from data
    for i in range(len(hdu_data)):
        if np.isnan(hdu_data[i][3])==False:
            time.append(hdu_data[i][0])
            sap_flux.append(hdu_data[i][3])
            pdc_flux.append(hdu_data[i][7])
        
    time =      np.array(time)
    sap_flux =  np.array(sap_flux)
    pdc_flux =  np.array(pdc_flux)
    
    # Plot raw data
    if print_fig==True:    
        plt.figure()
        plt.plot(time, sap_flux, ',')
        plt.title(path)
        plt.xlabel('time [days]')
        plt.ylabel('flux')
    
    return time, sap_flux
    
#%% Function that removes noise through a fine-mesh median filter

# times:    Array of arrays of time 
# fluxes:   Array of arrays of flux
# n_sigmas: List of number of standard deviations that should be included in data
def fine_mesh_filter_tess(times, fluxes, n_sigmas, TICs, print_fig=False, save_fig=False):
    med_fluxes = []
    for i in range(len(fluxes)):
        med_flux = []
        for j in range(len(fluxes[i])-2):
            med_flux.append(np.median(fluxes[i][j:j+3]))
        med_flux = np.array(med_flux)
        
        med_fluxes.append(med_flux)
    
    med_fluxes = np.array(med_fluxes)
    
    # Keep old data before filtering
    old_fluxes = np.copy(fluxes)
    old_times = np.copy(times)
    
    # Remove data further from median than n_sigmas specifies
    for i in range(len(med_fluxes)):
        fluxes[i] = fluxes[i][1:-1]
        times[i] = times[i][1:-1]  
        norm_flux = fluxes[i]/med_fluxes[i]
        sigma = np.std(norm_flux)
        fluxes[i] = fluxes[i][np.where(np.logical_and(norm_flux-1 < n_sigmas[i]*sigma, 
                                                      norm_flux-1 > -n_sigmas[i]*sigma))]
        times[i] = times[i][np.where(np.logical_and(norm_flux-1 < n_sigmas[i]*sigma,
                                                    norm_flux-1 > -n_sigmas[i]*sigma))]                             
        
        # Plot 3 subplots of data filtering
        if print_fig == True:
            fig = plt.figure()
            plt.title('Fine-mesh normalization')
            # Old data and median filter
            ax1 = plt.subplot(3,1,1)
            ax1.set_xticklabels([])
            plt.plot(old_times[i], old_fluxes[i], 'k,')
            plt.plot(old_times[i][1:-1], med_fluxes[i], 'r-')
            plt.ylabel('Flux [ADU]')
            plt.xlim(old_times[i][0], old_times[i][-1])
            plt.ylim(np.median(old_fluxes[i])-300, np.median(old_fluxes[i])+300)
            # Normalized data and data cutoff boundaries
            ax2 = plt.subplot(3,1,2)
            ax2.set_xticklabels([])            
            plt.plot(old_times[i][1:-1], norm_flux, 'k,')
            plt.plot([times[i][0], times[i][-1]], np.array([n_sigmas[i]*sigma, n_sigmas[i]*sigma])+1, 'r--')            
            plt.plot([times[i][0], times[i][-1]], np.array([-n_sigmas[i]*sigma, -n_sigmas[i]*sigma])+1, 'r--')
            plt.ylabel('Normalized flux')
            plt.xlim(old_times[i][0], old_times[i][-1])
            plt.ylim(1-0.7*sigma, 1+0.7*sigma)
            # Comparison between old and new data
            plt.subplot(3,1,3)
            plt.plot(old_times[i], old_fluxes[i], 'k,')
            plt.plot(times[i], fluxes[i], 'r,')
            plt.xlabel('Time [days]')
            plt.ylabel('Flux [ADU]')
            plt.xlim(old_times[i][0], old_times[i][-1])
            plt.ylim(np.median(old_fluxes[i])-300, np.median(old_fluxes[i])+300)            
            fig.subplots_adjust(wspace=0, hspace=0)            
            plt.tight_layout()
            plt.show()
            if save_fig == True:
                plt.savefig('figures/' + timestamp + '_finemesh_TIC' + TICs[i] + '.pdf')
            
    
    return times, fluxes

#%% Script that can normalize light curve and return normalized flux

# times:        Array of arrays of time
# fluxes:       Array of arrays of flux
# intervals:    List of half-width of the time interval each point is calculated from in median filter
# cutoff:       Lower flux limit for inclusion of data
def normer_fluxes(times,fluxes,intervals,cutoff = 0.98,TICs,print_fig=False,save_fig=False):
    # Calculate median filter from all points within twice the "intervals"    
    med_fluxes = []
    for i in range(len(fluxes)):
        med_flux = []
        for j in range(len(fluxes[i])):
            med_flux.append(np.median(fluxes[i][np.where(np.logical_and(times[i] > times[i][j]-intervals[i],
            times[i] < times[i][j]+intervals[i]))]))
        med_flux = np.array(med_flux)    
        med_fluxes.append(med_flux)
    
    med_fluxes=np.array(med_fluxes)
    
    # Normalize data by dividing the data with the median filter
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
    
    # Plot 2 subplots
    if print_fig==True:        
        for i in range(len(norm_fluxes)):
            plt.figure()
            # Data and median filter
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
            # Normalized data
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
                plt.savefig('figures/' + timestamp + '_normcurve_TIC' + TICs[i] + '.pdf')

    return norm_times, norm_fluxes
    
#%% Script that can set normalized data in specific time periods to 1
    
# times:    Array of arrays of time
# fluxes:   Array of arrays of flux
# bad_data: 2D-array of times. Data between the two times in rows should be =1
def remove_bad_data(times, fluxes, bad_data):
    new_fluxes = []
    for i in range(len(times)):
        bad_index = np.array([])
        for j in range(len(bad_data[:,0])):
            bad_index = np.append(bad_index, np.where(np.logical_and(times[i]>bad_data[j,0], times[i]<bad_data[j,1])))       
        fluxes[i][bad_index.astype(int)] = 1
        new_fluxes.append(fluxes[i])
    new_fluxes = np.array(new_fluxes)
    
    return new_fluxes
    
#%% Script that can interpolate data to even time step and return interpolated data

# norm_times:   Array of arrays of time
# norm_fluxes:  Array of arrays of normalized flux
def interpolate_tess(norm_times, norm_fluxes, TICs, print_fig=False, save_fig=False):

    even_times = []
    even_fluxes = []
    
    # Interpolates data to 22000 points, equidistant along time axis
    for i in range(len(norm_times)):
        even_time = np.linspace(norm_times[i][0], norm_times[i][-1], 22000)
        flux_func = interp1d(norm_times[i], norm_fluxes[i], kind='linear')
        even_flux = flux_func(even_time)
        even_times.append(even_time)
        even_fluxes.append(even_flux)

    # Plot old data and interpolated data
    if print_fig==True:
        for i in range(len(norm_times)):        
            plt.figure()
            plt.plot(even_times[i], even_fluxes[i], 'r,')
            plt.plot(norm_times[i], norm_fluxes[i], 'b,')
            plt.show()
            if save_fig==True:
                plt.savefig('figures/' + timestamp + '_intercurve_TIC' + TICs[i] + '.pdf')
             
    even_times = np.array(even_times)
    even_fluxes = np.array(even_fluxes)

    return even_times, even_fluxes
    
#%% Script that can autocorrelate interpolated flux and return correlation spectra

# even_fluxes:  Array of arrays of interpolated flux
# time_steps:   List of distances between all points in time
def correlate_tess(even_fluxes, time_steps, TICs, print_fig=False, save_fig=False):

    # Autocorrelate all data below this normalized flux value
    zero_point = 1
    for i in range(len(even_fluxes)):
        even_fluxes[i][np.where(even_fluxes[i] > zero_point)] = zero_point
    
    correlation_spectra = []
    # Spectrum is inverted before correlation for positive signal
    for i in range(len(time_steps)):
        corr = np.correlate(zero_point-even_fluxes[i], zero_point-even_fluxes[i], mode='full')
        corr = corr[int((len(corr)-1)/2):-1]
        correlation_spectra.append(corr)
    correlation_spectra = np.array(correlation_spectra)
    
    correlation_x = np.linspace(0, (len(correlation_spectra[0])), len(correlation_spectra[0]))
    
    # Plot correlation spectrum
    if print_fig == True:
        for i in range(len(time_steps)):
            plt.figure()
            plt.plot(correlation_x*time_steps[i], correlation_spectra[i], 'b-')
            plt.xlabel('Period [days]')
            plt.ylabel('Correlation function')
            plt.tight_layout()
            plt.show()
            if save_fig == True:
                plt.savefig('figures/' + timestamp + '_correlation_TIC' + TICs[i] + '.pdf')
            
    return correlation_x, correlation_spectra
    
#%% Script that can find peaks in correlation spectrum and return a list of centroids

# Define gaussian function
def gaussian(x, a, b, c):
    return a * np.exp( - (x - b)**2 / (2*c**2) )

# correlation_x:    Array of arrays of steps taken in correlation
# correlation_y:    Array of arrays of correlation spectra
# thresholds:       List of vertical distance boundaries from peak to noise
def find_peaks(correlation_x, correlation_y, thresholds, TICs, print_fig=False, save_fig=False):
    peaks = []    
    # Find maxima that are farther from noise than the thresholds specify   
    for i in range(len(correlation_y[:,0])):
        peaks_temp = []        
        for j in range(len(correlation_y[i])-2):
            if correlation_y[i,j] <= correlation_y[i,j+1] >= correlation_y[i,j+2]:
                if correlation_y[i,j+1]-correlation_y[i,j-199] >= thresholds[i] and correlation_y[i,j+1]-correlation_y[i,j+201] >= thresholds[i]:
                    peaks_temp.append(j+1)
        peaks.append(np.array(peaks_temp))
    peaks = np.array(peaks)

    # Only allow one peak within 100 points; find the tallest one, delete all others
    for h in range(len(peaks)):
        if len(peaks[h]) > 0:
            delete_index = [1]
            while len(delete_index) > 0:
                delete_index = []
                for i in range(len(peaks[h])-1):
                    if abs(peaks[h][i]-peaks[h][i+1]) < 100:
                        if correlation_y[h][peaks[h][i]] < correlation_y[h][peaks[h][i+1]]:
                            delete_index.append(i)
                        elif correlation_y[h][peaks[h][i]] > correlation_y[h][peaks[h][i+1]]:
                            delete_index.append(i+1)
                delete_index = np.array(delete_index)
                peaks[h] = np.delete(peaks[h], delete_index)
                
    # Handy print of lists of peaks
    print('=========Array of found peaks:=========')
    print(peaks)
    print('=======================================')
            
    popts = []
    centroids = []
    # Fit gaussian functions to left-most peak (Redo to fit to all peaks and find period from all)
    for i in range(len(peaks)):
        if len(peaks[i]) > 0:
            popt, pcov = curve_fit(gaussian, correlation_x, correlation_y[i], 
                                   bounds = ([correlation_y[i][peaks[i][-1]]-0.01, peaks[i][-1]-20, 0],
                                             [correlation_y[i][peaks[i][-1]]+0.01, peaks[i][-1]+20, 50]))
            popts.append(popt)
            centroids.append(popt[1]/len(peaks[i]))
        else:
            popts.append(42)
            centroids.append(-1)
    
    # centroids aren't really centroids, but periods.            
    centroids = np.array(centroids)
    
    # Plot correlation spectra with gaussians and peaks
    if print_fig == True:
        for i in range(len(correlation_y[:,0])):
            if len(peaks[i]) > 0:
                plt.figure()
                plt.plot(correlation_x, correlation_y[i], 'b-')
                plt.plot(correlation_x[peaks[i]], correlation_y[i,peaks[i]], 'r*')
                plt.plot(correlation_x, gaussian(correlation_x, *popts[i]), 'k--')
                plt.xlabel('Points')
                plt.ylabel('Correlation function')
                plt.tight_layout()
                plt.show()
                if save_fig == True:
                    plt.savefig('figures/' + timestamp + '_peaks_TIC' + TICs[i] + '.pdf')
        
    return centroids
    
#%% Script that can bin and phasefold transit light curves and return binned data

# norm_times:   Array of arrays of time
# norm_fluxes:  Array of arrays of normalized flux
# periods:      List of periods
# offsets:      List of phase offsets
def bin_fluxes_and_times_tess(norm_times,norm_fluxes,periods,offsets,TICs,print_fig = False,save_fig = False):
    binned_times = []
    binned_fluxes = []
    for j in range(len(norm_times)):
        if periods[j] > 0: 
            phase_time = (norm_times[j]+offsets[j])%periods[j]
            # Bin length is set to length of the data set divided by the number of folds
            bins = np.linspace(0,periods[j],np.floor(len(norm_times[j])/((norm_times[j][-1]-norm_times[j][0])/periods[j])))
            
            binned_flux = []
            binned_time = []
            for i in range(len(bins)-1):
                binned_flux.append(np.median(norm_fluxes[j][np.where(np.logical_and(phase_time >= bins[i], phase_time < bins[i+1]))]))    
                binned_time.append(np.median(phase_time[np.where(np.logical_and(phase_time >= bins[i], phase_time < bins[i+1]))]))
    
            binned_flux = np.array(binned_flux)
            binned_time = np.array(binned_time)
        else:
            binned_time = -1
            binned_flux = -1
            
        binned_fluxes.append(binned_flux)
        binned_times.append(binned_time)
        
        # Plot normalized data and binned data
        if print_fig == True and periods[j] > 0:
            plt.figure()
            plt.plot(phase_time,norm_fluxes[j],',r')
            plt.plot(binned_time,binned_flux,'.k')
            plt.ylabel('Normalized median Flux')
            plt.xlabel('Time modulo Period [Days]')
            plt.title('Folded Light Curve')
            plt.xlim(0, periods[j])
            plt.tight_layout()
            plt.show()
            if save_fig == True:
                    plt.savefig('figures/' + timestamp + '_Folded_TIC' + TICs[j] + '.pdf')
    
    binned_fluxes = np.array(binned_fluxes)
    binned_times = np.array(binned_times)

    return binned_times, binned_fluxes
