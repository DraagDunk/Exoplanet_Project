# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 11:31:16 2018

@author: jesper
"""

import numpy as np
import matplotlib.pyplot as plt
import exo_functions as ex

plt.close('all')

paths = np.genfromtxt('fits.txt', dtype=str)
                 
for i in range(len(paths)):
    ex.import_tess_fits(paths[i])