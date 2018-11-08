# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 11:31:16 2018

@author: jesper
"""

import numpy as np
import matplotlib.pyplot as plt
import exo_functions as ex

plt.close('all')

paths = np.array(['tess2018206045859-s0001-0000000038846515-111-s_llc.fits',
                  'tess2018206045859-s0001-0000000089020549-111-s_llc.fits',
                  'tess2018234235059-s0002-0000000201248411-0121-s_lc.fits'])
                 
for i in range(len(paths)):
    ex.import_tess_fits(paths[i])