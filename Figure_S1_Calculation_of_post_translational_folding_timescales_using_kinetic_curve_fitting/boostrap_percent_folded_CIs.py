# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 13:04:13 2021

@author: nissley
"""

import sys, os
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

fdirs = ['fpt_analysis/fpt_analysis_1yta_wt/'  ,
         'fpt_analysis/fpt_analysis_1yta_fast/',
         'fpt_analysis/fpt_analysis_1yta_slow/']

iterations = 1000000 # 10^6

for fdir in fdirs:
    
    start = datetime.now()
    
    # determine number of unfolded trajectories at the final time point of the survival probability curve
    N_unfolded = int(np.loadtxt(fdir+'frac_U_vs_time_200trajs.txt')[:,2][-1]/0.005)
    
    print (fdir, N_unfolded)
    
    # populate 200-element 0D-array with a '1' for each unfolded trajectory
    data = np.zeros((200))
    for i in range (0, N_unfolded):
        data[i] = 1.0
        
    # compute true percent unfolded
    val = np.mean(data)

    val_list = []

    # bootstrap - selection WITH REPLACEMENT from the original data set
    for i in range (0, iterations):
    
        # make a random selection THE SAME SIZE AS THE ORIGINAL DATA SET with replacement
        new_data = np.random.choice(data, size=len(data), replace=True, p=None)
        new_val = np.mean(new_data)
        val_list.append(new_val)
    
    lower = np.percentile(val_list, 2.5)
    upper = np.percentile(val_list, 97.5)

    print (fdir+':', '%.5f' %val, '%.5f' %lower, '%.5f' %upper)

    print ('Done bootstrapping:', datetime.now()-start)