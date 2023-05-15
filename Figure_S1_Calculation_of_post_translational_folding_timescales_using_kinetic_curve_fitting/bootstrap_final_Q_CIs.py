# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 12:20:37 2021

@author: nissley
"""

import sys, os
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

fdirs = ['fpt_analysis/fpt_analysis_1yta_wt/props_ts/'  ,
         'fpt_analysis/fpt_analysis_1yta_fast/props_ts/',
         'fpt_analysis/fpt_analysis_1yta_slow/props_ts/',
         'fpt_analysis/fpt_analysis_2is3_wt/props_ts/'  ,
         'fpt_analysis/fpt_analysis_2is3_fast/props_ts/',
         'fpt_analysis/fpt_analysis_2is3_slow/props_ts/']

iterations = 1000000 # 10^6

for fdir in fdirs:

    start = datetime.now()    

    data = np.zeros((200))

    for j in range (0, 200):
        data[j] = np.loadtxt(fdir+'props_t'+str(j+1)+'.dat')[-1]
        #print (j+1, data[j])    

    print ('Done collecting values', datetime.now()-start)

    mean = np.mean(data)

    print ('Mean:', mean)
    means_list = []

    # bootstrap - selection WITH REPLACEMENT from the original data set
    for i in range (0, iterations):
    
        # make a random selection THE SAME SIZE AS THE ORIGINAL DATA SET
        new_data = np.random.choice(data, size=len(data), replace=True, p=None)
        new_mean = np.mean(new_data)
        means_list.append(new_mean)

    lower = np.percentile(means_list, 2.5)
    upper = np.percentile(means_list, 97.5)

    print (fdir+':', '%.5f' %mean, '%.5f' %lower, '%.5f' %upper)

    print ('Done bootsrapping', datetime.now()-start)
    print ()