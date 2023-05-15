# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 19:49:57 2021

@author: nissley
"""

import os, sys
import matplotlib.pyplot as plt
import numpy as np

from matplotlib.pyplot import figure
figure(num=None, figsize=(2, 1.75), dpi=300, facecolor='w', edgecolor='k')

fontsize = 5
lw = 0.7
alpha = 0.7
tick_width = 0.75
spine_width = 0.75

# function to take codon sequence and return translation schedule
def getTransTimeSchedule(codon_seq):
    
    # codon_seq: str, an input
    
    mapCodonToTime, avg_time = getCodonTimes('fluitt_trans_times_mean_840000.txt')
    
    # assume sequence is in one line (and the first line at that)
    codon_seq = readFile(codon_seq)[0].strip()
    
    trans_schedule = []
    
    ## NORMALIZE CODON TIMES SO THE AVERAGE IS 1.0???
    
    for i in range (0, len(codon_seq), 3):
        codon = codon_seq[i:i+3]
        #print (codon, mapCodonToTime[codon])
        trans_schedule.append(mapCodonToTime[codon]/avg_time)
        
    return trans_schedule

# read translation times from file
def getCodonTimes(inpFile):
    
    # inpFile: str, path to file with codon trans times; assumed to be formatted as " codon      time "
    
    # returns: dictionary; keys are mRNA codons, values are translation times
    
    temp_map = {}
    
    dfile = readFile(inpFile)
    total = 0.0
    for line in dfile:
        line = line.split()
        codon = line[0]
        ta = float(line[1])
        temp_map[codon] = ta
        total += ta
    avg = total/64.
        
    return temp_map, avg

# function to read a file into memory as a list of lines
def readFile(filePath):
    
    # filePath: str, path to file to be read in
    
    # returns: contents of filePath as a list of lines
    
    if os.path.exists(filePath):
        ftemp = open(filePath)
        temp = ftemp.readlines()
        ftemp.close()
        return temp
    else:
        print (filePath, 'does not exist.')
        sys.exit()

# function to compute moving average
def ma(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

# set of three colors with best contrast

#          blue      yellow     red
colors = ['#004488', '#DDAA33', '#BB5566']

# width for running average
w = 15

# 1YTA

wt_1yta   = getTransTimeSchedule(  '1yta_wt_mrna_sequence.txt')
fast_1yta = getTransTimeSchedule('1yta_fast_mrna_sequence.txt')
slow_1yta = getTransTimeSchedule('1yta_slow_mrna_sequence.txt')
xvals = np.arange(1.0, len(wt_1yta)+1.0, 1.0)

plt.clf()
ax = plt.subplot()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(spine_width)
ax.spines['bottom'].set_linewidth(spine_width)
ax.xaxis.set_tick_params(width=tick_width)
ax.yaxis.set_tick_params(width=tick_width)
plt.plot(ma(xvals, w), ma(fast_1yta, w), color=colors[2], label='FAST', linewidth=lw, alpha=1.0)
plt.plot(ma(xvals, w), ma(  wt_1yta, w), color=colors[1], label='WT', linewidth=lw, alpha=1.0)
plt.plot(ma(xvals, w), ma(slow_1yta, w), color=colors[0], label='SLOW', linewidth=lw, alpha=1.0)
plt.xlabel('Codon position', fontsize=fontsize)
plt.ylabel('Normalized mean\ntranslation time', fontsize=fontsize)
plt.ylim(0.0, 2.5)
plt.xticks([0, 60, 120, 181], fontsize=fontsize)
plt.yticks([0.0, 0.5, 1.0, 1.5, 2.0, 2.5], fontsize=fontsize)
plt.legend(framealpha=0.0, loc='upper left', fontsize=4.5, ncol=3)
plt.tight_layout()
plt.savefig('1yta_trans_schedule.png', dpi=300)

# 2IS3

wt_2is3   = getTransTimeSchedule('2is3_wt_mrna_sequence.txt')
fast_2is3 = getTransTimeSchedule('2is3_fast_mrna_sequence.txt')
slow_2is3 = getTransTimeSchedule('2is3_slow_mrna_sequence.txt')
xvals = np.arange(1.0, len(wt_2is3)+1.0, 1.0)

plt.clf()
ax = plt.subplot()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(spine_width)
ax.spines['bottom'].set_linewidth(spine_width)
ax.xaxis.set_tick_params(width=tick_width)
ax.yaxis.set_tick_params(width=tick_width)
plt.plot(ma(xvals, w), ma(  wt_2is3, w), color=colors[1], label='WT', linewidth=lw, alpha=1.0)
plt.plot(ma(xvals, w), ma(fast_2is3, w), color=colors[2], label='FAST', linewidth=lw, alpha=1.0)
plt.plot(ma(xvals, w), ma(slow_2is3, w), color=colors[0], label='SLOW', linewidth=lw, alpha=1.0)
plt.xlabel('Codon position', fontsize=fontsize)
plt.ylabel('Normalized mean\ntranslation time', fontsize=fontsize)
plt.ylim(0.0, 2.5)
plt.xticks([0, 55, 110, 165, 215], fontsize=fontsize)
plt.yticks([0.0, 0.5, 1.0, 1.5, 2.0, 2.5], fontsize=fontsize)
plt.tight_layout()
plt.savefig('2is3_trans_schedule.png', dpi=300)