# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 19:49:59 2021

@author: nissley
"""

import os, sys
import matplotlib.pyplot as plt
import numpy as np

from matplotlib.pyplot import figure
figure(num=None, figsize=(2,1.75), dpi=300, facecolor='w', edgecolor='k')

wt_1yta = np.loadtxt('1yta_wt_ave_Q.dat')[:,1]
fast_1yta = np.loadtxt('1yta_fast_ave_Q.dat')[:,1]
slow_1yta = np.loadtxt('1yta_slow_ave_Q.dat')[:,1]

wt_2is3 = np.loadtxt('2is3_wt_ave_Q.dat')[:,1]
fast_2is3 = np.loadtxt('2is3_fast_ave_Q.dat')[:,1]
slow_2is3 = np.loadtxt('2is3_slow_ave_Q.dat')[:,1]

times = np.loadtxt('1yta_wt_ave_Q.dat')[:,0]

#          blue      yellow     red
colors = ['#004488', '#DDAA33', '#BB5566']

fontsize = 5
lw = 0.7
alpha = 0.7
tick_width = 0.75
spine_width = 0.75
alpha = 1
width = 0.5
capsize = 2.0
ewidth = 1.
msize = 4
plt.clf()

ax = plt.subplot()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(spine_width)
ax.spines['bottom'].set_linewidth(spine_width)
ax.xaxis.set_tick_params(width=tick_width)
ax.yaxis.set_tick_params(width=tick_width)

plt.plot(times, wt_1yta, color=colors[1], linewidth=width, alpha=0.85, label='WT')
plt.plot(times, fast_1yta, color=colors[2], linewidth=width, alpha=0.85, label='FAST')
plt.plot(times, slow_1yta, color=colors[0], linewidth=width, alpha=0.85, label='SLOW')

#plt.ylim(0.0, 1.0)
plt.yticks([0.75, 0.80, 0.85, 0.9])
plt.ylabel('Fraction of native contacts, Q', fontsize=fontsize)

plt.xlim(-0.1, 5.05)
plt.xticks([0, 1, 2, 3, 4, 5], fontsize=fontsize)
plt.xlabel('Time since release, '+r'$\mu$s', fontsize=fontsize)
plt.yticks(fontsize=fontsize)

#plt.legend(loc='lower right', fancybox=True, framealpha=0.5, fontsize=fontsize)
plt.tight_layout()
plt.savefig('1yta_Q.png', dpi=300)

plt.clf()

ax = plt.subplot()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(spine_width)
ax.spines['bottom'].set_linewidth(spine_width)
ax.xaxis.set_tick_params(width=tick_width)
ax.yaxis.set_tick_params(width=tick_width)

plt.plot(times, wt_2is3, color=colors[1], linewidth=width, alpha=0.85, label='WT')
plt.plot(times, fast_2is3, color=colors[2], linewidth=width, alpha=0.85, label='FAST')
plt.plot(times, slow_2is3, color=colors[0], linewidth=width, alpha=0.85, label='SLOW')

#plt.ylim(0.0, 1.0)
plt.yticks([0.50, 0.60, 0.70, 0.80, 0.90], ['0.50', '0.60', '0.70', '0.80', '0.90'])
plt.ylabel('Fraction of native contacts, Q', fontsize=fontsize)

plt.xlim(-0.1, 5.05)
plt.xticks([0, 1, 2, 3, 4, 5], fontsize=fontsize)
plt.xlabel('Time since release, '+r'$\mu$s', fontsize=fontsize)
plt.yticks(fontsize=fontsize)

plt.tight_layout()
plt.savefig('2is3_Q.png', dpi=300)