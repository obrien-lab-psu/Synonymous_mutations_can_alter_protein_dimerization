# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 11:39:30 2021

@author: nissley
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
figure(num=None, figsize=(2, 1.75), dpi=300, facecolor='w', edgecolor='k')

fontsize = 7
lw = 0.7
alpha = 0.7
tick_width = 0.75
spine_width = 0.75
alpha = 1
width = 0.15
capsize = 2.0
ewidth = 1.
msize = 4

fast_1yta       = -60.501
fast_1yta_upper = -58.350
fast_1yta_lower = -62.471
fast_1yta_yerr  = [[np.abs(fast_1yta) - np.abs(fast_1yta_upper)], [np.abs(fast_1yta_lower) - np.abs(fast_1yta)]]

slow_1yta       = -64.127
slow_1yta_upper = -62.386
slow_1yta_lower = -65.687
slow_1yta_yerr  = [[np.abs(slow_1yta) - np.abs(slow_1yta_upper)], [np.abs(slow_1yta_lower) - np.abs(slow_1yta)]]

wt_1yta         = -58.261
wt_1yta_upper   = -55.940
wt_1yta_lower   = -60.315
wt_1yta_yerr    = [[np.abs(wt_1yta) - np.abs(wt_1yta_upper)], [np.abs(wt_1yta_lower) - np.abs(wt_1yta)]]

fast_2is3       = -34.826
fast_2is3_upper = -31.740
fast_2is3_lower = -37.682 
fast_2is3_yerr  = [[np.abs(fast_2is3) - np.abs(fast_2is3_upper)], [np.abs(fast_2is3_lower) - np.abs(fast_2is3)]]

slow_2is3       = -35.221
slow_2is3_upper = -32.324
slow_2is3_lower = -37.923
slow_2is3_yerr  = [[np.abs(slow_2is3) - np.abs(slow_2is3_upper)], [np.abs(slow_2is3_lower) - np.abs(slow_2is3)]]

wt_2is3         = -36.132
wt_2is3_upper   = -33.208
wt_2is3_lower   = -38.894
wt_2is3_yerr    = [[np.abs(wt_2is3) - np.abs(wt_2is3_upper)], [np.abs(wt_2is3_lower) - np.abs(wt_2is3)]]

#          blue      yellow     red
colors = ['#004488', '#DDAA33', '#BB5566']

plt.clf()

ax = plt.subplot()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(spine_width)
ax.spines['bottom'].set_linewidth(spine_width)
ax.xaxis.set_tick_params(width=tick_width)
ax.yaxis.set_tick_params(width=tick_width)

plt.errorbar(2, fast_1yta, marker='o', label='FAST', color=colors[2], alpha=alpha, yerr=fast_1yta_yerr, capsize=capsize, ecolor='k', elinewidth=ewidth, markersize=msize)
plt.errorbar(1,   wt_1yta, marker='o', label=  'WT', color=colors[1], alpha=alpha, yerr=  wt_1yta_yerr, capsize=capsize, ecolor='k', elinewidth=ewidth, markersize=msize)
plt.errorbar(3, slow_1yta, marker='o', label='SLOW', color=colors[0], alpha=alpha, yerr=slow_1yta_yerr, capsize=capsize, ecolor='k', elinewidth=ewidth, markersize=msize)

arrow_props = dict(arrowstyle="-", connectionstyle="bar")

# fast and slow
plt.annotate('', xy=(0.95, wt_1yta+6), xytext=(3.05, wt_1yta+6), arrowprops=arrow_props)
plt.annotate('****', xy=(1.75, wt_1yta+10.2), xytext=(1.62, wt_1yta+10.2), fontsize=8)

# wt and slow
plt.annotate('', xy=(2.05, wt_1yta+3.0), xytext=(2.95, wt_1yta+3.0), arrowprops=arrow_props)
plt.annotate('**', xy=(2.09, slow_1yta+10.5), xytext=(2.35, slow_1yta+10.5), fontsize=8)

plt.xticks([1,2,3], ['WT', 'FAST', 'SLOW'], fontsize=fontsize)
plt.yticks([-45.0, -50.0, -55.0, -60.0, -65.0, -70.0], [-45.0, -50.0, -55.0, -60.0, -65.0, -70.0],  fontsize=fontsize)
plt.xlim(0.3, 3.7)
plt.ylim(-70, -45)

plt.ylabel('Interface interaction energy,\nkcal/mol', fontsize=fontsize)
plt.tight_layout()

plt.savefig('1yta_interface_energy.png', dpi=300)

plt.clf()

ax = plt.subplot()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(spine_width)
ax.spines['bottom'].set_linewidth(spine_width)
ax.xaxis.set_tick_params(width=tick_width)
ax.yaxis.set_tick_params(width=tick_width)
plt.errorbar(1,   wt_2is3, marker='o', label=  'WT', color=colors[1], alpha=alpha, yerr=  wt_2is3_yerr, capsize=capsize, ecolor='k', elinewidth=ewidth, markersize=msize)
plt.errorbar(2, fast_2is3, marker='o', label='FAST', color=colors[2], alpha=alpha, yerr=fast_2is3_yerr, capsize=capsize, ecolor='k', elinewidth=ewidth, markersize=msize)
plt.errorbar(3, slow_2is3, marker='o', label='SLOW', color=colors[0], alpha=alpha, yerr=slow_2is3_yerr, capsize=capsize, ecolor='k', elinewidth=ewidth, markersize=msize)

plt.ylabel('Interface interaction energy,\n kcal/mol', fontsize=fontsize)
plt.xticks([1,2,3], ['WT', 'FAST', 'SLOW'], fontsize=fontsize)
plt.yticks([-30.0, -32.5, -35.0, -37.5, -40], [-30.0, -32.5, -35.0, -37.5, -40.0],  fontsize=fontsize)
plt.xlim(0.3, 3.7)
plt.ylim(-40, -30)

#plt.xlim(0, 4)
#plt.legend()
plt.tight_layout()
plt.savefig('2is3_interface_energy.png', dpi=300)