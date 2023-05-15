# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 14:56:55 2021

@author: nissley
"""

import os, sys
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats.stats import pearsonr
import matplotlib.pyplot as plt

from matplotlib.pyplot import figure
figure(num=None, figsize=(6,6), dpi=300, facecolor='w', edgecolor='k')

# the fitting function - note that we constrain the two coefficient to sum to unity
def func(t, k1, k2, f1):
    return (f1*np.exp(-k1*t))+((1.0-f1)*np.exp(-k2*t))

### 1YTA / oligoribonuclease WT

data = np.loadtxt('fpt_analysis/fpt_analysis_1yta_wt/frac_U_vs_time_200trajs.txt')
times = data[:,1]
surv_prob_U = data[:,2]

# do the fitting
sigma = np.ones(len(times))
sigma[0] = 0.01 # gives higher weights to first point, so it will go almost exactly through it
popt, pcov = curve_fit(func, times, surv_prob_U, p0=[1E-5, 1E-8, 0.2], bounds=(0., [np.inf, np.inf, 1.]))
k1 = popt[0]
k2 = popt[1]
f1 = popt[2]
f2 = 1.0 - popt[2]

print ('f1:', f1)
print ('k1:', k1)
print ('t1:', 1./k1)
print ('f2:', f2)
print ('k2:', k2)
print ('t2:', 1./k2)

fit = (f1*np.exp(-k1*times))+((f2)*np.exp(-k2*times))
r = pearsonr(fit, surv_prob_U)[0]
print (r**2.)

plt.clf()

ax = plt.subplot(321)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.plot(data[:,1], data[:,2], color='dodgerblue', alpha=0.7, linewidth=3.0, label='Simulation')
plt.plot(times, fit, color='magenta', alpha=0.7, linewidth=1.0, linestyle='--', label='Fit')
plt.xlim(-0.1, 5.1)
plt.xticks([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
plt.ylim(-0.05, 1.05)
#plt.xlabel('Time, '+r'$\mu$s')
plt.legend()
plt.ylabel(r'$S_{\rm U}(t)$')

### 1YTA / oligoribonuclease FAST

data = np.loadtxt('fpt_analysis/fpt_analysis_1yta_fast/frac_U_vs_time_200trajs.txt')
times = data[:,1]
surv_prob_U = data[:,2]

# do the fitting
sigma = np.ones(len(times))
sigma[0] = 0.01 # gives higher weights to first point, so it will go almost exactly through it
popt, pcov = curve_fit(func, times, surv_prob_U, p0=[1E-5, 1E-8, 0.2], bounds=(0., [np.inf, np.inf, 1.]))
k1 = popt[0]
k2 = popt[1]
f1 = popt[2]
f2 = 1.0 - popt[2]

print ('f1:', f1)
print ('k1:', k1)
print ('t1:', 1./k1)
print ('f2:', f2)
print ('k2:', k2)
print ('t2:', 1./k2)

fit = (f1*np.exp(-k1*times))+((f2)*np.exp(-k2*times))
r = pearsonr(fit, surv_prob_U)[0]
print (r**2.)

ax = plt.subplot(323)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.plot(data[:,1], data[:,2], color='dodgerblue', alpha=0.7, linewidth=3.0)
plt.plot(times, fit, color='magenta', alpha=0.7, linewidth=1.0, linestyle='--')
plt.xlim(-0.1, 5.1)
plt.xticks([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
plt.ylim(-0.05, 1.05)
#plt.xlabel('Time, '+r'$\mu$s')
plt.ylabel(r'$S_{\rm U}(t)$')

### 1YTA / oligoribonuclease SLOW

data = np.loadtxt('fpt_analysis/fpt_analysis_1yta_slow/frac_U_vs_time_200trajs.txt')
times = data[:,1]
surv_prob_U = data[:,2]

# do the fitting
sigma = np.ones(len(times))
sigma[0] = 0.01 # gives higher weights to first point, so it will go almost exactly through it
popt, pcov = curve_fit(func, times, surv_prob_U, p0=[1E-5, 1E-8, 0.2], bounds=(0., [np.inf, np.inf, 1.]))
k1 = popt[0]
k2 = popt[1]
f1 = popt[2]
f2 = 1.0 - popt[2]

print ('f1:', f1)
print ('k1:', k1)
print ('t1:', 1./k1)
print ('f2:', f2)
print ('k2:', k2)
print ('t2:', 1./k2)

fit = (f1*np.exp(-k1*times))+((f2)*np.exp(-k2*times))
r = pearsonr(fit, surv_prob_U)[0]
print (r**2.)

ax = plt.subplot(325)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.plot(data[:,1], data[:,2], color='dodgerblue', alpha=0.7, linewidth=3.0)
plt.plot(times, fit, color='magenta', alpha=0.7, linewidth=1.0, linestyle='--')
plt.xlim(-0.1, 5.1)
plt.xticks([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
plt.ylim(-0.05, 1.05)
plt.xlabel('Time, '+r'$\mu$s')
plt.ylabel(r'$S_{\rm U}(t)$')

#### 2IS3 / ribonuclease T WT

data = np.loadtxt('fpt_analysis/fpt_analysis_2is3_wt/frac_U_vs_time_200trajs.txt')
times = data[:,1]
surv_prob_U = data[:,2]

# do the fitting
sigma = np.ones(len(times))
sigma[0] = 0.01 # gives higher weights to first point, so it will go almost exactly through it
popt, pcov = curve_fit(func, times, surv_prob_U, p0=[1E-5, 1E-8, 0.2], bounds=(0., [np.inf, np.inf, 1.]))
k1 = popt[0]
k2 = popt[1]
f1 = popt[2]
f2 = 1.0 - popt[2]

print ('f1:', f1)
print ('k1:', k1)
print ('t1:', 1./k1)
print ('f2:', f2)
print ('k2:', k2)
print ('t2:', 1./k2)

fit = (f1*np.exp(-k1*times))+((f2)*np.exp(-k2*times))
r = pearsonr(fit, surv_prob_U)[0]
print (r**2.)

ax = plt.subplot(322)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.plot(data[:,1], data[:,2], color='dodgerblue', alpha=0.7, linewidth=3.0)
plt.plot(times, fit, color='magenta', alpha=0.7, linewidth=1.0, linestyle='--')
plt.xlim(-0.02, 0.82)
plt.xticks([0.0, 0.2, 0.4, 0.6, 0.8])
plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1], [])
plt.ylim(-0.05, 1.05)
#plt.xlabel('Time, '+r'$\mu$s')

#### 2IS3 / ribonuclease T FAST

data = np.loadtxt('fpt_analysis/fpt_analysis_2is3_fast/frac_U_vs_time_200trajs.txt')
times = data[:,1]
surv_prob_U = data[:,2]

# do the fitting
sigma = np.ones(len(times))
sigma[0] = 0.01 # gives higher weights to first point, so it will go almost exactly through it
popt, pcov = curve_fit(func, times, surv_prob_U, p0=[1E-5, 1E-8, 0.2], bounds=(0., [np.inf, np.inf, 1.]))
k1 = popt[0]
k2 = popt[1]
f1 = popt[2]
f2 = 1.0 - popt[2]

print ('f1:', f1)
print ('k1:', k1)
print ('t1:', 1./k1)
print ('f2:', f2)
print ('k2:', k2)
print ('t2:', 1./k2)

fit = (f1*np.exp(-k1*times))+((f2)*np.exp(-k2*times))
r = pearsonr(fit, surv_prob_U)[0]
print (r**2.)

ax = plt.subplot(324)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.plot(data[:,1], data[:,2], color='dodgerblue', alpha=0.7, linewidth=3.0)
plt.plot(times, fit, color='magenta', alpha=0.7, linewidth=1.0, linestyle='--')
plt.xlim(-0.02, 0.5)
plt.xticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1], [])
plt.ylim(-0.05, 1.05)
#plt.xlabel('Time, '+r'$\mu$s')

#### 2IS3 / ribonuclease T SLOW

data = np.loadtxt('fpt_analysis/fpt_analysis_2is3_slow/frac_U_vs_time_200trajs.txt')
times = data[:,1]
surv_prob_U = data[:,2]

# do the fitting
sigma = np.ones(len(times))
sigma[0] = 0.01 # gives higher weights to first point, so it will go almost exactly through it
popt, pcov = curve_fit(func, times, surv_prob_U, p0=[1E-5, 1E-8, 0.2], bounds=(0., [np.inf, np.inf, 1.]))
k1 = popt[0]
k2 = popt[1]
f1 = popt[2]
f2 = 1.0 - popt[2]

print ('f1:', f1)
print ('k1:', k1)
print ('t1:', 1./k1)
print ('f2:', f2)
print ('k2:', k2)
print ('t2:', 1./k2)

fit = (f1*np.exp(-k1*times))+((f2)*np.exp(-k2*times))
r = pearsonr(fit, surv_prob_U)[0]
print (r**2.)

ax = plt.subplot(326)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.plot(data[:,1], data[:,2], color='dodgerblue', alpha=0.7, linewidth=3.0)
plt.plot(times, fit, color='magenta', alpha=0.7, linewidth=1.0, linestyle='--')
plt.xlim(-0.02, 0.5)
plt.xticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1], [])
plt.ylim(-0.05, 1.05)
plt.xlabel('Time, '+r'$\mu$s')

# final stuff and making the plot
#plt.tight_layout()
plt.subplots_adjust(hspace=0.5)
plt.savefig('figure_s1.png', dpi=300)
