#!/usr/bin/env python3

import sys,os
import numpy as np
import glob
from mlxtend.evaluate import bootstrap
import matplotlib.pyplot as plt

if len(sys.argv) != 6:
    print(f'[1] path to MSM file')
    print(f'[2] path to SASA files')
    print(f'[3] path to output file')
    print(f'[4] mask indexed from 0')
    print(f'[5] pdb')
    quit()

#MSM model was build with trajs 1-200 from fast slow wt

##############################################################################################
#DEFS
##############################################################################################
def get_key(fp):
    return int(fp.split('/')[-1].split('_')[0])


def plot_hist(data, state, variant, axs, xlabel='', clear=0):

    #plt.hist(data, 50, range=[min(data), max(data)], density=True, rwidth=1, alpha=0.5)
    hist, bin_edges = np.histogram(data, 20)
    print('\n',hist, bin_edges, len(bin_edges), len(hist), len(data))
    hist = hist / float(len(data))
    print(hist, bin_edges, len(bin_edges), len(hist))
    bin_middles = [bin_edges[i]+(bin_edges[i+1]-bin_edges[i])/2 for i in range(len(bin_edges)) if i+1 < len(bin_edges)]
    print(bin_middles, len(bin_middles))

    bar_width = bin_edges[1] - bin_edges[0]
    axs[state].bar(bin_middles, hist, alpha=0.5, width=bar_width, label=f'{variant}_{state}')
    #axs[state].set_ylim(0,0.3)
    axs[state].legend(loc='upper right')
    if xlabel != '':
        axs[state].set_xlabel(xlabel)
    #plt.savefig(f'Hist_{variant}_state{state}.png')

    return

##############################################################################################
#MAIN
##############################################################################################
inpfiles = sys.argv[1]

#load msm data
pdb = sys.argv[5]
print(f'PDB: {pdb}')
print(f'')
msm = np.load(f'{inpfiles}{pdb}_msm_data.npz', allow_pickle=True)
msm_data = msm['meta_dtrajs']
print(msm_data, msm_data.shape, list(msm.keys()))

#get number of unique meta stable states
num_unique_msm_states = np.unique(msm_data)
print(f'num_unique_msm_states: {num_unique_msm_states}')

#parse mask to residue index starting from 0
mask_str = sys.argv[4]
mask = np.arange(int(sys.argv[4].split('-')[0]), int(sys.argv[4].split('-')[1])+1)
print(mask)


#for each meta stable state calcualte the relative SASA deviation
msm_outdata_base = []
raw_msm_outdata_base = []
percent_exposed_outdata_base = []

##############################################################################################
#load SASA data into same shape array as msm_data (traj, frame, res)
disordered_files = sorted(glob.glob(f'{sys.argv[2]}/native/*native_800k*{pdb}_ca_synth_charmm_sasa.txt*'), key=get_key)
print(disordered_files)
#disordered_data = [np.loadtxt(x, usecols=mask) for x in disordered_files]
disordered_data = []
for x in disordered_files:
    print(x)
    disordered_data.append(np.loadtxt(x, usecols=mask))

disordered_data = np.dstack(disordered_data)
disordered_data = np.rollaxis(disordered_data,-1)
disordered_data = disordered_data.sum(axis=2)
print(disordered_data.shape)
#plot_hist(disordered_data.flatten(), 'N', 'disordered', 'SASA')


#get average disordered data
print(f'\nAnalyzsing disordered state:\n')
disordered_mean = np.mean(disordered_data)
disordered_original, disordered_std_err, disordered_ci_bounds = bootstrap(disordered_data.flatten(), num_rounds=10000, func=np.mean, ci=0.95, seed=np.random.randint(100))
print( disordered_original,disordered_std_err, disordered_ci_bounds )
disordered_var = (np.sqrt(disordered_data.size)* disordered_std_err)**2
lower_delta = disordered_original - disordered_ci_bounds[0]
upper_delta = disordered_ci_bounds[1] - disordered_original
raw_msm_outdata_base.append([ -2 ,disordered_original,  disordered_std_err, disordered_var,  disordered_ci_bounds[0],  disordered_ci_bounds[1],  lower_delta, upper_delta])

##############################################################################################
#load SASA data into same shape array as msm_data (traj, frame, res)
native_files = sorted(glob.glob(f'{sys.argv[2]}/native/*native_{pdb}_ca_synth_charmm_sasa.txt*'), key=get_key)
print(native_files)
native_data = []
for x in native_files:
    print(x)
    native_data.append(np.loadtxt(x, usecols=mask))
native_data = np.stack(native_data)
native_data = native_data.sum(axis=2)
print(native_data.shape)
#plot_hist(native_data.flatten(), 'N', 'Native', 'SASA')


#get average native data
print(f'\nAnalyzsing native state:\n')
native_mean = np.mean(native_data)
native_original, native_std_err, native_ci_bounds = bootstrap(native_data.flatten(), num_rounds=10000, func=np.mean, ci=0.95, seed=np.random.randint(100))
print( native_original,native_std_err, native_ci_bounds )
native_var = (np.sqrt(native_data.size)* native_std_err)**2
lower_delta = native_original - native_ci_bounds[0]
upper_delta = native_ci_bounds[1] - native_original
raw_msm_outdata_base.append([ -1 ,native_original,  native_std_err, native_var,  native_ci_bounds[0],  native_ci_bounds[1],  lower_delta, upper_delta])

native_percent_exposed_data = ((native_data / disordered_original))*100
original, std_err, ci_bounds = bootstrap(native_percent_exposed_data.flatten(), num_rounds=10000, func=np.mean, ci=0.95, seed=np.random.randint(100))
var = (np.sqrt(native_data.size)*std_err)**2
lower_delta = original - ci_bounds[0]
upper_delta = ci_bounds[1] - original
#print( original, std_err, ci_bounds )
percent_exposed_outdata_base.append([ -1 ,original, std_err, var, ci_bounds[0], ci_bounds[1], lower_delta, upper_delta])


##############################################################################################
#for each meta stable state calcualte the relative SASA deviation
msm_outdata = []
raw_msm_outdata = []
percent_exposed_outdata = []
METS_outdata = []

#get master equation stead state probability
METS_data = np.loadtxt(f'{inpfiles}/fast_METS.dat')[:,1:]
print(METS_data, METS_data.shape)

fast_files = sorted(glob.glob(f'{sys.argv[2]}/synthesis/*fast*{pdb}_ca_synth_charmm_sasa.txt*'), key=get_key)
#print(fast_files, len(fast_files))
fast_data = []
for x in fast_files:
    print(x)
    fast_data.append(np.loadtxt(x, usecols=mask))
fast_data = np.dstack(fast_data)
fast_data = np.rollaxis(fast_data,-1)
fast_data = fast_data.sum(axis=2)
#print(fast_data, fast_data.shape)

fast_fig, fast_axs = plt.subplots(len(num_unique_msm_states)+1, sharex=True)
for msm_state in num_unique_msm_states:
    print(f'\nAnalyzsing FAST meta stable state: {msm_state}\n')
    idxs = np.where(msm_data[0:200,:] == msm_state)

    raw_fast_msm_sasa_data = fast_data[idxs[0],idxs[1]]
    plot_hist(raw_fast_msm_sasa_data.flatten(), msm_state, 'fast', fast_axs)
    original, std_err, ci_bounds = bootstrap(raw_fast_msm_sasa_data.flatten(), num_rounds=10000, func=np.mean, ci=0.95, seed=np.random.randint(100))
    var = (np.sqrt(raw_fast_msm_sasa_data.size)*std_err)**2
    #print( original, std_err, ci_bounds )
    lower_delta = original - ci_bounds[0]
    upper_delta = ci_bounds[1] - original
    raw_msm_outdata.append([ msm_state ,original, std_err, var, ci_bounds[0], ci_bounds[1], lower_delta, upper_delta])


    fast_msm_sasa_data = ((fast_data[idxs[0],idxs[1]] / native_original) - 1)*100
    original, std_err, ci_bounds = bootstrap(fast_msm_sasa_data.flatten(), num_rounds=10000, func=np.mean, ci=0.95, seed=np.random.randint(100))
    var = (np.sqrt(fast_msm_sasa_data.size)*std_err)**2
    lower_delta = original - ci_bounds[0]
    upper_delta = ci_bounds[1] - original
    #print( original, std_err, ci_bounds )
    msm_outdata.append([ msm_state ,original, std_err, var, ci_bounds[0], ci_bounds[1], lower_delta, upper_delta])

    fast_msm_percent_exposed_data = ((fast_data[idxs[0],idxs[1]] / disordered_original))*100
    original, std_err, ci_bounds = bootstrap(fast_msm_percent_exposed_data.flatten(), num_rounds=10000, func=np.mean, ci=0.95, seed=np.random.randint(100))
    var = (np.sqrt(fast_msm_sasa_data.size)*std_err)**2
    lower_delta = original - ci_bounds[0]
    upper_delta = ci_bounds[1] - original
    #print( original, std_err, ci_bounds )
    percent_exposed_outdata.append([ msm_state ,original, std_err, var, ci_bounds[0], ci_bounds[1], lower_delta, upper_delta])

    #plot histogram

    #Master equeation state_prob
    original, std_err, ci_bounds = bootstrap(METS_data[500:,msm_state].flatten(), num_rounds=10000, func=np.mean, ci=0.95, seed=np.random.randint(100))
    var = (np.sqrt(METS_data[500:,msm_state].size)*std_err)**2
    lower_delta = original - ci_bounds[0]
    upper_delta = ci_bounds[1] - original
    METS_outdata.append([msm_state ,original, std_err, var, ci_bounds[0], ci_bounds[1], lower_delta, upper_delta])


plt.suptitle = f'fast_SASA_histograms'
plot_hist(native_data.flatten(), -1, 'Native' ,fast_axs, xlabel='SASA')
fast_fig.tight_layout()
fast_fig.savefig(f'{sys.argv[3]}_fast_SASA_histograms_resids_{mask_str}.png')

plt.clf()


msm_outdata = msm_outdata_base + msm_outdata
raw_msm_outdata = raw_msm_outdata_base + raw_msm_outdata
percent_exposed_outdata = percent_exposed_outdata + percent_exposed_outdata_base

#np.savetxt(f'{sys.argv[3]}_zeta_msm_{sys.argv[4]}_SASA.txt',msm_outdata, header=f'state, mean, stderr, ci_lower, ci_upper, variant', fmt='%i, %.8f,%.8f,%.8f,%.8f,%s')
np.savetxt(f'{sys.argv[3]}raw_{pdb}_fast_msm_{sys.argv[4]}_SASA.txt',raw_msm_outdata, header=f'state_raw_fast_{pdb}_{sys.argv[4]}, mean_raw_fast_{pdb}_{sys.argv[4]}, stderr_raw_fast_{pdb}_{sys.argv[4]}, var_raw_fast_{pdb}_{sys.argv[4]}, ci_lower_raw_fast_{pdb}_{sys.argv[4]}, ci_upper_raw_fast_{pdb}_{sys.argv[4]}, lower_delta_raw_fast_{pdb}_{sys.argv[4]}, upper_delta_raw_fast_{pdb}_{sys.argv[4]}', fmt='%s,%s,%s,%s,%s,%s,%s,%s')
print(f'Saved: {sys.argv[3]}raw_{pdb}_fast_msm_{sys.argv[4]}_SASA.txt')

#np.savetxt(f'{sys.argv[3]}_zeta_{pdb}_msm_{sys.argv[4]}_SASA.txt',msm_outdata_fast_{pdb}_{sys.argv[4]}, header=f'state_fast_{pdb}_{sys.argv[4]}, mean_fast_{pdb}_{sys.argv[4]}, stderr_fast_{pdb}_{sys.argv[4]}, ci_lower_fast_{pdb}_{sys.argv[4]}, ci_upper_fast_{pdb}_{sys.argv[4]}, variant'_fast_{pdb}_{sys.argv[4]}, fmt='%i_fast_{pdb}_{sys.argv[4]}, %.8f,%.8f,%.8f,%.8f,%s')
np.savetxt(f'{sys.argv[3]}zeta_{pdb}_fast_msm_{sys.argv[4]}_SASA.txt',msm_outdata, header=f'state_zeta_fast_{pdb}_{sys.argv[4]}, mean_zeta_fast_{pdb}_{sys.argv[4]}, stderr_zeta_fast_{pdb}_{sys.argv[4]}, var_zeta_fast_{pdb}_{sys.argv[4]}, ci_lower_zeta_fast_{pdb}_{sys.argv[4]}, ci_upper_zeta_fast_{pdb}_{sys.argv[4]}, lower_delta_zeta_fast_{pdb}_{sys.argv[4]}, upper_delta_zeta_fast_{pdb}_{sys.argv[4]},  variant_zeta_fast', fmt='%s,%s,%s,%s,%s,%s,%s,%s')
print(f'Saved: {sys.argv[3]}zeta_{pdb}_fast_msm_{sys.argv[4]}_SASA.txt')

#np.savetxt(f'{sys.argv[3]}_zeta_{pdb}_msm_{sys.argv[4]}_SASA.txt',msm_outdata_fast_{pdb}_{sys.argv[4]}, header=f'state_fast_{pdb}_{sys.argv[4]}, mean_fast_{pdb}_{sys.argv[4]}, stderr_fast_{pdb}_{sys.argv[4]}, ci_lower_fast_{pdb}_{sys.argv[4]}, ci_upper_fast_{pdb}_{sys.argv[4]}, variant'_fast_{pdb}_{sys.argv[4]}, fmt='%i_fast_{pdb}_{sys.argv[4]}, %.8f,%.8f,%.8f,%.8f,%s')
np.savetxt(f'{sys.argv[3]}percent_exposed_{pdb}_fast_msm_{sys.argv[4]}_SASA.txt',percent_exposed_outdata, header=f'state_perce_fast_{pdb}_{sys.argv[4]}, mean_perce_fast_{pdb}_{sys.argv[4]}, stderr_perce_fast_{pdb}_{sys.argv[4]}, var_perce_fast_{pdb}_{sys.argv[4]}, ci_lower_perce_fast_{pdb}_{sys.argv[4]}, ci_upper_perce_fast_{pdb}_{sys.argv[4]}, lower_delta_perce_fast_{pdb}_{sys.argv[4]}, upper_delta_perce_fast_{pdb}_{sys.argv[4]}', fmt='%s,%s,%s,%s,%s,%s,%s,%s')
print(f'Saved: {sys.argv[3]}percent_exposed_{pdb}_fast_msm_{sys.argv[4]}_SASA.txt')

np.savetxt(f'{sys.argv[3]}METS_{pdb}_fast_msm_{sys.argv[4]}_SASA.txt',METS_outdata, header=f'state_prob_fast_{pdb}_{sys.argv[4]}, mean_prob_fast_{pdb}_{sys.argv[4]}, stderr_prob_fast_{pdb}_{sys.argv[4]}, var_prob_fast_{pdb}_{sys.argv[4]}, ci_lower_prob_fast_{pdb}_{sys.argv[4]}, ci_upper_prob_fast_{pdb}_{sys.argv[4]}, lower_delta_prob_fast_{pdb}_{sys.argv[4]}, upper_delta_prob_fast_{pdb}_{sys.argv[4]}', fmt='%s,%s,%s,%s,%s,%s,%s,%s')
print(f'Saved: {sys.argv[3]}METS_{pdb}_fast_msm_{sys.argv[4]}_SASA.txt')
#np.savetxt(f'{sys.argv[3]}_slow_msm_{sys.argv[4]}_SASA.txt',msm_outdata, header=f'state, mean, stderr, ci_lower, ci_upper, variant')
#print(f'Saved: {sys.argv[3]}_slow_msm_{sys.argv[4]}_SASA.txt')

##############################################################################################
msm_outdata = []
raw_msm_outdata = []
percent_exposed_outdata = []
METS_outdata = []

#get master equation stead state probability
METS_data = np.loadtxt(f'{inpfiles}/slow_METS.dat')[:,1:]
print(METS_data, METS_data.shape)

slow_files = sorted(glob.glob(f'{sys.argv[2]}/synthesis/*slow*{pdb}_ca_synth_charmm_sasa.txt*'), key=get_key)
#print(slow_files, len(slow_files))
slow_data = [np.loadtxt(x, usecols=mask) for x in slow_files]
slow_data = np.dstack(slow_data)
slow_data = np.rollaxis(slow_data,-1)
slow_data = slow_data.sum(axis=2)
#print(slow_data, slow_data.shape)

slow_fig, slow_axs = plt.subplots(len(num_unique_msm_states)+1, sharex=True)
#for each meta stable state calcualte the relative SASA deviation
#msm_outdata = []
for msm_state in num_unique_msm_states:
    print(f'\nAnalyzsing slow meta stable state: {msm_state}\n')
    idxs = np.where(msm_data[200:400,:] == msm_state)

    raw_slow_msm_sasa_data = slow_data[idxs[0],idxs[1]]
    plot_hist(raw_slow_msm_sasa_data.flatten(), msm_state, 'slow', slow_axs)
    original, std_err, ci_bounds = bootstrap(raw_slow_msm_sasa_data.flatten(), num_rounds=10000, func=np.mean, ci=0.95, seed=np.random.randint(100))
    var = (np.sqrt(raw_slow_msm_sasa_data.size)*std_err)**2
    lower_delta = original - ci_bounds[0]
    upper_delta = ci_bounds[1] - original
    #print( original, std_err, ci_bounds )
    raw_msm_outdata.append([ msm_state ,original, std_err, var, ci_bounds[0], ci_bounds[1], lower_delta, upper_delta])

    slow_msm_sasa_data = ((slow_data[idxs[0],idxs[1]] / native_original) - 1)*100
    original, std_err, ci_bounds = bootstrap(slow_msm_sasa_data.flatten(), num_rounds=10000, func=np.mean, ci=0.95, seed=np.random.randint(100))
    var = (np.sqrt(slow_msm_sasa_data.size)*std_err)**2
    lower_delta = original - ci_bounds[0]
    upper_delta = ci_bounds[1] - original
    #print( original, std_err, ci_bounds )
    msm_outdata.append([ msm_state ,original, std_err, var, ci_bounds[0], ci_bounds[1], lower_delta, upper_delta])

    slow_msm_percent_exposed_data = ((slow_data[idxs[0],idxs[1]] / disordered_original))*100
    original, std_err, ci_bounds = bootstrap(slow_msm_percent_exposed_data.flatten(), num_rounds=10000, func=np.mean, ci=0.95, seed=np.random.randint(100))
    var = (np.sqrt(slow_msm_sasa_data.size)*std_err)**2
    lower_delta = original - ci_bounds[0]
    upper_delta = ci_bounds[1] - original
    #print( original, std_err, ci_bounds )
    percent_exposed_outdata.append([ msm_state ,original, std_err, var, ci_bounds[0], ci_bounds[1], lower_delta, upper_delta])

    #Master equeation state_prob
    original, std_err, ci_bounds = bootstrap(METS_data[500:,msm_state].flatten(), num_rounds=10000, func=np.mean, ci=0.95, seed=np.random.randint(100))
    var = (np.sqrt(METS_data[500:,msm_state].size)*std_err)**2
    lower_delta = original - ci_bounds[0]
    upper_delta = ci_bounds[1] - original
    METS_outdata.append([msm_state ,original, std_err, var, ci_bounds[0], ci_bounds[1], lower_delta, upper_delta])

plt.suptitle = f'slow_SASA_histograms'
plot_hist(native_data.flatten(), -1, 'Native' ,slow_axs, xlabel='SASA')
slow_fig.tight_layout()
slow_fig.savefig(f'{sys.argv[3]}_slow_SASA_histograms_resids_{mask_str}.png')

plt.clf()

#np.savetxt(f'{sys.argv[3]}_slow_msm_{sys.argv[4]}_SASA.txt',msm_outdata)
#print(f'Saved: {sys.argv[3]}_slow_msm_{sys.argv[4]}_SASA.txt')
msm_outdata = msm_outdata_base + msm_outdata
raw_msm_outdata = raw_msm_outdata_base + raw_msm_outdata
percent_exposed_outdata = percent_exposed_outdata + percent_exposed_outdata_base


np.savetxt(f'{sys.argv[3]}_raw_{pdb}_slow_msm_{sys.argv[4]}_SASA.txt',raw_msm_outdata, header=f'state_raw_slow_{pdb}_{sys.argv[4]}, mean_raw_slow_{pdb}_{sys.argv[4]}, stderr_raw_slow_{pdb}_{sys.argv[4]}, var_raw_slow_{pdb}_{sys.argv[4]}, ci_lower_raw_slow_{pdb}_{sys.argv[4]}, ci_upper_raw_slow_{pdb}_{sys.argv[4]}, lower_delta_raw_slow_{pdb}_{sys.argv[4]}, upper_delta_raw_slow_{pdb}_{sys.argv[4]}', fmt='%s,%s,%s,%s,%s,%s,%s,%s')
print(f'Saved: {sys.argv[3]}_raw_{pdb}_slow_msm_{sys.argv[4]}_SASA.txt')

#np.savetxt(f'{sys.argv[3]}_zeta_{pdb}_msm_{sys.argv[4]}_SASA.txt',msm_outdata_slow_{pdb}_{sys.argv[4]}, header=f'state_slow_{pdb}_{sys.argv[4]}, mean_slow_{pdb}_{sys.argv[4]}, stderr_slow_{pdb}_{sys.argv[4]}, ci_lower_slow_{pdb}_{sys.argv[4]}, ci_upper_slow_{pdb}_{sys.argv[4]}, variant'_slow_{pdb}_{sys.argv[4]}, fmt='%i_slow_{pdb}_{sys.argv[4]}, %.8f,%.8f,%.8f,%.8f,%s')
np.savetxt(f'{sys.argv[3]}_zeta_{pdb}_slow_msm_{sys.argv[4]}_SASA.txt',msm_outdata, header=f'state_zeta_slow_{pdb}_{sys.argv[4]}, mean_zeta_slow_{pdb}_{sys.argv[4]}, stderr_zeta_slow_{pdb}_{sys.argv[4]}, var_zeta_slow_{pdb}_{sys.argv[4]}, ci_lower_zeta_slow_{pdb}_{sys.argv[4]}, ci_upper_zeta_slow_{pdb}_{sys.argv[4]}, lower_delta_zeta_slow_{pdb}_{sys.argv[4]}, upper_delta_zeta_slow_{pdb}_{sys.argv[4]},  variant_zeta_slow', fmt='%s,%s,%s,%s,%s,%s,%s,%s')
print(f'Saved: {sys.argv[3]}_zeta_{pdb}_slow_msm_{sys.argv[4]}_SASA.txt')

#np.savetxt(f'{sys.argv[3]}_zeta_{pdb}_msm_{sys.argv[4]}_SASA.txt',msm_outdata_slow_{pdb}_{sys.argv[4]}, header=f'state_slow_{pdb}_{sys.argv[4]}, mean_slow_{pdb}_{sys.argv[4]}, stderr_slow_{pdb}_{sys.argv[4]}, ci_lower_slow_{pdb}_{sys.argv[4]}, ci_upper_slow_{pdb}_{sys.argv[4]}, variant'_slow_{pdb}_{sys.argv[4]}, fmt='%i_slow_{pdb}_{sys.argv[4]}, %.8f,%.8f,%.8f,%.8f,%s')
np.savetxt(f'{sys.argv[3]}_percent_exposed_{pdb}_slow_msm_{sys.argv[4]}_SASA.txt',percent_exposed_outdata, header=f'state_perce_slow_{pdb}_{sys.argv[4]}, mean_perce_slow_{pdb}_{sys.argv[4]}, stderr_perce_slow_{pdb}_{sys.argv[4]}, var_perce_slow_{pdb}_{sys.argv[4]}, ci_lower_perce_slow_{pdb}_{sys.argv[4]}, ci_upper_perce_slow_{pdb}_{sys.argv[4]}, lower_delta_perce_slow_{pdb}_{sys.argv[4]}, upper_delta_perce_slow_{pdb}_{sys.argv[4]}', fmt='%s,%s,%s,%s,%s,%s,%s,%s')
print(f'Saved: {sys.argv[3]}_percent_exposed_{pdb}_slow_msm_{sys.argv[4]}_SASA.txt')

np.savetxt(f'{sys.argv[3]}_METS_{pdb}_slow_msm_{sys.argv[4]}_SASA.txt',METS_outdata, header=f'state_prob_slow_{pdb}_{sys.argv[4]}, mean_prob_slow_{pdb}_{sys.argv[4]}, stderr_prob_slow_{pdb}_{sys.argv[4]}, var_prob_slow_{pdb}_{sys.argv[4]}, ci_lower_prob_slow_{pdb}_{sys.argv[4]}, ci_upper_prob_slow_{pdb}_{sys.argv[4]}, lower_delta_prob_slow_{pdb}_{sys.argv[4]}, upper_delta_prob_slow_{pdb}_{sys.argv[4]}', fmt='%s,%s,%s,%s,%s,%s,%s,%s')
print(f'Saved: {sys.argv[3]}_METS_{pdb}_slow_msm_{sys.argv[4]}_SASA.txt')


##############################################################################################
msm_outdata = []
raw_msm_outdata = []
percent_exposed_outdata = []
METS_outdata = []

#get master equation stead state probability
METS_data = np.loadtxt(f'{inpfiles}/wt_METS.dat')[:,1:]
print(METS_data, METS_data.shape)

wt_files = sorted(glob.glob(f'{sys.argv[2]}/synthesis/*wt*{pdb}_ca_synth_charmm_sasa.txt*'), key=get_key)
#print(wt_files, len(wt_files))
wt_data = [np.loadtxt(x, usecols=mask) for x in wt_files]
wt_data = np.dstack(wt_data)
wt_data = np.rollaxis(wt_data,-1)
wt_data = wt_data.sum(axis=2)
#print(wt_data, wt_data.shape)

wt_fig, wt_axs = plt.subplots(len(num_unique_msm_states)+1, sharex=True)
#for each meta stable state calcualte the relative SASA deviation
#msm_outdata = []
for msm_state in num_unique_msm_states:
    print(f'\nAnalyzsing wt meta stable state: {msm_state}\n')
    idxs = np.where(msm_data[400:600,:] == msm_state)

    raw_wt_msm_sasa_data = wt_data[idxs[0],idxs[1]]
    plot_hist(raw_wt_msm_sasa_data.flatten(), msm_state, 'wt', wt_axs)
    original, std_err, ci_bounds = bootstrap(raw_wt_msm_sasa_data.flatten(), num_rounds=10000, func=np.mean, ci=0.95, seed=np.random.randint(100))
    var = (np.sqrt(raw_wt_msm_sasa_data.size)*std_err)**2
    lower_delta = original - ci_bounds[0]
    upper_delta = ci_bounds[1] - original
    #print( original, std_err, ci_bounds )
    raw_msm_outdata.append([ msm_state ,original, std_err, var, ci_bounds[0], ci_bounds[1], lower_delta, upper_delta])

    wt_msm_sasa_data = ((wt_data[idxs[0],idxs[1]] / native_original) - 1)*100
    original, std_err, ci_bounds = bootstrap(wt_msm_sasa_data.flatten(), num_rounds=10000, func=np.mean, ci=0.95, seed=np.random.randint(100))
    var = (np.sqrt(wt_msm_sasa_data.size)*std_err)**2
    lower_delta = original - ci_bounds[0]
    upper_delta = ci_bounds[1] - original
    #print( original, std_err, ci_bounds )
    msm_outdata.append([ msm_state ,original, std_err, var, ci_bounds[0], ci_bounds[1], lower_delta, upper_delta])

    wt_msm_percent_exposed_data = ((wt_data[idxs[0],idxs[1]] / disordered_original))*100
    original, std_err, ci_bounds = bootstrap(wt_msm_percent_exposed_data.flatten(), num_rounds=10000, func=np.mean, ci=0.95, seed=np.random.randint(100))
    var = (np.sqrt(wt_msm_sasa_data.size)*std_err)**2
    lower_delta = original - ci_bounds[0]
    upper_delta = ci_bounds[1] - original
    #print( original, std_err, ci_bounds )
    percent_exposed_outdata.append([ msm_state ,original, std_err, var, ci_bounds[0], ci_bounds[1], lower_delta, upper_delta])

    #Master equeation state_prob
    original, std_err, ci_bounds = bootstrap(METS_data[500:,msm_state].flatten(), num_rounds=10000, func=np.mean, ci=0.95, seed=np.random.randint(100))
    var = (np.sqrt(METS_data[500:,msm_state].size)*std_err)**2
    lower_delta = original - ci_bounds[0]
    upper_delta = ci_bounds[1] - original
    METS_outdata.append([msm_state ,original, std_err, var, ci_bounds[0], ci_bounds[1], lower_delta, upper_delta])

plt.suptitle = f'wt_SASA_histograms'
plot_hist(native_data.flatten(), -1, 'Native' ,wt_axs, xlabel='SASA')
wt_fig.tight_layout()
wt_fig.savefig(f'{sys.argv[3]}_wt_SASA_histograms_resids_{mask_str}.png')

plt.clf()

#np.savetxt(f'{sys.argv[3]}_wt_msm_{sys.argv[4]}_SASA.txt',msm_outdata)
#print(f'Saved: {sys.argv[3]}_wt_msm_{sys.argv[4]}_SASA.txt')
msm_outdata = msm_outdata_base + msm_outdata
raw_msm_outdata = raw_msm_outdata_base + raw_msm_outdata
percent_exposed_outdata = percent_exposed_outdata + percent_exposed_outdata_base


np.savetxt(f'{sys.argv[3]}_raw_{pdb}_wt_msm_{sys.argv[4]}_SASA.txt',raw_msm_outdata, header=f'state_raw_wt_{pdb}_{sys.argv[4]}, mean_raw_wt_{pdb}_{sys.argv[4]}, stderr_raw_wt_{pdb}_{sys.argv[4]}, var_raw_wt_{pdb}_{sys.argv[4]}, ci_lower_raw_wt_{pdb}_{sys.argv[4]}, ci_upper_raw_wt_{pdb}_{sys.argv[4]}, lower_delta_raw_wt_{pdb}_{sys.argv[4]}, upper_delta_raw_wt_{pdb}_{sys.argv[4]}', fmt='%s,%s,%s,%s,%s,%s,%s,%s')
print(f'Saved: {sys.argv[3]}_raw_{pdb}_wt_msm_{sys.argv[4]}_SASA.txt')

#np.savetxt(f'{sys.argv[3]}_zeta_{pdb}_msm_{sys.argv[4]}_SASA.txt',msm_outdata_wt_{pdb}_{sys.argv[4]}, header=f'state_wt_{pdb}_{sys.argv[4]}, mean_wt_{pdb}_{sys.argv[4]}, stderr_wt_{pdb}_{sys.argv[4]}, ci_lower_wt_{pdb}_{sys.argv[4]}, ci_upper_wt_{pdb}_{sys.argv[4]}, variant'_wt_{pdb}_{sys.argv[4]}, fmt='%i_wt_{pdb}_{sys.argv[4]}, %.8f,%.8f,%.8f,%.8f,%s')
np.savetxt(f'{sys.argv[3]}_zeta_{pdb}_wt_msm_{sys.argv[4]}_SASA.txt',msm_outdata, header=f'state_zeta_wt_{pdb}_{sys.argv[4]}, mean_zeta_wt_{pdb}_{sys.argv[4]}, stderr_zeta_wt_{pdb}_{sys.argv[4]}, var_zeta_wt_{pdb}_{sys.argv[4]}, ci_lower_zeta_wt_{pdb}_{sys.argv[4]}, ci_upper_zeta_wt_{pdb}_{sys.argv[4]}, lower_delta_zeta_wt_{pdb}_{sys.argv[4]}, upper_delta_zeta_wt_{pdb}_{sys.argv[4]},  variant_zeta_wt', fmt='%s,%s,%s,%s,%s,%s,%s,%s')
print(f'Saved: {sys.argv[3]}_zeta_{pdb}_wt_msm_{sys.argv[4]}_SASA.txt')

#np.savetxt(f'{sys.argv[3]}_zeta_{pdb}_msm_{sys.argv[4]}_SASA.txt',msm_outdata_wt_{pdb}_{sys.argv[4]}, header=f'state_wt_{pdb}_{sys.argv[4]}, mean_wt_{pdb}_{sys.argv[4]}, stderr_wt_{pdb}_{sys.argv[4]}, ci_lower_wt_{pdb}_{sys.argv[4]}, ci_upper_wt_{pdb}_{sys.argv[4]}, variant'_wt_{pdb}_{sys.argv[4]}, fmt='%i_wt_{pdb}_{sys.argv[4]}, %.8f,%.8f,%.8f,%.8f,%s')
np.savetxt(f'{sys.argv[3]}_percent_exposed_{pdb}_wt_msm_{sys.argv[4]}_SASA.txt',percent_exposed_outdata, header=f'state_perce_wt_{pdb}_{sys.argv[4]}, mean_perce_wt_{pdb}_{sys.argv[4]}, stderr_perce_wt_{pdb}_{sys.argv[4]}, var_perce_wt_{pdb}_{sys.argv[4]}, ci_lower_perce_wt_{pdb}_{sys.argv[4]}, ci_upper_perce_wt_{pdb}_{sys.argv[4]}, lower_delta_perce_wt_{pdb}_{sys.argv[4]}, upper_delta_perce_wt_{pdb}_{sys.argv[4]}', fmt='%s,%s,%s,%s,%s,%s,%s,%s')
print(f'Saved: {sys.argv[3]}_percent_exposed_{pdb}_wt_msm_{sys.argv[4]}_SASA.txt')

np.savetxt(f'{sys.argv[3]}_METS_{pdb}_wt_msm_{sys.argv[4]}_SASA.txt',METS_outdata, header=f'state_prob_wt_{pdb}_{sys.argv[4]}, mean_prob_wt_{pdb}_{sys.argv[4]}, stderr_prob_wt_{pdb}_{sys.argv[4]}, var_prob_wt_{pdb}_{sys.argv[4]}, ci_lower_prob_wt_{pdb}_{sys.argv[4]}, ci_upper_prob_wt_{pdb}_{sys.argv[4]}, lower_delta_prob_wt_{pdb}_{sys.argv[4]}, upper_delta_prob_wt_{pdb}_{sys.argv[4]}', fmt='%s,%s,%s,%s,%s,%s,%s,%s')
print(f'Saved: {sys.argv[3]}_METS_{pdb}_wt_msm_{sys.argv[4]}_SASA.txt')
