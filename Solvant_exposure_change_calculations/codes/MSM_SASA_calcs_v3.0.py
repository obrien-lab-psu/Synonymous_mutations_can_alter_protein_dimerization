#!/usr/bin/env python3

import sys,os
import numpy as np
import glob
from mlxtend.evaluate import bootstrap
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact, chi2_contingency

if len(sys.argv) != 6:
    print(f'[1] path to MSM file')
    print(f'[2] path to SASA files')
    print(f'[3] path to output file')
    print(f'[4] mask indexed from 0')
    print(f'[5] pdb')
    quit()

#MSM model was build with trajs 1-200 from fast slow wt

script_title = 'MSM_SASA_calcs'
script_version = 3.0
out_path = sys.argv[3]
print(f'out_path: {out_path}')

##################################################################################################################
### START initial loading of structure files and qualtiy control ###
##################################################################################################################

### START dir declaration ###

if os.path.exists(f'{out_path}/'):
    print(f'{out_path}/ does exists and will be used')
    pass
else:
    os.makedirs(f'{out_path}/')

if os.path.exists(f'{out_path}{script_title}_{script_version}/'):
    print(f'{out_path}{script_title}_{script_version}/ does exists and will be used')
    pass
else:
    os.makedirs(f'{out_path}{script_title}_{script_version}/')

if os.path.exists(f'{out_path}{script_title}_{script_version}/logs/'):
    print(f'{out_path}{script_title}_{script_version}/logs/ does exists and will be used')
    pass
else:
    os.makedirs(f'{out_path}{script_title}_{script_version}/logs/')

if os.path.exists(f'{out_path}{script_title}_{script_version}/output/'):
    print(f'{out_path}{script_title}_{script_version}/output/ does exists and will be used')
    pass
else:
    os.makedirs(f'{out_path}{script_title}_{script_version}/output/')

### END dir declaration ###

### START preference setting ###


np.set_printoptions(precision=4, suppress=True, linewidth=sys.maxsize, threshold=sys.maxsize)
np.seterr(divide='ignore')

### END preference setting ###

######################################################################################################################
# USER DEFINED FUNCTIONS                                                                                             #
######################################################################################################################
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


def ctable_anal(change_data, op_data, state_idxs, native_mean, tail):
    print(change_data)
    print(op_data)
    state_idxs = np.hstack((state_idxs[0][:,None]+1, state_idxs[1][:,None]))
    print(state_idxs)


    outdata = []
    for cutoff in np.linspace(0,100, num=101):
        #[[change + upsasa, upsasa],
        # [change, NEITHER]]
        c_table = np.zeros((2,2))

        print(cutoff)
        for i,(traj, frame) in enumerate(state_idxs):

            #determine if a change was present 1 or not 0

            change_flag = (frame in change_data[traj])

            diff_flag = (op_data[i] - native_mean > cutoff)

            if change_flag and diff_flag:
                c_table[0,0] += 1
                continue

            if not change_flag and not diff_flag:
                c_table[1,1] += 1
                continue

            if not change_flag and diff_flag:
                c_table[1,0] += 1
                continue

            if change_flag and not diff_flag:
                c_table[0,1] += 1
                continue

        print(c_table)

        if 0 in c_table:
            outdata.append([cutoff, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan])
            continue

        chi2, chi2p, dof, ex = chi2_contingency(c_table, correction=True)
        log_chi2p = -np.log(chi2p)
        oddsr, fisherp = fisher_exact(c_table, alternative=f'{tail}')
        log_fisherp = -np.log(fisherp)
        outdata.append([cutoff, oddsr, fisherp, log_fisherp, chi2, dof, chi2p, log_chi2p])

    outdata = np.stack(outdata)
    return outdata

##############################################################################################
#MAIN
##############################################################################################

#load msm data
msm = np.load(sys.argv[1], allow_pickle=True)
msm_data = msm['meta_dtrajs']
#print(msm_data, msm_data.shape, list(msm.keys()))

#get number of unique meta stable states
num_unique_msm_states = np.unique(msm_data)
print(f'num_unique_msm_states: {num_unique_msm_states}')

#parse mask to residue index starting from 0
mask_str = sys.argv[4]
if mask_str.isdigit():
    mask = np.arange(int(mask_str), int(mask_str) + 1)
else:
    mask = np.arange(int(sys.argv[4].split('-')[0]), int(sys.argv[4].split('-')[1])+1)
print(mask)

pdb = sys.argv[5]
print(f'PDB: {pdb}')
#for each meta stable state calcualte the relative SASA deviation
msm_outdata_base = []
raw_msm_outdata_base = []
percent_exposed_outdata_base = []

##############################################################################################
#load SASA data into same shape array as msm_data (traj, frame, res)
print(f'\nAnalyzsing disordered state:\n')
disordered_files = sorted(glob.glob(f'{sys.argv[2]}/native/*native_800k*{pdb}_ca_synth_charmm_sasa.txt*'), key=get_key)
print(disordered_files)
#disordered_data = [np.loadtxt(x, usecols=mask) for x in disordered_files]
disordered_data = []
for x in disordered_files:
    print(x)
    disordered_data.append(np.loadtxt(x, usecols=mask))

disordered_data = np.dstack(disordered_data)
disordered_data = np.rollaxis(disordered_data,-1)
if len(mask) == 1:
    disordered_data = disordered_data.sum(axis=1)
else:
    disordered_data = disordered_data.sum(axis=2)
print(disordered_data.shape)
#plot_hist(disordered_data.flatten(), 'N', 'disordered', 'SASA')

#get average disordered data
disordered_mean = np.mean(disordered_data)
disordered_original, disordered_std_err, disordered_ci_bounds = bootstrap(disordered_data.flatten(),  num_rounds=100000, func=np.mean, ci=0.95, seed=np.random.randint(100))
print( disordered_original,disordered_std_err, disordered_ci_bounds )
disordered_var = (np.sqrt(disordered_data.size)* disordered_std_err)**2
lower_delta = disordered_original - disordered_ci_bounds[0]
upper_delta = disordered_ci_bounds[1] - disordered_original
raw_msm_outdata_base.append([ -2 ,disordered_original,  disordered_std_err, disordered_var,  disordered_ci_bounds[0],  disordered_ci_bounds[1],  lower_delta, upper_delta])

##############################################################################################
#load SASA data into same shape array as msm_data (traj, frame, res)
print(f'\nAnalyzsing native state:\n')
native_files = sorted(glob.glob(f'{sys.argv[2]}/native/*native_{pdb}_ca_synth_charmm_sasa.txt*'), key=get_key)
print(native_files)
native_data = []
for x in native_files:
    print(x)
    native_data.append(np.loadtxt(x, usecols=mask))
native_data = np.dstack(native_data)
native_data = np.rollaxis(native_data,-1)
if len(mask) == 1:
    native_data = native_data.sum(axis=1)
else:
    native_data = native_data.sum(axis=2)
print(native_data.shape)
#plot_hist(native_data.flatten(), 'N', 'Native', 'SASA')

#get average native data
native_mean = np.mean(native_data)
native_original, native_std_err, native_ci_bounds = bootstrap(native_data.flatten(),  num_rounds=100000, func=np.mean, ci=0.95, seed=np.random.randint(100))
print( native_original,native_std_err, native_ci_bounds )
native_var = (np.sqrt(native_data.size)* native_std_err)**2
lower_delta = native_original - native_ci_bounds[0]
upper_delta = native_ci_bounds[1] - native_original
raw_msm_outdata_base.append([ -1 ,native_original,  native_std_err, native_var,  native_ci_bounds[0],  native_ci_bounds[1],  lower_delta, upper_delta])

native_percent_exposed_data = ((native_data / disordered_original))*100
original, std_err, ci_bounds = bootstrap(native_percent_exposed_data.flatten(),  num_rounds=100000, func=np.mean, ci=0.95, seed=np.random.randint(100))
native_percent_exposed_mean = original
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

for vidx, variant in enumerate(['fast', 'slow', 'wt']):
    #get master equation stead state probability
    METS_data = np.loadtxt(f'msm_{pdb}/{variant}_METS.dat')[:,1:]
    print(METS_data, METS_data.shape)

    #get change data
    #msm/post_translation/fast/1-200/analysis/qbb/G_48_i33_post_trans_1yta.dat
    change_files = glob.glob(f'msm_{pdb}/post_translation/{variant}/1-200/analysis/qbb/G_*')
    change_dictionary = {}
    for change_file in change_files:

        traj = int(change_file.split('/')[-1].split('_')[1])
        print(change_file, traj)

        data = np.loadtxt(change_file, skiprows=8, usecols=(6))

        change_frames = np.where(data != 0)[0]
        print(f'Number of change frames: {len(change_frames)}')
        change_dictionary[traj] = change_frames


    #get SASA files for variant
    variant_files = sorted(glob.glob(f'{sys.argv[2]}/synthesis/*{variant}*{pdb}_ca_synth_charmm_sasa.txt*'), key=get_key)
    #print(variant_files, len(variant_files))
    variant_data = []
    for x in variant_files:
        print(x)
        variant_data.append(np.loadtxt(x, usecols=mask))
    variant_data = np.dstack(variant_data)
    variant_data = np.rollaxis(variant_data,-1)
    if len(mask) == 1:
        variant_data = variant_data.sum(axis=1)
    else:
        variant_data = variant_data.sum(axis=2)
    print(variant_data, variant_data.shape)
    variant_fig, variant_axs = plt.subplots(len(num_unique_msm_states)+1, sharex=True)
    for msm_state in num_unique_msm_states:
        print(f'\nAnalyzsing {variant} meta stable state: {msm_state}\n')
        idxs = np.where(msm_data[200*vidx:200*vidx+200,:] == msm_state)

        print(f'idxs: {idxs}', len(idxs))
        variant_state_data = variant_data[idxs]
        print(variant_state_data.shape)

        raw_variant_msm_sasa_data = variant_state_data
        plot_hist(raw_variant_msm_sasa_data.flatten(), msm_state, '{variant}', variant_axs)
        original, std_err, ci_bounds = bootstrap(raw_variant_msm_sasa_data.flatten(),  num_rounds=100000, func=np.mean, ci=0.95, seed=np.random.randint(100))
        var = (np.sqrt(raw_variant_msm_sasa_data.size)*std_err)**2
        #print( original, std_err, ci_bounds )
        lower_delta = original - ci_bounds[0]
        upper_delta = ci_bounds[1] - original
        raw_msm_outdata.append([ msm_state ,original, std_err, var, ci_bounds[0], ci_bounds[1], lower_delta, upper_delta])

        variant_msm_sasa_data = ((variant_state_data / native_original) - 1)*100
        original, std_err, ci_bounds = bootstrap(variant_msm_sasa_data.flatten(),  num_rounds=100000, func=np.mean, ci=0.95, seed=np.random.randint(100))
        var = (np.sqrt(variant_msm_sasa_data.size)*std_err)**2
        lower_delta = original - ci_bounds[0]
        upper_delta = ci_bounds[1] - original
        #print( original, std_err, ci_bounds )
        msm_outdata.append([ msm_state ,original, std_err, var, ci_bounds[0], ci_bounds[1], lower_delta, upper_delta])

        variant_msm_percent_exposed_data = ((variant_state_data / disordered_original))*100
        original, std_err, ci_bounds = bootstrap(variant_msm_percent_exposed_data.flatten(),  num_rounds=100000, func=np.mean, ci=0.95, seed=np.random.randint(100))
        var = (np.sqrt(variant_msm_sasa_data.size)*std_err)**2
        lower_delta = original - ci_bounds[0]
        upper_delta = ci_bounds[1] - original
        #print( original, std_err, ci_bounds )
        percent_exposed_outdata.append([ msm_state ,original, std_err, var, ci_bounds[0], ci_bounds[1], lower_delta, upper_delta])

        # contingency table analysis for the presence of a increase in the percent exposed and change in entanglement
        for tail in ['less','greater']:
            ctable_data = ctable_anal(change_dictionary, variant_msm_percent_exposed_data, idxs, native_percent_exposed_mean, tail)

            outfile = f'{out_path}{script_title}_{script_version}/output/{variant}_state{msm_state}_{mask_str}_ctable_anal_{tail}.out'
            np.savetxt(outfile, ctable_data, header=f'cutoff_{variant}_state{msm_state}_{mask_str}{tail}, oddsr_{variant}_state{msm_state}_{mask_str}{tail}, fisherp_{variant}_state{msm_state}_{mask_str}{tail}, log_fisherp_{variant}_state{msm_state}_{mask_str}{tail}, chi2_{variant}_state{msm_state}_{mask_str}{tail}, dof_{variant}_state{msm_state}_{mask_str}{tail}, chi2p_{variant}_state{msm_state}_{mask_str}{tail}, log_chi2p_{variant}_state{msm_state}_{mask_str}')
            print(f'SAVED: {outfile}')

        #Master equeation state_prob
        original, std_err, ci_bounds = bootstrap(METS_data[500:,msm_state].flatten(),  num_rounds=100000, func=np.mean, ci=0.95, seed=np.random.randint(100))
        var = (np.sqrt(METS_data[500:,msm_state].size)*std_err)**2
        lower_delta = original - ci_bounds[0]
        upper_delta = ci_bounds[1] - original
        METS_outdata.append([msm_state ,original, std_err, var, ci_bounds[0], ci_bounds[1], lower_delta, upper_delta])


    plt.suptitle = f'{variant}_SASA_histograms'
    plot_hist(native_data.flatten(), -1, 'Native' ,variant_axs, xlabel='SASA')
    variant_fig.tight_layout()
    variant_fig.savefig(f'{out_path}{script_title}_{script_version}/output/_{variant}_SASA_histograms_resids_{mask_str}.png')

    plt.clf()

    #data output
    msm_outdata = msm_outdata_base + msm_outdata
    raw_msm_outdata = raw_msm_outdata_base + raw_msm_outdata
    percent_exposed_outdata = percent_exposed_outdata + percent_exposed_outdata_base

    np.savetxt(f'{out_path}{script_title}_{script_version}/output/raw_{pdb}_{variant}_msm_{sys.argv[4]}_SASA.txt',raw_msm_outdata, header=f'state_raw_{variant}_{pdb}_{sys.argv[4]}, mean_raw_{variant}_{pdb}_{sys.argv[4]}, stderr_raw_{variant}_{pdb}_{sys.argv[4]}, var_raw_{variant}_{pdb}_{sys.argv[4]}, ci_lower_raw_{variant}_{pdb}_{sys.argv[4]}, ci_upper_raw_{variant}_{pdb}_{sys.argv[4]}, lower_delta_raw_{variant}_{pdb}_{sys.argv[4]}, upper_delta_raw_{variant}_{pdb}_{sys.argv[4]}', fmt='%s,%s,%s,%s,%s,%s,%s,%s')
    print(f'Saved: {out_path}{script_title}_{script_version}/output/raw_{pdb}_{variant}_msm_{sys.argv[4]}_SASA.txt')

    np.savetxt(f'{out_path}{script_title}_{script_version}/output/zeta_{pdb}_{variant}_msm_{sys.argv[4]}_SASA.txt',msm_outdata, header=f'state_zeta_{variant}_{pdb}_{sys.argv[4]}, mean_zeta_{variant}_{pdb}_{sys.argv[4]}, stderr_zeta_{variant}_{pdb}_{sys.argv[4]}, var_zeta_{variant}_{pdb}_{sys.argv[4]}, ci_lower_zeta_{variant}_{pdb}_{sys.argv[4]}, ci_upper_zeta_{variant}_{pdb}_{sys.argv[4]}, lower_delta_zeta_{variant}_{pdb}_{sys.argv[4]}, upper_delta_zeta_{variant}_{pdb}_{sys.argv[4]},  {variant}_zeta_{variant}', fmt='%s,%s,%s,%s,%s,%s,%s,%s')
    print(f'Saved: {out_path}{script_title}_{script_version}/output/zeta_{pdb}_{variant}_msm_{sys.argv[4]}_SASA.txt')

    np.savetxt(f'{out_path}{script_title}_{script_version}/output/percent_exposed_{pdb}_{variant}_msm_{sys.argv[4]}_SASA.txt',percent_exposed_outdata, header=f'state_perce_{variant}_{pdb}_{sys.argv[4]}, mean_perce_{variant}_{pdb}_{sys.argv[4]}, stderr_perce_{variant}_{pdb}_{sys.argv[4]}, var_perce_{variant}_{pdb}_{sys.argv[4]}, ci_lower_perce_{variant}_{pdb}_{sys.argv[4]}, ci_upper_perce_{variant}_{pdb}_{sys.argv[4]}, lower_delta_perce_{variant}_{pdb}_{sys.argv[4]}, upper_delta_perce_{variant}_{pdb}_{sys.argv[4]}', fmt='%s,%s,%s,%s,%s,%s,%s,%s')
    print(f'Saved: {out_path}{script_title}_{script_version}/output/percent_exposed_{pdb}_{variant}_msm_{sys.argv[4]}_SASA.txt')

    np.savetxt(f'{out_path}{script_title}_{script_version}/output/METS_{pdb}_{variant}_msm_{sys.argv[4]}_SASA.txt',METS_outdata, header=f'state_prob_{variant}_{pdb}_{sys.argv[4]}, mean_prob_{variant}_{pdb}_{sys.argv[4]}, stderr_prob_{variant}_{pdb}_{sys.argv[4]}, var_prob_{variant}_{pdb}_{sys.argv[4]}, ci_lower_prob_{variant}_{pdb}_{sys.argv[4]}, ci_upper_prob_{variant}_{pdb}_{sys.argv[4]}, lower_delta_prob_{variant}_{pdb}_{sys.argv[4]}, upper_delta_prob_{variant}_{pdb}_{sys.argv[4]}', fmt='%s,%s,%s,%s,%s,%s,%s,%s')
    print(f'Saved: {out_path}{script_title}_{script_version}/output/METS_{pdb}_{variant}_msm_{sys.argv[4]}_SASA.txt')

