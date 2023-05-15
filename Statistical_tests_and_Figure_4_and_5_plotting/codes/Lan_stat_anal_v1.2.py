#!/usr/bin/env python3
import sys,os
import pickle
import numpy as np
import glob
import pandas as pd
from scipy.stats import bootstrap
import itertools
from scipy.stats import permutation_test, fisher_exact
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec

if len(sys.argv) != 3:
    print('[1] path to directory containing processed csvs')
    print('[2] outpath')
    quit()

# Permanently changes the pandas settings
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', -1)

def statistic(x, y, axis):
    return np.mean(x, axis=axis) - np.mean(y, axis=axis)


files = glob.glob(sys.argv[1])
print(files)

for f_i, f in enumerate(files):
    print(f_i, f)
    if f_i == 0:
        df = pd.read_csv(f)
    else:
        df = df.append(pd.read_csv(f), ignore_index=True)
print(df.shape)

pdbs = np.unique(df['PDB'])
variants = np.unique(df['variant'])
subunits = np.unique(df['subunit'])

stat_dict = {}
stat_dict['PDB'] = []
stat_dict['variant'] = []
stat_dict['subunit'] = []
stat_dict['frac_ent_pres'] = []
stat_dict['frac_ent_pres_cil'] = []
stat_dict['frac_ent_pres_ciu'] = []
stat_dict['frac_ent_pres@int'] = []
stat_dict['frac_ent_pres@int_cil'] = []
stat_dict['frac_ent_pres@int_ciu'] = []

permutation_data_dict = {}
ent_present_ctable_dict = {}
ent_present_at_int_ctable_dict = {}
#            EntY   EntN
# weakB Y    [[],    [],
# weak  N    [],    []]
#is there an association between weak binding and Ent being present at the interface
ctable_output_dict = {}
for pdb in pdbs:
    for variant in variants:
        for subunit in subunits:
            print(pdb, variant, subunit)

            stat_dict['PDB'] += [pdb]
            stat_dict['variant'] += [variant]
            stat_dict['subunit'] += [subunit]


            if pdb not in ent_present_ctable_dict:
                ent_present_ctable_dict[pdb] = {}
            if variant not in ent_present_ctable_dict[pdb]:
                ent_present_ctable_dict[pdb][variant] = {}

            if pdb not in ent_present_at_int_ctable_dict:
                ent_present_at_int_ctable_dict[pdb] = {}
            if variant not in ent_present_at_int_ctable_dict[pdb]:
                ent_present_at_int_ctable_dict[pdb][variant] = {}

            loc_df = df[df['PDB'] == pdb]
            loc_df = loc_df[loc_df['variant'] == variant]
            loc_df = loc_df[loc_df['subunit'] == subunit]

            IDs = np.unique(loc_df['ID'])

            if subunit == 'dimer':
                print('DIMER')
                ent_present = []
                ent_present_at_int = []
                for ID in IDs:
                    loc_ent_df = loc_df[loc_df['ID'] == ID]
                    ent_pres = any(loc_ent_df['ent_present'])
                    ent_presint = any(loc_ent_df['ent_present@int_method2'])

                    if ent_pres:
                        ent_present += [1]
                        ent_present_flag = 1
                    else:
                        ent_present += [0]
                        ent_present_flag = 0

                    if ent_presint:
                        ent_present_at_int += [1]
                        ent_present_at_int_flag = 1
                    else:
                        ent_present_at_int += [0]
                        ent_present_at_int_flag = 0

                    #ctable
                    for cutoff in np.arange(5, 100, 5):

                        if cutoff not in ent_present_ctable_dict[pdb][variant]:
                            ent_present_ctable_dict[pdb][variant][cutoff] = {'ctable':np.asarray([[0,0], [0,0]]), 'pvalue':[], 'OR':[]}
                        if cutoff not in ent_present_at_int_ctable_dict[pdb][variant]:
                            ent_present_at_int_ctable_dict[pdb][variant][cutoff] = {'ctable':np.asarray([[0,0], [0,0]]), 'pvalue':[], 'OR':[]}

                        weak_binding = any(loc_ent_df[f'weak_binder@{cutoff}%'])

                        if ent_present_flag == 1 and weak_binding == True:
                            ent_present_ctable_dict[pdb][variant][cutoff]['ctable'] += np.asarray([[1,0], [0,0]])
                        if ent_present_flag == 0 and weak_binding == True:
                            ent_present_ctable_dict[pdb][variant][cutoff]['ctable'] += np.asarray([[0,1], [0,0]])
                        if ent_present_flag == 1 and weak_binding == False:
                            ent_present_ctable_dict[pdb][variant][cutoff]['ctable'] += np.asarray([[0,0], [1,0]])
                        if ent_present_flag == 0 and weak_binding == False:
                            ent_present_ctable_dict[pdb][variant][cutoff]['ctable'] += np.asarray([[0,0], [0,1]])

                        if ent_present_at_int_flag == 1 and weak_binding == True:
                            ent_present_at_int_ctable_dict[pdb][variant][cutoff]['ctable'] += np.asarray([[1,0], [0,0]])
                        if ent_present_at_int_flag == 0 and weak_binding == True:
                            ent_present_at_int_ctable_dict[pdb][variant][cutoff]['ctable'] += np.asarray([[0,1], [0,0]])
                        if ent_present_at_int_flag == 1 and weak_binding == False:
                            ent_present_at_int_ctable_dict[pdb][variant][cutoff]['ctable'] += np.asarray([[0,0], [1,0]])
                        if ent_present_at_int_flag == 0 and weak_binding == False:
                            ent_present_at_int_ctable_dict[pdb][variant][cutoff]['ctable'] += np.asarray([[0,0], [0,1]])



                for cutoff in np.arange(5, 100, 5):

                    if cutoff not in ctable_output_dict:
                        ctable_output_dict[cutoff] = {}
                    if pdb not in ctable_output_dict[cutoff]:
                        ctable_output_dict[cutoff][pdb] = {}
                    if variant not in ctable_output_dict[cutoff][pdb]:
                        ctable_output_dict[cutoff][pdb][variant] = {}
                    if 'entanglement_present' not in ctable_output_dict[cutoff][pdb][variant]:
                        ctable_output_dict[cutoff][pdb][variant]['entanglement_present'] = ''
                    if 'entanglement_present@int' not in ctable_output_dict[cutoff][pdb][variant]:
                        ctable_output_dict[cutoff][pdb][variant]['entanglement_present@int'] = ''


                    #for ent present at all
                    ctable = ent_present_ctable_dict[pdb][variant][cutoff]['ctable']
                    OR, pvalue = fisher_exact(ctable, alternative='two-sided')

                    ctable_str = ""
                    ctable_str += f"{'Entangled?': >40}"
                    ctable_str += f"\n{'Y    N': >37}"
                    ctable_str += f"\nWeak Binding Y  {ctable[0,0]: <3}  {ctable[0,1]}"
                    ctable_str += f"\nWeak Binding N  {ctable[1,0]: <3}  {ctable[1,1]}"
                    ctable_str += f"\nOR  {OR:.3E}"
                    ctable_str += f"\npvalue  {pvalue:.3E}"
                    ctable_output_dict[cutoff][pdb][variant]['entanglement_present'] = ctable_str

                    ent_present_ctable_dict[pdb][variant][cutoff]['pvalue'] += [pvalue]
                    ent_present_ctable_dict[pdb][variant][cutoff]['OR'] += [OR]

                    #for ent present at interface
                    ctable = ent_present_at_int_ctable_dict[pdb][variant][cutoff]['ctable']
                    OR, pvalue = fisher_exact(ctable, alternative='two-sided')

                    ctable_str = ""
                    ctable_str += f"{'Entangled@int?': >40}"
                    ctable_str += f"\n{'Y    N': >37}"
                    ctable_str += f"\nWeak Binding Y  {ctable[0,0]: <3}  {ctable[0,1]}"
                    ctable_str += f"\nWeak Binding N  {ctable[1,0]: <3}  {ctable[1,1]}"
                    ctable_str += f"\nOR  {OR:.3E}"
                    ctable_str += f"\npvalue  {pvalue:.3E}"
                    ctable_output_dict[cutoff][pdb][variant]['entanglement_present@int'] = ctable_str

                    ent_present_at_int_ctable_dict[pdb][variant][cutoff]['pvalue'] += [pvalue]
                    ent_present_at_int_ctable_dict[pdb][variant][cutoff]['OR'] += [OR]

                #get bootstrapped confidence interfvals
                rng = np.random.default_rng()
                data = (ent_present, )
                res = bootstrap(data, np.mean, confidence_level=0.95, random_state=rng, n_resamples=100000)
                ci_l, ci_u = res.confidence_interval
                mean = np.mean(ent_present)

                stat_dict['frac_ent_pres'] += [mean]
                stat_dict['frac_ent_pres_cil'] += [ci_l]
                stat_dict['frac_ent_pres_ciu'] += [ci_u]

                rng = np.random.default_rng()
                data = (ent_present_at_int,)
                res = bootstrap(data, np.mean, confidence_level=0.95, n_resamples=100000)
                ci_l, ci_u = res.confidence_interval
                mean = np.mean(ent_present_at_int)

                stat_dict['frac_ent_pres@int'] += [mean]
                stat_dict['frac_ent_pres@int_cil'] += [ci_l]
                stat_dict['frac_ent_pres@int_ciu'] += [ci_u]


            if subunit == 'mono':
                print('MONO')
                ent_present = loc_df['ent_present']
                ent_present_at_int = loc_df['ent_present@int_method2']

                rng = np.random.default_rng()
                data = (ent_present, )
                res = bootstrap(data, np.mean, confidence_level=0.95, random_state=rng, n_resamples=100000)
                ci_l, ci_u = res.confidence_interval
                mean = np.mean(ent_present)
                print(ci_l, mean, ci_u)

                stat_dict['frac_ent_pres'] += [mean]
                stat_dict['frac_ent_pres_cil'] += [ci_l]
                stat_dict['frac_ent_pres_ciu'] += [ci_u]

                rng = np.random.default_rng()
                data = (ent_present_at_int,)
                res = bootstrap(data, np.mean, confidence_level=0.95, random_state=rng, n_resamples=100000)
                ci_l, ci_u = res.confidence_interval
                mean = np.mean(ent_present_at_int)
                print(ci_l, mean, ci_u)

                stat_dict['frac_ent_pres@int'] += [mean]
                stat_dict['frac_ent_pres@int_cil'] += [ci_l]
                stat_dict['frac_ent_pres@int_ciu'] += [ci_u]

            #save data for permutateion
            key = (pdb, variant, subunit)
            if key not in permutation_data_dict:
                permutation_data_dict[key] = [ent_present, ent_present_at_int]


stat_df = pd.DataFrame(stat_dict)
stat_df = stat_df.sort_values(['PDB', 'subunit'])
print(stat_df)
stat_df.to_csv('state_fraction_df.csv', index=False)
print(f'SAVED: state_fraction_df.csv')

#output ctable .csv
ctable_df = {'PDB': [], 'variant': [], 'cutoff': [], 'entanglement': [], 'ctable': []}
writer = pd.ExcelWriter('Contingency_tables.xlsx', engine='xlsxwriter')
workbook  = writer.book
wrap_format = workbook.add_format({'text_wrap': True})
for pdb in pdbs:
    for ent_loc in ['entanglement_present', 'entanglement_present@int']:
        sheet_df = {'cutoff': []}
        for variant in variants:
            if variant not in sheet_df:
                sheet_df[f'{variant}'] = []

            for cutoff in np.arange(5, 100, 5):

                if cutoff not in sheet_df['cutoff']:
                    sheet_df['cutoff'] += [cutoff]
                sheet_df[f'{variant}'] += [ctable_output_dict[cutoff][pdb][variant][ent_loc]]


                ctable_df['PDB'] += [pdb]
                ctable_df['variant'] += [variant]
                ctable_df['cutoff'] += [cutoff]
                ctable_df['entanglement'] += [ent_loc]
                ctable_df['ctable'] += [ctable_output_dict[cutoff][pdb][variant][ent_loc]]

        sheet_df = pd.DataFrame(sheet_df)

        if ent_loc == 'entanglement_present':
            sheet_df.to_excel(writer, sheet_name=f'{pdb}_ent', index=False)
            writer.sheets[f'{pdb}_ent'].set_column(5, 0, 25, wrap_format)


        if ent_loc == 'entanglement_present@int':
            sheet_df.to_excel(writer, sheet_name=f'{pdb}_ent@int', index=False)
            writer.sheets[f'{pdb}_ent@int'].set_column(5, 0, 25, wrap_format)


writer.save()
print(f'SAVED: Contingency_tables.xlsx')

ctable_df = pd.DataFrame(ctable_df)
ctable_df.to_csv('Contingency_tables.csv')
print(f'SAVED: Contingency_tables.csv')

### permutation test for pvalues in difference between state fractions of entanglements at the interface or not.
permutation_pvalue_dict = {}
for pdb in pdbs:
    for subunit in subunits:

        permut_keys = [k for k in permutation_data_dict.keys() if pdb in k and subunit in k]

        for pair in itertools.combinations(permut_keys, 2):

            ent_present1, ent_present_at_int1 = permutation_data_dict[pair[0]]
            ent_present2, ent_present_at_int2 = permutation_data_dict[pair[1]]

            data = (ent_present1, ent_present2)
            rng = np.random.default_rng()

            mean1 = np.mean(ent_present1)
            mean2 = np.mean(ent_present2)

            if mean1 > mean2:
                res1 = permutation_test(data, statistic, n_resamples=100000, vectorized=True, alternative='greater', random_state=rng)
            if mean1 < mean2:
                res1 = permutation_test(data, statistic, n_resamples=100000, vectorized=True, alternative='less', random_state=rng)

            data = (ent_present_at_int1, ent_present_at_int2)
            rng = np.random.default_rng()

            mean1 = np.mean(ent_present_at_int1)
            mean2 = np.mean(ent_present_at_int2)

            if mean1 > mean2:
                res2 = permutation_test(data, statistic, n_resamples=100000, vectorized=True, alternative='greater', random_state=rng)
            if mean1 < mean2:
                res2 = permutation_test(data, statistic, n_resamples=100000, vectorized=True, alternative='less', random_state=rng)

            permutation_pvalue_dict[pair] = [res1.pvalue, res2.pvalue]

print('pair1, pair2, pvalue for ent in general, pvalue for ent at int')
for pair,res in permutation_pvalue_dict.items():

    res1_str = 'ns'
    res2_str  = 'ns'
    if res[0] < 0.05:
        res1_str = '*'
    if res[0] < 0.01:
        res1_str = '**'
    if res[0] < 0.001:
        res1_str = '***'
    if res[1] < 0.05:
        res2_str = '*'
    if res[1] < 0.01:
        res2_str = '**'
    if res[1] < 0.001:
        res2_str = '***'

    print(pair, res[0], res1_str, res[1], res2_str)


#plot fraction ent figure
print('making fraction ent figure')
fig, axs = plt.subplots(nrows=2, ncols=2, sharex=False, sharey=False)
order_variants_xaxis = ['wt', 'fast', 'slow', 'wt', 'fast', 'slow']
order_pdbs_xaxis = ['1yta', '1yta', '1yta', '2is3', '2is3', '2is3']
xlabel_pos = [0, 5, 10, 16, 21, 26]
panel_labels = np.asarray([['a', 'b'], ['c', 'd']])
print(panel_labels)
colors = ['#DDAA32', '#BB5566', '#004488', '#DDAA32', '#BB5566', '#004488']
for subunit_idx, subunit in enumerate(['mono', 'dimer']):
    print(subunit_idx, subunit)
    loc_df = stat_df[stat_df['subunit']==subunit]
    print(loc_df)

    for ent_type_idx, ent_type in enumerate(['Overall gain in entanglement', 'Gain in entanglement at interface']):
        print(ent_type_idx, ent_type)
        ax = axs[subunit_idx, ent_type_idx]

        if ent_type_idx == 0:
            for v,p,i,c in zip(order_variants_xaxis, order_pdbs_xaxis, xlabel_pos, colors):
                print(v,p,i)
                dp_data = loc_df[loc_df['PDB'] == p]
                dp_data = dp_data[dp_data['variant'] == v]
                print(dp_data)
                cil_delta = dp_data['frac_ent_pres'] - dp_data['frac_ent_pres_cil']
                ciu_delta = dp_data['frac_ent_pres_ciu'] - dp_data['frac_ent_pres']
                ax.errorbar(i, dp_data['frac_ent_pres'], yerr=[cil_delta, ciu_delta], ls='None', marker='o', color=c, ecolor='k', capsize=2, markersize=4)

            ax.set_xticks(ticks=xlabel_pos, labels=[x.upper() for x in order_variants_xaxis], fontsize=6)
            ax.tick_params(axis='x', labelsize=6)
            ax.tick_params(axis='y', labelsize=6)
            ax.set_ylabel('State fraction\nentanglements present', fontsize=8)
            ax.set_xlabel('Oligoribonuclease   Ribonuclease T', fontsize=8)
            ax.set_xlim([-1, 27])
            ax.set_ylim([-0.05,1])
            ax.get_yaxis().set_label_coords(-0.2,0.5)
            ax.text(-0.2,1.15, panel_labels[subunit_idx, ent_type_idx], transform=ax.transAxes, fontsize=10, va='top', ha='right')

            if subunit == 'mono':
                ax.text(0.5, 1.1, 'Monomers resulting from synthesis', ha='center', va='center', transform=ax.transAxes, fontsize=7)
            if subunit == 'dimer':
                ax.text(0.5, 1.1, 'Annealed dimers', ha='center', va='center', transform=ax.transAxes, fontsize=7)

        if ent_type_idx == 1:
            for v,p,i,c in zip(order_variants_xaxis, order_pdbs_xaxis, xlabel_pos, colors):
                print(v,p,i)
                dp_data = loc_df[loc_df['PDB'] == p]
                dp_data = dp_data[dp_data['variant'] == v]
                print(dp_data)
                cil_delta = dp_data['frac_ent_pres@int'] - dp_data['frac_ent_pres@int_cil']
                ciu_delta = dp_data['frac_ent_pres@int_ciu'] - dp_data['frac_ent_pres@int']
                print(cil_delta, ciu_delta)
                ax.errorbar(i, dp_data['frac_ent_pres@int'], yerr=[cil_delta, ciu_delta], ls='None', marker='o', color=c, ecolor='k', capsize=2, markersize=4)

            ax.set_xticks(ticks=xlabel_pos, labels=[x.upper() for x in order_variants_xaxis], fontsize=6)
            ax.tick_params(axis='x', labelsize=6)
            ax.tick_params(axis='y', labelsize=6)
            ax.set_ylabel('State fraction\nentanglements at interface', fontsize=8)
            ax.set_xlabel('Oligoribonuclease   Ribonuclease T', fontsize=8)
            ax.set_xlim([-1, 27])
            ax.set_ylim([-0.05,1])
            ax.get_yaxis().set_label_coords(-0.2,0.5)
            ax.text(-0.2,1.15, panel_labels[subunit_idx, ent_type_idx], transform=ax.transAxes, fontsize=10, va='top', ha='right')

            if subunit == 'mono':
                ax.text(0.5, 1.1, 'Monomers resulting from synthesis', ha='center', va='center', transform=ax.transAxes, fontsize=7)
            if subunit == 'dimer':
                ax.text(0.5, 1.1, 'Annealed dimers', ha='center', va='center', transform=ax.transAxes, fontsize=7)

plt.subplots_adjust(wspace = 0.5)
plt.subplots_adjust(hspace = 0.5)
plt.savefig('state_fraction_fig.png', dpi=600)
print('SAVED: state_fraction_fig.png')
plt.clf()

#plot contingency figure with OR and Pvalues from Fisher  Exact test
colors = {'wt':'#DDAA32', 'fast':'#BB5566', 'slow':'#004488'}
#fig, axs = plt.subplots(nrows=5, ncols=2, sharex=False, sharey=False, figsize=(6,6))
fig = plt.figure(figsize=(6,6))
gs = fig.add_gridspec(3, 1, hspace=0.65)
gs0 = gs[0].subgridspec(2, 2, hspace=0.3, wspace=0.4)
gs1 = gs[1].subgridspec(2, 2, hspace=0.3, wspace=0.4)
gs2 = gs[2].subgridspec(1, 2, wspace=0.4)
set1, set1_sub = gs0.subplots()
set2, set2_sub = gs1.subplots()
set3 = gs2.subplots()
panel_labels = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j']
panel_idx = 0

for edict_idx, edict in {0:ent_present_ctable_dict, 2:ent_present_at_int_ctable_dict}.items():
    print(edict_idx)

    for pdb_idx, pdb in enumerate(edict.keys()):
        print(pdb_idx, pdb)

        if edict_idx == 0:
            pvalue_axis = set1[pdb_idx]
        else:
            pvalue_axis = set2[pdb_idx]

        pvalue_axis.set_xlim([0, 100])

        if panel_idx in [0]:
            pvalue_axis.set_ylim([0, 30])
            pvalue_axis.set_yticks(np.arange(0, 31, 10))
        elif panel_idx in [1]:
            pvalue_axis.set_ylim([0, 3])
            pvalue_axis.set_yticks(np.arange(0, 3.1, 0.5))
        elif panel_idx in [2]:
            pvalue_axis.set_ylim([0, 20])
            pvalue_axis.set_yticks(np.arange(0, 21, 5))
        elif panel_idx in [3]:
            pvalue_axis.set_ylim([0, 3])
            pvalue_axis.set_yticks(np.arange(0, 3.1, 0.5))

        pvalue_axis.set_ylabel(r'$-\mathrm{Log_{10}}\mathrm{(}p_{value}\mathrm{)}$', fontsize=6, va='center', ha='center')
        pvalue_axis.axhline(y=1.3, color='r', ls='--', linewidth=1)

        if edict_idx == 0:
            if pdb =='1yta':
                pvalue_axis.set_title(f'Oligoribonuclease', fontsize=10, pad=15)
            if pdb =='2is3':
                pvalue_axis.set_title(f'Ribonuclease T', fontsize=10, pad=15)

        if edict_idx == 0:
            OR_axis = set1_sub[pdb_idx]
        else:
            OR_axis = set2_sub[pdb_idx]

        OR_axis.set_xlim([0, 100])
        OR_axis.set_ylabel(f'OR', fontsize=6, va='center', ha='center')
        OR_axis.axhline(y=1, color='r', ls='--', linewidth=1)

        if panel_idx in [0]:
            OR_axis.set_ylim([0.5,110])
        elif panel_idx in [1]:
            OR_axis.set_ylim([0.5,10])
        elif panel_idx in [2]:
            OR_axis.set_ylim([0.5,50])
        elif panel_idx in [3]:
            OR_axis.set_ylim([0.5,10])

        OR_axis.set_xlabel('Weak binding region lower bound,\nPercentile of WT CDF', fontsize=6)
        OR_axis.set_yscale("log")

        #intE_axis = axs[4, pdb_idx]
        intE_axis = set3[pdb_idx]
        intE_axis.set_xlabel('Interaction Energy, kcal/mol', fontsize=6)
        intE_axis.set_ylabel('CDF', fontsize=6, va='center', ha='center')

        #align yaxis labels
        pvalue_axis.get_yaxis().set_label_coords(-0.2,0.5)
        OR_axis.get_yaxis().set_label_coords(-0.2,0.5)
        intE_axis.get_yaxis().set_label_coords(-0.2,0.5)
        pvalue_axis.text(-0.2,1.3, panel_labels[panel_idx], transform=pvalue_axis.transAxes, fontsize=10, va='center', ha='center')
        panel_idx += 1

        #get those cutoffs where no variant has an invalid OR
        check_list = []
        for variant in edict[pdb].keys():
            cutoff_list = []

            var_check_list = []
            for cutoff in edict[pdb][variant].keys():
                pvalue = (-1)*np.log10(edict[pdb][variant][cutoff]['pvalue'])
                OR = edict[pdb][variant][cutoff]['OR']

                cutoff_list += [cutoff]
                if OR[0] == np.inf:
                    var_check_list += [0]
                else:
                    var_check_list += [1]

            check_list += [var_check_list]
        sum_passes = np.sum(np.vstack(check_list).T, axis=1)
        check_list = np.vstack((cutoff_list, check_list, sum_passes)).T


        for variant in edict[pdb].keys():

            for cutoff in edict[pdb][variant].keys():
                pvalue = (-1)*np.log10(edict[pdb][variant][cutoff]['pvalue'])
                OR = edict[pdb][variant][cutoff]['OR']

                sum_pass_for_cutoff = check_list[np.where(check_list[:,0]==cutoff)][0][-1]

                if sum_pass_for_cutoff == 3:

                    c = colors[variant]
                    pvalue_axis.plot(cutoff, pvalue, c=c, ls='None', marker='o', markersize=2)
                    OR_axis.plot(cutoff, OR, c=c, ls='None', marker='o', markersize=2)


            Eint_file = f'interaction_energy_files/{pdb}_{variant}_inte.dat'
            Eint_data = sorted(np.loadtxt(Eint_file))

            N = len(Eint_data)
            cum = []
            for d_i, d in enumerate(Eint_data):
                cum += [[d, (d_i+1)/N]]
            cum = np.vstack(cum)
            intE_axis.plot(cum[:,0], cum[:,1], c=colors[variant])


        labels = [item.get_text() for item in pvalue_axis.get_xticklabels()]

        empty_string_labels = ['']*len(labels)
        pvalue_axis.set_xticklabels(empty_string_labels)
        pvalue_axis.tick_params(axis='both', which='minor', labelsize=6)
        OR_axis.tick_params(axis='both', which='minor', labelsize=6)
        intE_axis.tick_params(axis='both', which='minor', labelsize=6)

set3[0].text(-0.2,1.1, panel_labels[panel_idx], transform=set3[0].transAxes, fontsize=10, va='center', ha='center')
set3[1].text(-0.2,1.1, panel_labels[panel_idx+1], transform=set3[1].transAxes, fontsize=10, va='center', ha='center')

set1[0].tick_params(axis='x', labelsize=6)
set1[0].tick_params(axis='y', labelsize=6)
set1_sub[0].tick_params(axis='x', labelsize=6)
set1_sub[0].tick_params(axis='y', labelsize=6)
set2[0].tick_params(axis='x', labelsize=6)
set2[0].tick_params(axis='y', labelsize=6)
set2_sub[0].tick_params(axis='x', labelsize=6)
set2_sub[0].tick_params(axis='y', labelsize=6)
set3[0].tick_params(axis='x', labelsize=6)
set3[0].tick_params(axis='y', labelsize=6)
set1[1].tick_params(axis='x', labelsize=6)
set1[1].tick_params(axis='y', labelsize=6)
set1_sub[1].tick_params(axis='x', labelsize=6)
set1_sub[1].tick_params(axis='y', labelsize=6)
set2[1].tick_params(axis='x', labelsize=6)
set2[1].tick_params(axis='y', labelsize=6)
set2_sub[1].tick_params(axis='x', labelsize=6)
set2_sub[1].tick_params(axis='y', labelsize=6)
set3[1].tick_params(axis='x', labelsize=6)
set3[1].tick_params(axis='y', labelsize=6)
plt.savefig('contingency_anal.svg', dpi=600)
print(f'SAVED: contingency_anal.svg')
