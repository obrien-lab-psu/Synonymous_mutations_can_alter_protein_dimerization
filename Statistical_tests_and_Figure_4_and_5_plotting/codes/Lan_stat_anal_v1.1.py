#!/usr/bin/env python3
import sys,os
import pickle
import numpy as np
import glob
import pandas as pd
from scipy.stats import bootstrap

if len(sys.argv) != 3:
    print('[1] path to directory containing processed csvs')
    print('[2] outpath')
    quit()

# Permanently changes the pandas settings
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', -1)


#processes data and make a master output dataframe
#master_dict = {}
#master_dict['PDB'] = []
#master_dict['variant'] = []
#master_dict['subunit'] = []
#master_dict['ID'] = []
#master_dict['label'] = []
#master_dict['ent_present'] = []
#master_dict['ent_present@int_method1'] = []
#master_dict['ent_present@int_method2'] = []
#master_dict['Eint'] = []
#for cutoff in np.arange(5, 100, 5):
#    print(cutoff)
#    master_dict[f'weak_binder@{cutoff}%'] = []
#print(master_dict)

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
print(pdbs)
print(variants)
print(subunits)

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

for pdb in pdbs:
    for variant in variants:
        for subunit in subunits:
            #print(pdb, variant, subunit)

            stat_dict['PDB'] += [pdb]
            stat_dict['variant'] += [variant]
            stat_dict['subunit'] += [subunit]


            loc_df = df[df['PDB'] == pdb]
            loc_df = loc_df[loc_df['variant'] == variant]
            loc_df = loc_df[loc_df['subunit'] == subunit]
            #print(loc_df)

            IDs = np.unique(loc_df['ID'])
            #print(IDs)

            if subunit == 'dimer':
                #print('DIMER')
                ent_present = []
                ent_present_at_int = []
                for ID in IDs:
                    loc_ent_df = loc_df[loc_df['ID'] == ID]
                    #print(loc_ent_df)

                    #print(loc_ent_df['ent_present'])
                    #print(loc_ent_df['ent_present@int_method2'])
                    ent_pres = any(loc_ent_df['ent_present'])
                    ent_presint = any(loc_ent_df['ent_present@int_method2'])
                    #print(ent_pres, ent_presint)

                    if ent_pres:
                        ent_present += [1]
                    else:
                        ent_present += [0]

                    if ent_presint:
                        ent_present_at_int += [1]
                    else:
                        ent_present_at_int += [0]

                #print(ent_present)
                #print(ent_present_at_int)

                rng = np.random.default_rng()
                data = (ent_present, )
                res = bootstrap(data, np.mean, confidence_level=0.95, random_state=rng)
                ci_l, ci_u = res.confidence_interval
                mean = np.mean(ent_present)
                #print(ci_l, mean, ci_u)

                stat_dict['frac_ent_pres'] += [mean]
                stat_dict['frac_ent_pres_cil'] += [ci_l]
                stat_dict['frac_ent_pres_ciu'] += [ci_u]

                rng = np.random.default_rng()
                data = (ent_present_at_int,)
                res = bootstrap(data, np.mean, confidence_level=0.95, random_state=rng)
                mean = np.mean(ent_present_at_int)
                #print(ci_l, mean, ci_u)

                stat_dict['frac_ent_pres@int'] += [mean]
                stat_dict['frac_ent_pres@int_cil'] += [ci_l]
                stat_dict['frac_ent_pres@int_ciu'] += [ci_u]


            if subunit == 'mono':
                #print('MONO')
                ent_present = loc_df['ent_present']
                ent_present_at_int = loc_df['ent_present@int_method2']

                rng = np.random.default_rng()
                data = (ent_present, )
                res = bootstrap(data, np.mean, confidence_level=0.95, random_state=rng)
                ci_l, ci_u = res.confidence_interval
                mean = np.mean(ent_present)
                #print(ci_l, mean, ci_u)

                stat_dict['frac_ent_pres'] += [mean]
                stat_dict['frac_ent_pres_cil'] += [ci_l]
                stat_dict['frac_ent_pres_ciu'] += [ci_u]

                rng = np.random.default_rng()
                data = (ent_present_at_int,)
                res = bootstrap(data, np.mean, confidence_level=0.95, random_state=rng)
                mean = np.mean(ent_present_at_int)
                #print(ci_l, mean, ci_u)

                stat_dict['frac_ent_pres@int'] += [mean]
                stat_dict['frac_ent_pres@int_cil'] += [ci_l]
                stat_dict['frac_ent_pres@int_ciu'] += [ci_u]


stat_df = pd.DataFrame(stat_dict)
stat_df = stat_df.sort_values(['PDB', 'subunit'])
print(stat_df)
quit()
#calculate
#create array of [pdb, variant, monomer or dimerID, ent, ent@interface, Eint(if dimer), weak binder @ various cutoffs binary flag]
#(1) fraction of dimers with entanglements
#(2) fraction of dimers with entanglements at the interface


#create array of [pdb, variant, monomerID, ent, ent@interface]
#(3) fraction of monomers with entanglements
#(4) fraction of monomers with entanglements at the interface

#contingency table analysis for varying Eint percentiles
