#!/usr/bin/env python3
import sys, getopt, math, os, multiprocessing, time, traceback
import numpy as np
import parmed as pmd
import mdtraj as mdt

################################# Arguments ###################################
# Default parameters
n_traj = 100
mutant_type_list = ['fast', 'slow']
act_mask = ''
prefix_dir = ''

# read control file
ctrlfile = ''

usage = 'python get_post_trans_order_parameters_v1.1.py [1] [2]\n[1] control file path\n[2] PDBID'
if len(sys.argv) != 3:
    print(usage)
    sys.exit()

ctrlfile = sys.argv[1]
pdb = sys.argv[2]

if not os.path.exists(ctrlfile):
    print('Error: cannot find control file ' + ctrlfile + '.')
    sys.exit()

file_object = open(ctrlfile,'r')
try:
    for line in file_object:
        line = line.strip()
        if not line:
            # This is a blank line
            continue
        if line.startswith('#'):
            # This is a comment line
            continue
        if line.startswith('n_traj'):
            words = line.split('=')
            n_traj = int(words[1].strip())
            continue
        if line.startswith('mutant_type_list'):
            words = line.split('=')
            mutant_type_list = words[1].strip().split()
            continue
        if line.startswith('act_mask'):
            words = line.split('=')
            act_mask = words[1].strip()
            continue
        if line.startswith('prefix_dir'):
            words = line.split('=')
            prefix_dir = words[1].strip()
            continue

finally:
     file_object.close()

################################# Functions ###################################
def get_Qbb_act():
    global n_traj, prefix_dir, act_mask, mutant_type, pdb
    print(mutant_type+' QBB_act:')
    qbb_list = [];
    for i in range(n_traj):
        Q_file = f'{prefix_dir}Q/{mutant_type}/qbb_{i+1}_i33_post_trans_{pdb}.dat'
        print(Q_file)
        f = open(f'{Q_file}')
        C = f.readlines()
        f.close()
        os.system('rm -rf tmp/')
        C = [C[k].strip().split() for k in range(1,len(C))]
        C = np.array(C, dtype=np.float32)
        Q_ts_0 = C[:,-1].reshape((len(C[:,-1]),1))
        qbb_list.append(Q_ts_0)

    np.save('%s_QBB_act.npy'%mutant_type, qbb_list)
    print(f'Saved: {mutant_type}_QBB_act.npy')

def get_entanglement():
    global n_traj, prefix_dir, act_mask, mutant_type, pdb
    print(mutant_type+' G:')
    G_number_list = [];
    for i in range(n_traj):
        G_file = f'{prefix_dir}G/{mutant_type}/G_{i+1}_i33_post_trans_{pdb}.dat'
        print(G_file)
        f = open(f'{G_file}')
        C = f.readlines()
        f.close()
        os.system('rm -rf tmp/')
        C0 = [C[k].strip().split() for k in range(8,len(C))]
        for CC in C0:
            gmax = np.array(CC[0].split(','), dtype=int)
            CC[0] = np.abs(gmax).max()
        C0 = np.array(C0, dtype=np.float32)
        G_number_list.append(C0)
    np.save('%s_Entanglement.npy'%mutant_type, G_number_list)
    print(f'Saved: {mutant_type}_Entanglement.npy')

################################## MAIN #######################################
for mutant_type in mutant_type_list:
    #outfiles/G/wt/G_1_i33_post_trans_1yta.dat
    # Get Qbb trajectory
    get_Qbb_act()

    # Get Entanglement trajectory
    get_entanglement()

print('NORMAL TERMINATION')
