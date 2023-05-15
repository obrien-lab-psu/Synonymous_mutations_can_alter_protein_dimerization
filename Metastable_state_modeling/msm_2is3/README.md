# Metastable state analysis  
1. Q analysis
2. ent analysis
3. get OP data
4. build post-trans model

NOTE: This workflow was applied to 200 statistically independant MD trajectories for each synonymous variant. Due to size contraints only a single DCD for each variant is present but the files in the outfiles directory are those from the analysis applied to the full 200. 


### Fraction of native contact analysis 

#### Usage: perl calc_native_contact_fraction.pl
~~~
              --input | -i <INPUT.COR> for identify native contacts
              --domain | -d <DOMAIN.DAT> for domain defination
              [--secs | -s] <SECONDARY STRUCTURE> for secondary 
                          structure defination. If no file specified,
                          it will calculate all native contacts regardless
                          of the secondary structure.
              --traj | -t <TRAJ.DCD> for simulation trajectory
              [--begin | -b] <START STEP> to split trajectory. 
                                  Default is 1
              [--end | -e] <END STEP> to split trajectory. 
                                Default is the end frame.
              [--meaningful | -m] <1 or 0> 1 to only calculate contacts
                                  in the interface has native contacts;
                                  0 to calculate all possible interfaces.
                                  The interface has no native contacts will
                                  always has Q = 1.
                                  Default is 0. 
              [--mask | -k] <MASK> for a subset of residues to calculate their
                                   native contacts with all other residues. 
                                   Default is 'all'.
                                   Example MASK: 
                                   '1-45' for selection from resid 1 to 45;
                                   '1,3,5-9' for selection of resid 1, 3 
                                   and from 5 to 9.
              [--outdir | -o] <DIRECTORY> for the outputs. Default is the directory
                                          of your trajectories
              [--help | -h]

  Example for DOMAIN.DAT:
  1:96 302:350 a #Domain 1 is from resid 1 to 96 and 302 to 350, and in
                 #alpha-helix class
  97:155 b #Domain 2 is from resid 97 to 155 and in beta-sheet class
  156:301 c #Domain 3 is from resid 156 to 301 and in alpha-beta class
~~~

##### Example command  
~~~ 
perl codes/calc_native_contact_fraction.pl -i inpfiles/2is3_xray_a_ca.cor -d inpfiles/domain.dat -s inpfiles/dimer_2is3A.secstructs -t post_translation/wt/1-200/1/1_i33_post_trans_2is3.dcd  -k 1-181 -o post_translation/wt/
~~~ 


### Fraction of native contact with changes in entanglement analysis

#### Usage: perl calc_entanglement_number.pl
~~~
              --input | -i <INPUT.COR> for identify native contacts
              --traj | -t <TRAJ.DCD> for simulation trajectory
              [--begin | -b] <START STEP> to split trajectory. 
                                  Default is 1
              [--end | -e] <END STEP> to split trajectory. 
                                Default is the end frame.
              [--mask | -k] <MASK> for a subset of residues to calculate their
                                   native contacts with all other residues. 
                                   Default is 'all'.
                                   Example MASK: 
                                   '1-45' for selection from resid 1 to 45;
                                   '1,3,5-9' for selection of resid 1, 3 
                                   and from 5 to 9.
              [--outdir | -o] <DIRECTORY> for the outputs. Default is the directory
                                          of your trajectories
              [--help | -h]

~~~

##### Example command
~~~
perl codes/calc_entanglement_number.pl -i inpfiles/2is3_xray_a_ca.cor -t post_translation/wt/1-200/1/1_i33_post_trans_2is3.dcd -o outfiles/G/wt/
~~~


### Get order parameter files organized and ready for building meta stable states

#### Usage: python get_post_trans_order_parameters_v1.1.py 

python get_post_trans_order_parameters_v1.1.py [1] [2]
[1] control file path
[2] PDBID

Control File Contents:
n_traj = # of trajectories in each variant
mutant_type_list = fast slow wt
act_mask = Amber style mask for atom selection
prefix_dir = directory containing OP directories G and Q

##### Example command
~~~
python codes/get_post_trans_order_parameters_v1.1.py inpfiles/get_pt_op.cntrl 2is3
~~~


### Build the metastable state model across the post translational simulations 

#### Usage: python codes/build_post_trans_kinetic_model.py 

python codes/build_post_trans_kinetic_model.py [1] [2]
[1] control file path
[2] PDBID

##### 3. Parameters in Control file
| Keywords | Type | Intepretation |
| ------ | ------ | ------ |
| end_t | float | The end time (in seconds) for the Master equation model estimation. Default is 60 s. |
| dt | float | Time step (in nanoseconds) of the CG simulation. Default is 0.000015 ns. |
| nsave | int | Number of steps of saving one frame in the trajectory. Default is 5000. |
| alpha | float | The scaling factor ($`\alpha`$) used to convert the experimental time scale (*in vivo*, $`\tau^{\text{exp}}`$ in s) to simulation time scale ($`\tau^{\text{sim}}`$ in ns) by $`\tau^{\text{sim}} = \frac{1}{\alpha}\tau^{\text{exp}} \times 10^9`$. Default is 4331293.0. |
| n_window | int | Number of the frames to calculate the mode value. Default is 200. |
| n_traj | int | Total number of co-translational trajectories for each mutant. Default is 100. |
| mutant_type_list | list of string | Name of the mutants. Must be seperated by a white space. Default is 'fast slow'. |
| n_cluster | int | Number of k-means clusters to group. Default is 400. |
| stride | int | Stride of reading trajectory frame when clustring by k-means. Default is 10. |
| n_large_states | int | Maximum number of metastable states to be assigned for the largest connective subgraph. Default is 10. |
| n_small_states | int | Maxumum number of metastable states to be assigned for the rest subgraph (non-connective). Default is 2. |
| lag_t | int | Lag frame for building the markov state model when assigning the metastable states. Default is 1. |
| start_idx | int | Start trajectory index of replicated trajectories. Default is 1.  |
| end_idx | int | End trajectory index of replicated trajectories. Default is 10. |
| sample_size | int | Number of representative structures to be sampled for each metastable state. Default is 5. |
| native_AA_pdb | string | Absolute path to the all-atom native pdb file of the CG protein for backmapping. |
| prefix_dir | string | relative path to the directory that contains the `continuous_synthesis/` and `post_translation/` folders. |
| OP_dir | string | relative path to the directory that contains the OP .npy files created by codes/get_post_trans_order_parameters_v1.1.py |
| visualiz_threshold | float | The threshold of G value to determine whether or not to generate the entanglement visualiztion in the `vmd tcl` scripts. Default is 0.02 (recommended). |
| if_cluster | 0 or 1 | 1: Run k-means clustering; 0: Skip running k-means clustering and read cluster data from `msm_data.npz`. (Must have this file created.) Default is 1. |
| if_visualize | 0 or 1 | 1: Backmap and Generate the `vmd tcl` scripts to visualize the first representative structure of each metastable states; 0: Skip backmapping and generating the `vmd tcl` scripts. Default is 1. |
| if_sample | 0 or 1 | 1: Sample the trajetories to get representative structures; 0: Skip sampling the representative structures and read the representative structure data from `msm_data.npz`. (Must have this file created.) Default is 1. |

##### 4. Example of Control file
```
end_t = 60
dt = 0.000015
nsave = 5000
alpha = 3967486.0 
n_window = 200
n_traj = 200
mutant_type_list = fast slow wt
n_cluster = 400
stride = 10
n_large_states = 6
n_small_states = 2
lag_t = 1
start_idx = 1
end_idx = 1
sample_size = 5
nbins = 15 15
native_AA_pdb = inpfiles/mini_2is3_A.pdb
prefix_dir = ./
OP_dir = outfiles/
visualiz_threshold = 0.02
if_cluster = 1
if_visualize = 1
if_sample = 1              

```

##### Example command
~~~
python codes/build_post_trans_kinetic_model.py inpfiles/build_pt_kmodel.cntrl 2is3
~~~
