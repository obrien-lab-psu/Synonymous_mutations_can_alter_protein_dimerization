# Caclulate the relative percent change in the solvant expsure of key residues 

### Usage: MSM_SASA_calcs_v3.0.py 
~~~
python MSM_SASA_calcs_v3.0.py [1] [2] [3] [4]
[1] path to MSM file
[2] path to SASA files
[3] path to output file
[4] mask indexed from 0
[5] pdb
~~~

##### Example commands
~~~
python codes/MSM_SASA_calcs_v3.0.py inpfiles/msm_1yta/msm_data.npz sasa_data/ ./ 84-94 1yta
python codes/MSM_SASA_calcs_v3.0.py inpfiles/msm_1yta/msm_data.npz sasa_data/ ./ 166-175 1yta
python codes/MSM_SASA_calcs_v3.0.py inpfiles/msm_2is3/msm_data.npz sasa_data/ ./ 2-13 2is3
python codes/MSM_SASA_calcs_v3.0.py inpfiles/msm_2is3/msm_data.npz sasa_data/ ./ 204-3.0 2is3
~~~
