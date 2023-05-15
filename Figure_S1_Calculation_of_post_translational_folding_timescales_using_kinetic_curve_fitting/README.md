The directory fpt_analysis/ contains six subdirectories, one for each 1yta and 2is3 translation schedule used:

fpt_analysis_1yta_fast/
fpt_analysis_1yta_slow/
fpt_analysis_1yta_wt/
fpt_analysis_2is3_fast/
fpt_analysis_2is3_slow/
fpt_analysis_2is3_wt/

as well as four files containing confidence interval information:

1yta_final_Q_mean_and_CIs.txt
1yta_percent_unfolded_and_CIs.txt

2is3_final_Q_mean_and_CIs.txt
2is3_percent_unfolded_and_CIs.txt

First-passage time analysis was performed for each set of 200 trajectories run for a given protein and mRNA template. 
For example, in fpt_analysis_1yta_fast/ the command 'python fpt_analysisV4.1_new.py list_props.txt 66000 0.69080 5000 0.015 > frac_U_vs_time_200trajs.txt' 
was run to compute the survival probability of the unfolded state as a function of time. See the code for details. The individual time series of 
Q for each trajectory (e.g., fpt_analysis/fpt_analysis_1yta_fast/props_ts/props_t1.dat) were computed using XXX. The raw trajectory files are
large, so we include an example calculation for one frame of YYY in ZZZ. 

Each of the resulting survival probability curves was then plotted and fit to a double-exponential function with plot_time_series.py

The results of this plotting program were used to generate Figure S1, with the originals provided here as "figure_s1.png".

A copy of the output printed to screen by plot_time_series.py, which contains the fit results reported in Table S1, is provided as "fit_results.txt"





