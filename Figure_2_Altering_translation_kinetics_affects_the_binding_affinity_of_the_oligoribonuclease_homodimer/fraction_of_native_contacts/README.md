The six files

1yta_wt_ave_Q.dat, 1yta_slow_ave_Q.dat, 1yta_fast_ave_Q.dat

and

2is3_wt_ave_Q.dat, 2is3_slow_ave_Q.dat, 2is3_fast_ave_Q.dat

contain the ensemble average fraction of native contacts as a function of time since release from the ribosome
across the set of 200 trajectories for each mRNA. Each file contains two columns. The first column is time in microseconds, and
the second column is the average fraction of native contacts at that time. 

The Python3 program plot_meanQ.py can be used to plot the results in these files to generate 1yta_Q.png and 2is3_Q.png, which are
displayed as Figures 2 c and g in the main text, respectively. 
