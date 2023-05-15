The mRNA sequences for the wildtype, fastest-translating synonymous mutant, and slowest-translating synonymous
mutant for 1yta and 2is3 are provided in the files:

1yta_wt_mrna_sequence.txt
1yta_fast_mrna_sequence.txt
1yta_slow_mrna_sequence.txt

2is3_wt_mrna_sequence.txt
2is3_fast_mrna_sequence.txt
2is3_slow_mrna_sequence.txt

The codon-specific translation times (taken from Supplementary Table 7 of Nissley et al. JACS 2020) in units
of 0.015-ps integration time steps, are in the file fluitt_trans_times_mean_840000.txt

The Python3 program plot_ta.py reads in the mRNA sequences, matches translation times to codons, and then plots
a moving average of the normalized translation time. The raw values in units of integration times steps are divided by
the average translation time of 840000 to set the average time to 1. The results are saved in 

1yta_trans_schedule.png

and

2is3_trans_schedule.png

These images appear in the main text as Figure 2 panels b and f, respectively. 


