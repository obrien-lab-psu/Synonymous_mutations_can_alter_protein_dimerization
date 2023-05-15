----------------------------------------------
Explainations of continuous synthesis protocol
----------------------------------------------

The "charmm_scripts" folder contains:
  -rnc_syn_dyna_flex.inp : run dynamics of nascent chain at a current length
  -rnc_syn_elong_mini_flex.inp: insert new aa into nascent chain length and minimize the chain 
The "inpfiles" folder contains:
  -initial structure information of the nascent chain with one amino acid located at PTC: 2is3_r1_rnc_complex.psf, 2is3_r1_rnc_complex.cor
  -topology files for combination of protein 2is3 and ribosome: 2is3_rnc_complex.prm, 2is3_rnc_complex.top
  -coordinate and protein structure information of ribosome: ribo_no_l22.cor, ribo_no_l22.psf  
  -sequence information: 2is3_mrna_sequence.txt
  -translation speeds of triple codons : fluitt_trans_times_mean_840000.txt
 
