-----------------------------------------------------
Explainations of interaction energy calculation
-----------------------------------------------------

The charmm_scripts/inte_ene.inp will read in topology and protein structure informations of dimer complex: *.psf, *.top, *.prm . A loop over 500 dimer samplings from 500 annealing simulations of each dimer complex and write out informations: trajactory index, the distance between two monomers, total potential energy, interaction energy (line 50)   

The folder /inpfiles contains a sample trajactory of the annealing dimer structure 100th: anneal_0_complex_100_100_traj1.cor

For each dimer structure, the lowest interaction energy value from the set of 500 simulation trajactories will be collected.
The averaged interaction energy of 200 lowest-energy values corresponding to 200 dimer structures for each mRNA sequence is plotted on Figure2d and 2h. 
