All-atom simulations GROMACS
------------------------------

This README file explains the input structures and mdp files for running molecular dynamics simulations of protein dimers and monomers in GROMACS.

---------
* The `input_structures` folder contains eight pdb files: five for dimer simulations (d1.pdb-d5.pdb) and three for monomer simulations (m1.pdb-m3.pdb). These are the initial structures for the simulations.
---------
* The `mdp_files` folder contains five configuration files for different steps of the simulations in GROMACS. 
The files are used in the following order:

  1. `em.mdp`: This file performs energy minimization on the input structures to remove any steric clashes or bad contacts.
  2. `nvt.mdp`: This file runs a 1ns equilibrium simulation in the NVT ensemble with position restraints applied to all heavy atoms of the protein. This step ensures that the temperature of the system is stable.
  3. `npt.mdp`: This file runs a 1ns equilibrium simulation in the NPT ensemble with position restraints applied to all heavy atoms of the protein. This step ensures that the pressure and density of the system are stable. The simulation continues from the NVT step.
  4. `npt2.mdp`: This file runs a second 1ns equilibrium simulation in the NPT ensemble with position restraints applied to only the C-alpha atoms of the protein. This step allows more flexibility for the protein backbone and side chains. The simulation continues from the NPT step.
  5. `md.mdp`: This file runs a 500ns production simulation in the NPT ensemble without any position restraints. This is the final step of the simulation where the protein dynamics are observed.
