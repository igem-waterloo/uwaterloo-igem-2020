# UWaterloo iGEM 2020 Project Files

In this project we are attempting to simulate protein-transition metal binding, with the
ultimate goal of calculating forward and reverse rate coefficients of binding in the 
mass-action model as well as to compare binding affinities of different candidate proteins
which have been selected via rational design principles.

In md.mdp you will find the production level code. We have followed the "Protein-Ligand Binding"
tutorial from gromacs that can be found at http://www.mdtutorials.com/gmx/complex/index.html.
This setup does not incorporate quantum effects directly, we use the Charmm36 parameters found
at https://www.charmm.org/charmm/resources/charmm-force-fields/ and http://www.charmm-gui.org/?doc=input/nanomaterial
which incorporate parameters found in https://pubs.acs.org/doi/full/10.1021/jp505875v and
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4306492/ and https://pubs.acs.org/doi/pdf/10.1021/ct400146w.
These parameters are calculated by the Lennard Jones 12-6-4 method.

# Todo List

MAYBE TODO: Use Gromacs QM/MM package which can be found at [1] http://wwwuser.gwdg.de/~ggroenh/qmmm.html.
This package hooks Gromacs up to Gaussian, which simulates the quantum mechanics of a portion
of the system. 
1. Select a number of residues which make up the active site where we want to simulate quantum effects.
This will likely just be the proteins which make up the active site geometry, but what about the atoms 
further away?
2. Turn this into a group by running `gmx make_ndx` which allows you to specify groups of residues and
name them. We can just name our group "Active_Site" and then use that in the commands.
