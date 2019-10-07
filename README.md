# tight-binding
TB model for halide perovskites

This is a set of programs for converting moving between MD, TB, and DFT / wannier functions for MAPbI3 tight binding. 
The most important files are probably lammps_analysis_tools.py and cell_geometry.py

lammps_structure_tools.py
Tools for initializing structures and lammps, and defining structures for analysis later on. Key feature
is the structure abstract data type, and several functions which allow for creating lattices of MAPbI3 crystals. 

lammps_script_tools.py
Allows for writing the MD input files for runs of various types. These runs can be initiated with a 
structure, restart, or (I think) a dump file, and runs generally include an annealing phase and a phase 
which generates the trajectories. There is a sloppy trick to find the coulomb potential which involves 
setting the van der waals interactions to 0 and running a 0.0 duration "timestep". 

expiso.py 
Creates structures for different temperatures with fixed experimental lattice parameter ratios. Uses 
lammps_structure_tools.py and lammps_script_tools.py heavily. 

lammps_analysis_tools.py
Implements some tools for interpretting lammps ouput data for structure trajectories. Currently only uses 
data from dumps, but this could be changed. Implements the dump class and timestep class, which allow for 
extracting the position and potential (if charges change, they should be included here too) of an atom at 
an arbitrary timestep. Also allows for extracting box size / tilting, if these change.

index_conversions.py 
A bunch of tools for converting between stoichiometric unit cell geometry and orthorhombic unit cell 
geometry. Allows for retrieving the atom (ie in a dump) at some geometry, for example, finds the I atom
in the x direction relative to the Pb atom in the unit cell (4,1,2). Note: this would be a good place
to allow for structures which aren't periodic, since calculated neighbors uses the periodicity. 
Unfortunately, I think some of the functions that use index_conversion.py pass it already periodic 
values (ie, they do nx mod nx_max, in addition to the index conversion function doing this), so this 
might need to be cleaned up. 

tb_input_tools.py
Given a dump, writes the tight binding hopping parameters used to construct the hamiltonian. However, it
is hard coded for the old model and therefore not very useful. 

Note: cell_geometry.py is the main program for getting tight binding information, but it does not 
currently write the hopping parameters, only the features used to find them. The features can be converted 
to hopping parameters and written into the file using tb_parameter_analysis.jl. This system works but 
should be updated for simultaneous MD / TB. 

cell_geometry.py
Given a dump, calculates structural features relevant to each hopping parameter. For example, bond lengths,
angles, potentials, etc. Currently outputs these in a sequence for each parameter, which is the format 
used for fitting the model. 

tb_parameter_analysis.jl
(Note, this is in julia). Tools for reading in hopping data, performing statistical analysis, and
calculating spectral functions. Also allows for converting lists of features (outputed by cell_geometry.py)
to hopping parameters given a fit, and writing these into the format for the tight binding input. 

tb_analysis_tools.py
Some tools for analyzing and plotting tb output band structures (finding gaps, vbm, cbm, masses, etc).

Tools related to DFT, wannier, and fitting the model, but not necessary to run it: 

dft_input_tools.py
Takes a MD dump (from analysis_tools) and writes files to run dft and wannierization on the structure,
including scf input, nscf input, all the wannier inputs, and job scripts for dft and wannier. 

wan_hopping_tools.py
Allows for interpretting wannier output files in terms of the structure and various hopping parameters. 

run_wan_hop.py
Not sure if I uploaded this to github, but it runs the wannier_hopping_tools.py on wannier hamiltonian files
to produce files of hopping values, onsite values, and next nearest neighbor hopping. 

ml_tools.jl
Allows for doing various linear fits of features or functions of features. Not that related to ML, 
but it was used in conjunction. I want to implement some simple nonlinear fits. 

struct_to_dump.py
Uses the indices of the wannier orbitals to match dft wannier data with structural data from 
cell_geometry.py so that the wannier values can be fit with the structural features. (Also contains 
tools for turning Liang's old DFT input files into dumps which can be read by the rest of my programs)


