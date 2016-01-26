# eLIFE-NPM_NMRrelaxation
All NMR related data analysis scripts for the eLife publication titled NPM1 integrates within the nucleolus via multi-modal interactions with proteins displaying R-rich linear motifs and rRNA by Mitrea DM et al., (2016)

Within the enclosed zipped file is a readme containing file descriptors and all scripts utilized for analyzing the NMR relaxation data recorded at 600, 800, and 1000 MHz. 

Also included is a python script for querying R-motifs from amino acid sequences (Multivalent_RMotifs).

###David Ban  Jan. 26, 2016 #####

Source_Code 1:

Python scripts used to determine individual and Global Kd (Kd_fits.py & Kd_lib.py) from NMR based TROSY-HSQC spectra.
Mathematica (global_KD_2state.nb) script provided to determine global Kd. 

Source_Code 2:

Determination of parameters that describe NMR fast (pico-nanosecond) relaxation using cross-correlated relaxation
and longitudinal relaxation rates. 

tC_indres_grid.py:  Initial grid search for individual residues determining tauC and S2 and tauF locally.
Library file called fastNMR_multifield_lib.py

global_tC_determination.nb:  Determines global tC for all residues and minimizes local S2 and tF.

local_param_determination.nb:  Determines local S2 and tF when tC is known from global minimzation.  Also performs MC error calculation.

flexible_param_determination.nb:  Determines tC_local and S2 for disordered tails on an individual residue basis.  Also performs MC error calculation.
