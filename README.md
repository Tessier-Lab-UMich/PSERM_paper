# PSERM Paper
### Matt Smith
### 5/2/2023

This repository contains all the datsets and analyses used to create the figures for PSERM paper. The numbered folders contain Jupyter Notebooks that were used to generate data and create the main figures. A conda environment is included which should easily allow for running any of the analyses. 

The main code for analyzing the raw fastq files from the sequencer and creating the sequence:frequency csv files (`ngs.py`) and the code for computing all PSERMs and scoring (`pserm.py`) is included withing this repository. `ngs.py` contains several classes used to compute and store the sequence and frequency data. These calsses are used as input for the `ngs_analysis` object of `pserm.py`. To use `ngs.py` and `pserm.py` with these datasets, refer to the Jupyter Notebooks to see examples of how these objects are called and initiated.

The deep sequencing datasets used for this paper can all be found within the `Input_Datasets` folder under each of the respective projects. These datasets only show the mutated positions of the antibody and not the full antibody sequences.

The `Project_PPMs` and `Project_PSSMs` folders contain the PPMs or PSSMs computed wihtin the notebooks. You can create them every time, however it's much quicker to save them and use the `load_ppm` or `load_pserm` functions.