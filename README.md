# Ising_programs

A suite of python programs to analyze repeat-protein unfolding data with a 1D Ising model.

This repository contains four main folders. The two main folders (```heteropolymer_fit``` and ```homopolymer_fit```) contain scripts to run 1D Ising analysis.
One is for "capped homopolymeric" NRC-type repeats, and the other is for capped heteropolymeric NRXC-type repeats.
Both perform nonlinear least-squares fits generated plots, and perform and statistical analysis using boostrap resampling.

In addition there is a folder (```data_conversion```) with some additional scripts to convert and merge data files of different types.

Finally, there is a folder for calculation of partial correlation coefficients from bootstrapped thermodynamic parameters (```partial_correlation```).

All folders contain data files which contain protein folding data on which the scripts can be run.  All programs were written to run in python 3.8.

A detailed description of this suite of programs and its applications will soon be submitted to the journal Protein Science
for publication.  A preprint can be obtained from Doug Barrick (barrick@jhu.edu).

## License
[License](LICENSE.txt)

## Setup
To run the code, follow the instructions below:
1. git clone the repository:

    ```git clone https://github.com/barricklab-at-jhu/Ising_programs.git```
1. Create and activate a new environment (requires conda or miniconda):
    ```bash
    conda env create -f environment.yml
    conda activate ising_py3
    ```
If one is more comfortable working in ```virtualenv```, a ```requirements.txt``` also provided.

## Quikstart
1. To run the code via the self-contained notebooks (after setup above):
   * ```jupyter notebook```
   * navigate to```homopolymer_fit/ising_fitter_homopolymer_fit.ipynb``` OR ```heteropolymer_fit/ising_fitter_heteropolymer_fit.ipynb```

2. To run each program independently (after setup above):
   * Each script can be executed independently in the directory it resides:
      * e.g. ```python data_conversion.py```
      