# Ising_programs

A suite of python programs to analyze repeat-protein unfolding data with a 1D Ising model.

This repository contains four folders.  The two main folders contain scripts to run 1D Ising analysis.
One is for "capped homopolymeric" NRC-type repeats, and the other is for capped heteropolymeric NRXC-type repeats.
Both perform nonlinear least-squares fits generated plots, and perform and statistical analysis using boostrap resampling.

In addition there is a folder with some additional scripts to convert and merge data files of different types.

Finally, there is a folder for calculation of partial correlation coefficients from bootstrapped thermodynamic parameters.

All folders contain data files on which the scripts can be run.  All programs were written to run in python 3.8.

A detailed description of this suite of programs and its applications will soon be submitted to the journal Protein Science
for publication.  A preprint can be obtained from Doug Barrick (barrick@jhu.edu).

## License
[License](LICENSE.txt)

## Setup
```
git clone https://github.com/barricklab-at-jhu/Ising_programs.git
```
1. Create and activate a new environment (requires conda or miniconda):
    ```bash
    conda env create -n "ising_py3" python=3.8
    conda activate ising_py3
    ```
2. Install from ```setup.py```:
   ```bash
   pip install .
   ```

## Quikstart
1. To run the code via the self-contained notebooks (after setup above):
   * ```jupyter notebook```
   * navigate to```examples/homopolymer_fit/Ising_fitter_heteropolymer.ipynb``` OR ```examples/heteropolymer_fit/Ising_fitter_homopolymer.ipynb```
