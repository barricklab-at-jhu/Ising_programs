# Ising_programs

May 23, 2020
Doug Barrick

A suite of python programs to analyze repeat-protein unfolding data with a 1D Ising model.

This repository contains four folders.  The two main folders contain scripts to run 1D Ising analysis.
One is for "capped homopolymeric" NRC-type repeats, and the other is for capped heteropolymeric NRXC-type repeats.
Both perform nonlinear least-squares fits generated plots, and perform and statistical analysis using boostrap resampling.

In addition there is a folder with some additional scripts to convert and merge data files of different types.

Finally, there is a folder for calculation of partial correlation coefficients from bootstrapped thermodynamic parameters.

All folders contain data files on which the scripts can be run.  All programs were written to run in python 2.7 except 
for the partial correlation program, which requires python 3.7.

A detailed description of this suite of programs and its applications will soon be submitted to the journal Protein Sciene
for publication.  A preprint can be obtained from Doug Barrick (barrick@jhu.edu).
