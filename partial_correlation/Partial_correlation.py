'''
Calculates partial correlation coefficeints.  

Most recent version 2020/05/21 D. Barrick

This needs to be run in python 3 because Pingouin will not run in python 2.

'''


import numpy as np
import pandas as pd
import pingouin as pg

path = '/Users/dougbarrick/OneDrive - Johns Hopkins/Manuscripts/Ising_program/\
Scripts/scripts_mutations/T4V_two_mi_values_2020_05_10/'

proj_name = 'T4V_NRC_2mi'

bs_params_df = pd.read_csv('{}{}_bootstrap_params.csv'.format(path, proj_name), index_col='Bootstrap Iter')
del bs_params_df['redchi**2'], bs_params_df['bestchi**2']

bs_params_partial_corr_df = bs_params_df.pcorr() # Gives matrix of partial correlations, where all other variables are factored out
bs_params_partial_corr_df.to_csv('{}{}_bootstrap_partial_corr.csv'.format(path, proj_name))