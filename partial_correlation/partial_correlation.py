"""
Calculates partial correlation coefficeints.
"""

import pandas as pd
import os

PATH = os.path.dirname(os.path.abspath(__file__))

proj_name = "T4V_NRC_2mi"

bs_params_df = pd.read_csv(
    os.path.join(PATH, f"{proj_name}_bootstrap_params.csv"),
    index_col="Bootstrap Iter",
)
del bs_params_df["redchi**2"], bs_params_df["bestchi**2"]

bs_params_partial_corr_df = (
    bs_params_df.pcorr()
)  # Gives matrix of partial correlations, where all other variables are factored out
bs_params_partial_corr_df.to_csv(
    os.path.join(PATH, f"{proj_name}_bootstrap_partial_corr.csv")
)
