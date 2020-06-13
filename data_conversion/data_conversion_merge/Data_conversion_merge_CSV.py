"""
Combines csv files for different types of capped repeat arrays (NRC homopolymer
arrays and NRXC heteropolymer arrays with point substitions (here T4V in
consensus ankyrin repeat arrays, creates numpy files for each melt, and 
generates a melt list and a construct list.
"""

import numpy as np
import pandas as pd
import json
import os

PATH = os.path.dirname(os.path.abspath(__file__))
# move to package ising.data imports when complete


proj_name = "T4V_NRC_2mi"

den_nsig_const_melt = []
constructs = (
    []
)  # List of constructs used to build partition functions in next script
melts = []  # List of melts to be used in fitting script

T4V_input_df = pd.read_csv(
    os.path.join(PATH, "T4Vdata_not_normalized.csv"),
    names=["denat", "signal", "construct_melt", "dataset"],
)
NRC_input_df = pd.read_csv(
    os.path.join(PATH, "NRC_data_dnmn.csv"),
    names=["denat", "signal", "construct_melt", "dataset"],
)
maxT4Vmelt = T4V_input_df[
    "dataset"
].max()  # Finds the maximum number of melts in the first df
NRC_input_df["dataset"] = (
    NRC_input_df["dataset"] + maxT4Vmelt
)  # Adds the max number to the melt numbers in second df
combined_input_df = pd.concat(
    [T4V_input_df, NRC_input_df],
    names=["denat", "signal", "construct_melt", "dataset"],
)
combined_input_df.to_csv(
    os.path.join(PATH, "T4V_NRC_dnmn.csv"), index=False, header=False
)

num_melts = combined_input_df["dataset"].max()
for melt in np.arange(num_melts) + 1:
    temp_df = combined_input_df.loc[
        combined_input_df.dataset == melt
    ]  # Pulls out just one melt

    # Normalizing the signal
    min = temp_df["signal"].min()
    max = temp_df["signal"].max()
    series = (temp_df["signal"] - min) / (max - min)
    temp_df["signal"] = series  # Overwrites un-normalized signal
    temp_df.rename(columns={"signal": "nsig"}, inplace=True)
    temp_list = temp_df.values.tolist()
    temp_nparray = np.array(temp_list)
    construct_melt = temp_df.iloc[0, 2]
    np.save(
        PATH + construct_melt, temp_nparray
    )  # Writes an npy file to disk for each melt.
    melts.append(construct_melt)

""" 
This loop puts melts in order of type (NRxC, NRx, RxC) and length.  This is useful for the
plotting script below, putting the by_melt legends in a sensible order
"""
NRClist = []
NRlist = []
RClist = []
melts.sort()
i = 0
for melt in melts:
    if melt[0] == "N":
        if melt[-3] == "C":
            NRClist.append(melt)
        else:
            NRlist.append(melt)
    else:
        RClist.append(melt)

NRClist.sort(key=len)
NRlist.sort(key=len)
RClist.sort(key=len)
melts = NRClist + NRlist + RClist

# Generate a list of just the constructs.  The loop removes duplicates.
for melt in melts:
    if melt[:-2] not in constructs:
        constructs.append(melt[:-2])

with open(os.path.join(PATH, f"{proj_name}_constructs.json"), "w") as file:
    json.dump(constructs, file)

with open(os.path.join(PATH, f"{proj_name}_melts.json"), "w") as file:
    json.dump(melts, file)
