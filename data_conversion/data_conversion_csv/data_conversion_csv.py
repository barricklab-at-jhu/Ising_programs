"""
Converts capped homopolymer (NRC) arrays from a single csv file to  
numpy files for each melt, and generates a melt list and a construct list.

Most recent revision: 05/21/2020, Doug Barrick

"""

import numpy as np
import pandas as pd
import json
import os

PATH = os.path.dirname(os.path.abspath(__file__))

proj_name = "AnkT4V"

den_nsig_const_melt = []
constructs = (
    []
)  # List of constructs used to build partition functions in next script
melts = []  # List of melts to be used in fitting script

csv_input_df = pd.read_csv(
    "T4Vdata_not_normalized.csv",
    names=["denat", "signal", "construct_melt", "dataset"],
)

# df_dict = {}  # Not sure I will use this.

num_melts = csv_input_df["dataset"].max()
for melt in np.arange(num_melts) + 1:
    temp_df = csv_input_df.loc[
        csv_input_df.dataset == melt
    ]  # Pulls out just one melt

    # Normalizing the signal
    min = temp_df["signal"].min()
    max = temp_df["signal"].max()
    series = (temp_df["signal"] - min) / (max - min)
    temp_df[
        "signal"
    ] = series  # Overwrites un-normalized signal (and generates a warning below)
    temp_df.rename(columns={"signal": "nsig"}, inplace=True)
    temp_list = temp_df.values.tolist()
    temp_nparray = np.array(temp_list)
    construct_melt = temp_df.iloc[0, 2]
    melts.append(construct_melt)
    np.save(
        PATH + construct_melt, temp_nparray
    )  # Writes an npy file to disk for each melt.

    # Generate a list of just the constructs.  The loop removes duplicates.
    for melt in melts:
        if melt[:-2] not in constructs:
            constructs.append(melt[:-2])

with open("{0}{1}_constructs.json".format(PATH, proj_name), "w") as r:
    json.dump(constructs, r)

with open("{0}{1}_melts.json".format(PATH, proj_name), "w") as s:
    json.dump(melts, s)
