"""
Converts Aviv.dat files into numpy files for NRC capped homopolymer repeats.

Most recent revision: 05/21/2020, Doug Barrick

"""

import numpy as np
import pandas as pd
import ntpath  # Good for path manipulations on a PC?
import glob  # Allows for unix-like specifications paths, using *, ?, etc.
import csv
import json
import time

start = time.time()

path = "/Users/dougbarrick/OneDrive - Johns Hopkins/Manuscripts/Ising_program/\
Scripts/scripts_vanilla_NRC/cANK_scripts/cANK_spyder_scripts/"
proj_name = "cANK"

den_nsig_const_melt = []

constructs = []  # List of constructs used to build partition functions and
# frac_folded expressions in next script.

melts = []  # List of melts to be used in fitting.

# Create an empty pandas dataframe to output a csv file from.
den_nsig_const_melt_df = pd.DataFrame(
    columns=["denat", "signal", "construct_melt", "dataset"]
)

# Gets file names, and extracts information including construct name, melt number.
num = 0
for filename in glob.glob("{}*.dat".format(path)):
    num = num + 1
    base = ntpath.basename(filename)
    melt = base.split(".")[0]
    construct = melt[:-2]

    # Reads the data portion of Aviv file, sticks it in a list, then normalizes
    # the y values and writes out a data file including construct and
    # a melt number to use as an ID for fitting in ising script.
    with open(filename, "r") as f:
        lines = (
            f.read().splitlines()
        )  # define the beginning and end of the data
        begin = 0
        end = 0
        while not lines[begin] == "$DATA":
            begin = begin + 1
        begin = begin + 4
        while not lines[end] == "$ENDDATA":
            end = end + 1
        xylist = []
        xyarray = []
        for row in range(begin, end - 1):  # extract the [denat] and CD signal
            line = lines[row]
            n = line.split()
            xylist.append([float(n[0]), float(n[1])])
        xyarray = np.array(xylist)

        # Below, the data is normalized.
        maxval = max(xyarray[:, 1])
        minval = min(xyarray[:, 1])
        normylist = []
        for i in range(0, len(xyarray)):
            normy = float(((xyarray[i, 1] - maxval) / (minval - maxval)))
            normylist.append(normy)
        for i in range(0, len(xylist)):
            den_nsig_const_melt.append(
                [xyarray[i, 0], normylist[i], construct, num]
            )

        # Build a numpy array for each melt and output for Ising fitter.
        # Columns are denaturant, normalized CD, construct, melt number.
        single_melt_dncm = []
        for i in range(0, len(xylist)):
            single_melt_dncm.append(
                [xyarray[i, 0], normylist[i], construct, num]
            )
        melt_array = np.array(single_melt_dncm)
        np.save(
            path + melt, melt_array
        )  # Writes an npy file to disk for each melt.
        temp_df = pd.DataFrame(melt_array)
        den_nsig_const_melt_df = den_nsig_const_melt_df.append(temp_df)
        if construct not in constructs:
            constructs.append(construct)
        melts.append(melt)

den_nsig_const_melt_df.to_csv(
    "{}{}_combined_data.csv".format(path, proj_name), index=False, header=False
)

# This loop puts melts in order of type (NRxC, NRx, RxC) and length.
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
melts = NRClist + NRlist + RClist

# This loop generates a construct list in the same order as the melts list.
for melt in melts:
    construct = melt[:-2]
    if construct not in constructs:
        constructs.append(construct)

# Write out the results.
with open("{0}{1}_constructs.txt".format(path, proj_name), "wb") as r:
    json.dump(constructs, r)

with open("{0}{1}_melts.txt".format(path, proj_name), "wb") as s:
    json.dump(melts, s)

stop = time.time()
runtime = stop - start
print("\nThe elapsed time was " + str(runtime) + " sec")
