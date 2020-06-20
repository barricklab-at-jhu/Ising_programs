"""
Generates a partition function and fraction folded expressions for an NRXC
capped heteropolymer repeat array.
"""

from __future__ import division
import sympy as sp
import numpy as np
import json
import time
import os

PATH = os.path.dirname(os.path.abspath(__file__))

proj_name = "T4V_NRC_2mi"

start = time.time()

print("\nGenerating partition functions and fraction folded expressions...")

# Parameters for partition function calculation.  Note these are sympy symbols.
RT = sp.Symbol("RT")
dGN = sp.Symbol("dGN")
dGR = sp.Symbol("dGR")
dGX = sp.Symbol("dGX")
dGC = sp.Symbol("dGC")
dGRR = sp.Symbol("dGRR")
dGXX = sp.Symbol("dGXX")
dGXR = sp.Symbol("dGXR")
dGRX = sp.Symbol("dGRX")
mR = sp.Symbol("mR")
mX = sp.Symbol("mX")
denat = sp.Symbol("denat")
KN = sp.Symbol("KN")
KR = sp.Symbol("KR")
KX = sp.Symbol("KX")
KC = sp.Symbol("KC")
TRR = sp.Symbol("TRR")
TXX = sp.Symbol("TXX")
TXR = sp.Symbol("TXR")
TRX = sp.Symbol("TRX")

exp = sp.Function("np.exp")

with open(os.path.join(PATH, f"{proj_name}_constructs.json"), "r") as cons:
    constructs = json.load(cons)

# define weight matricies and end vectors to be used to calculate partition functions
begin = sp.Matrix([[0, 1]])
woN = sp.Matrix(
    [[KN, 1], [KN, 1]]
)  # Leave off coupling term.  No zeroth repeat
woR = sp.Matrix(
    [[KR, 1], [KR, 1]]
)  # Leave off coupling term.  No zeroth repeat
woX = sp.Matrix(
    [[KX, 1], [KX, 1]]
)  # Leave off coupling term.  No zeroth repeat

wnR = sp.Matrix(
    [[KR * TRR, 1], [KR, 1]]
)  # Treat coupling same as an rR interface, as usual
wrR = sp.Matrix([[KR * TRR, 1], [KR, 1]])  # Same as wnR above
wxR = sp.Matrix([[KR * TXR, 1], [KR, 1]])

wnX = sp.Matrix([[KX * TRX, 1], [KX, 1]])  # Different coupling than
wrX = sp.Matrix([[KX * TRX, 1], [KX, 1]])  # Different coupling than wnX above
wxX = sp.Matrix([[KX * TXX, 1], [KX, 1]])

wrC = sp.Matrix(
    [[KC * TRR, 1], [KC, 1]]
)  # Treat copuling same as an rR interface, as usual
wxC = sp.Matrix([[KC * TXR, 1], [KC, 1]])

end = sp.Matrix([[1], [1]])

# Build a dictionary of these matrices
w_dict = {
    "oN": woN,
    "oR": woR,
    "oX": woX,
    "nR": wnR,
    "rR": wrR,
    "xR": wxR,
    "nX": wnX,
    "rX": wrX,
    "xX": wxX,
    "rC": wrC,
    "xC": wxC,
}

# Build dictionaries of partition functions, partial derivs with respect
# to K, and fraction folded.

q_dict = {}
dqdKN_dict = {}
dqdKR_dict = {}
dqdKX_dict = {}
dqdKC_dict = {}
frac_folded_dict = {}

# Number of repeats of each type.  Seems like they should be floats, but
# I get an error in the matrix multiplication (q_dict) if they are declared to be.


for construct in constructs:

    # Make partition function dictionary and expressions for fraction folded.
    # Note, only one pf is generated per construct, even when there are multiple melts.

    repeat_list = construct.split("_")
    repeat_list.insert(0, "o")

    pairs_list = []
    i = 1
    while i < len(repeat_list):
        pair = (
            repeat_list[i - 1].lower() + repeat_list[i]
        )  # Need to convert prev rept to lower case
        pairs_list.append(pair)
        i = i + 1

    q_dict[construct + "_q"] = begin

    for pair in pairs_list:
        q_dict[construct + "_q"] = q_dict[construct + "_q"] * w_dict[pair]

    q_dict[construct + "_q"] = q_dict[construct + "_q"] * end

    # Next two lines convert from sp.Matrix to np.array to something else.
    # Not sure the logic here, but it works.
    q_dict[construct + "_q"] = np.array(q_dict[construct + "_q"])
    q_dict[construct + "_q"] = q_dict[construct + "_q"].item(0)

    # Partial derivs wrt KN dictionary.
    dqdKN_dict[construct + "_dqdKN"] = sp.diff(q_dict[construct + "_q"], KN)

    # Partial derivs wrt KR dictionary.
    dqdKR_dict[construct + "_dqdKR"] = sp.diff(q_dict[construct + "_q"], KR)

    # Partial derivs wrt KX dictionary.
    dqdKX_dict[construct + "_dqdKX"] = sp.diff(q_dict[construct + "_q"], KX)

    # Partial derivs wrt KC dictionary.
    dqdKC_dict[construct + "_dqdKC"] = sp.diff(q_dict[construct + "_q"], KC)

    # Fraction folded dictionary.
    frac_folded_dict[construct + "_frac_folded"] = (
        KN / (q_dict[construct + "_q"]) * dqdKN_dict[construct + "_dqdKN"]
        + KR / (q_dict[construct + "_q"]) * dqdKR_dict[construct + "_dqdKR"]
        + KX / (q_dict[construct + "_q"]) * dqdKX_dict[construct + "_dqdKX"]
        + KC / (q_dict[construct + "_q"]) * dqdKC_dict[construct + "_dqdKC"]
    ) / (len(pairs_list))

# The loop below replaces K's and W's the fraction folded terms in the
# dictionary with DGs, ms, and denaturant concentrations.  The simplify line
# is really important for making compact expressions for fraction folded.
# This simplification greatly speeds up fitting.  The last line
# converts from a sympy object to a string, to allow for json dump.

for construct in frac_folded_dict:
    frac_folded_dict[construct] = frac_folded_dict[construct].subs(
        {
            KN: (exp(-((dGN + (mR * denat)) / RT))),
            KR: (exp(-((dGR + (mR * denat)) / RT))),
            KX: (exp(-((dGX + (mX * denat)) / RT))),
            KC: (exp(-((dGC + (mR * denat)) / RT))),
            TRR: (exp(-dGRR / RT)),
            TXX: (exp(-dGXX / RT)),
            TXR: (exp(-dGXR / RT)),
            TRX: (exp(-dGRX / RT)),
        }
    )
    frac_folded_dict[construct] = sp.simplify(frac_folded_dict[construct])
    frac_folded_dict[construct] = str(frac_folded_dict[construct])

with open(os.path.join(PATH, f"{proj_name}_frac_folded_dict.json"), "w") as f:
    json.dump(frac_folded_dict, f)

#  The next code block caluclates the rank of the coefficient matrix
#  and returns it to the user.

num_constructs = len(constructs)
thermo_param_list = ['dGN','dGR','dGX','dGC','dGRR','dGXX','dGXR','dGRX']
num_params = len(thermo_param_list)

coeff_matrix = np.zeros((num_constructs, num_params))

row = 0
for construct in constructs:
    repeats_list = construct.split('_')
    # For first repeat, don't need (or want) to worry about interface.
    i = 0
    if repeats_list[i] == 'N':
        coeff_matrix[row, 0] = 1
    elif repeats_list[i] == 'R':
        coeff_matrix[row, 1] = coeff_matrix[row, 1] + 1
    elif repeats_list[i] == 'X':
        coeff_matrix[row, 2] = coeff_matrix[row, 2] + 1
    else: 
        coeff_matrix[row, 3] = 1
    
    # For remaining repeats, need to worry about interfaces
    i = 1   
    while i < len(repeats_list):
        if repeats_list[i] == 'R':
            coeff_matrix[row, 1] = coeff_matrix[row, 1] + 1
        elif repeats_list[i] == 'X':
            coeff_matrix[row, 2] = coeff_matrix[row, 2] + 1
        else: 
            coeff_matrix[row, 3] = 1
        
        if repeats_list[i] == 'R':
            if (repeats_list[i-1] == 'N' or repeats_list[i-1] == 'R'):
                coeff_matrix[row, 4] = coeff_matrix[row, 4] + 1  # RR-type interface (includes NR)
            else:
                coeff_matrix[row, 6] = coeff_matrix[row, 6] + 1  # XR interface

        if repeats_list[i] == 'X':
            if (repeats_list[i-1] == 'N' or repeats_list[i-1] == 'R'):
                coeff_matrix[row, 7] = coeff_matrix[row, 7] + 1  # RX-type interface (includes NX)
            else:
                coeff_matrix[row, 5] = coeff_matrix[row, 5] + 1  # XX interface

        if repeats_list[i] == 'C':
            if (repeats_list[i-1] == 'R'):
                coeff_matrix[row, 4] = coeff_matrix[row, 4] + 1  # RR-type interface (RC)
            else:
                coeff_matrix[row, 6] = coeff_matrix[row, 6] + 1  # XR-type interface (XC)
            
        i =i + 1
    row = row + 1
        
rank = np.linalg.matrix_rank(coeff_matrix)

if rank == num_params:
    print("\nThe coefficeint matrix has full column rank (r=",rank,")") #leaves a space betw rank and ).  Not sure why.
else:
    print("\nThe coefficeint matrix has incomplete column rank (r=",rank,").") 
    print("You should revise your model or include the necessary constructs to obtain full rank.\n")

stop = time.time()
runtime = stop - start
print("\nThe elapsed time was " + str(runtime) + " sec")
