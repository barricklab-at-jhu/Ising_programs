"""
Generates partition function and fraction folded expressions for an NRC
capped homopolymer series, and generates a melts list and a constructs list.

Most recent revision: 05/21/2020, Doug Barrick

"""

from __future__ import division
import sympy as sp
import numpy as np
import json
import time

start = time.time()

print(
    "\nGenerating partition functions and fraction folded expressions.  This may take a minute..."
)

path = "/Users/dougbarrick/OneDrive - Johns Hopkins/Manuscripts/Ising_program/\
Scripts/scripts_vanilla_NRC/cANK_scripts/cANK_spyder_scripts/"
proj_name = "cANK"

# Parameters for partition function calculation.  Note these are sympy symbols.
RT = sp.Symbol("RT")
dGN = sp.Symbol("dGN")
dGR = sp.Symbol("dGR")
dGC = sp.Symbol("dGC")
mi = sp.Symbol("mi")
denat = sp.Symbol("denat")
Kn = sp.Symbol("Kn")
Kr = sp.Symbol("Kr")
Kc = sp.Symbol("Kc")
dGinter = sp.Symbol("dGinter")
W = sp.Symbol("W")

np.exp = sp.Function("np.exp")

with open("{0}{1}_constructs.txt".format(path, proj_name)) as cons:
    constructs = json.load(cons)

# define matricies  and end vectors to be used to calculate partition functions
begin = sp.Matrix([[0, 1]])
N = sp.Matrix([[(Kn * W), 1], [Kn, 1]])
R = sp.Matrix([[(Kr * W), 1], [Kr, 1]])
C = sp.Matrix([[(Kc * W), 1], [Kc, 1]])
end = sp.Matrix([[1], [1]])

# Build dictionaries of partition functions, partial derivs with respect
# to K, and fraction folded.

q_dict = {}
dqdKn_dict = {}
dqdKr_dict = {}
dqdKc_dict = {}
frac_folded_dict = {}

# Number of repeats of each type.  Seems like they should be floats, but
# I get an error in the matrix multiplication (q_dict) if they are declared to be.

for construct in constructs:

    # Make partition function dictionary and expressions for fraction folded.
    # Note, only one pf is generated per construct, even when there are multiple melts.

    matrixlist = construct.split("_")
    q_dict[construct + "_q"] = begin

    for i in range(0, len(matrixlist)):
        num_Ni = 0
        num_Ri = 0
        num_Ci = 0
        if matrixlist[i] == "N":
            num_Ni = 1
        if matrixlist[i] == "R":
            num_Ri = 1
        if matrixlist[i] == "C":
            num_Ci = 1

        q_dict[construct + "_q"] = (
            q_dict[construct + "_q"]
            * np.linalg.matrix_power(N, num_Ni)
            * np.linalg.matrix_power(R, num_Ri)
            * np.linalg.matrix_power(C, num_Ci)
        )

    q_dict[construct + "_q"] = q_dict[construct + "_q"] * end

    # Next two lines convert from sp.Matrix to np.array to something else.
    # Not sure the logic here, but it works.

    q_dict[construct + "_q"] = np.array(q_dict[construct + "_q"])
    q_dict[construct + "_q"] = q_dict[construct + "_q"].item(0)

    # Partial derivs wrt Kn dictionary.
    dqdKn_dict[construct + "_dqdKn"] = sp.diff(q_dict[construct + "_q"], Kn)

    # Partial derivs wrt Kr dictionary.
    dqdKr_dict[construct + "_dqdKr"] = sp.diff(q_dict[construct + "_q"], Kr)

    # Partial derivs wrt Kc dictionary.
    dqdKc_dict[construct + "_dqdKc"] = sp.diff(q_dict[construct + "_q"], Kc)

    # Fraction folded dictionary.
    frac_folded_dict[construct + "_frac_folded"] = (
        Kn / (q_dict[construct + "_q"]) * dqdKn_dict[construct + "_dqdKn"]
        + Kr / (q_dict[construct + "_q"]) * dqdKr_dict[construct + "_dqdKr"]
        + Kc / (q_dict[construct + "_q"]) * dqdKc_dict[construct + "_dqdKc"]
    ) / (len(matrixlist))

# The loop below replaces K's and W's the fraction folded terms in the
# dictionary with DGs, ms, and denaturant concentrations.  The simplify line
# is really important for making compact expressions for fraction folded.
# This simplification greatly speeds up fitting.  The last line
# converts from a sympy object to a string, to allow for json dump.

for construct in frac_folded_dict:
    frac_folded_dict[construct] = frac_folded_dict[construct].subs(
        {
            Kn: (np.exp(-((dGN - (mi * denat)) / RT))),
            Kr: (np.exp(-((dGR - (mi * denat)) / RT))),
            Kc: (np.exp(-((dGC - (mi * denat)) / RT))),
            W: (np.exp(-dGinter / RT)),
        }
    )
    frac_folded_dict[construct] = sp.simplify(frac_folded_dict[construct])
    frac_folded_dict[construct] = str(frac_folded_dict[construct])

with open("{0}{1}_frac_folded_dict.txt".format(path, proj_name), "wb") as f:
    json.dump(frac_folded_dict, f)


stop = time.time()
runtime = stop - start
print("\nThe elapsed time was " + str(runtime) + " sec")


print("\n")
print(
    """NOTE! You must quit the kernel before running the Ising fitter
      to clear the namespace. Otherwise you will get an error because
      numpy.exp has been reassigned."""
)
