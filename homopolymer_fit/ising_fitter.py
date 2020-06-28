"""
Fits an Ising model to NRC capped homopolymer repeats, plots data and fits,
performs boostraps, calculates boostrap statistics and absolute parameter
correlations, and plots histograms and parameter correlation plots.
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import lmfit, json
import csv
import time
import os
import sys

start = time.time()

PATH = os.path.dirname(os.path.abspath(__file__))

proj_name = "cANK"

plt.close()
plt.clf

RT = 0.001987 * 298.15  #  R in kcal/mol/K, T in Kelvin.

#  Dictionary of frac folded eqns from partition function generator script.
with open(
    os.path.join(PATH, f"{proj_name}_frac_folded_dict.json"), "r"
) as ffd:
    frac_folded_dict = json.load(ffd)

with open(
    os.path.join(PATH, f"{proj_name}_constructs.json"), "r"
) as construct:
    constructs = json.load(construct)

with open(os.path.join(PATH, f"{proj_name}_melts.json"), "r") as m:
    melts = json.load(m)

num_melts = len(melts)
num_constructs = len(constructs)

melt_data_dict = {}
for melt in melts:
    melt_data_dict[melt] = np.load(os.path.join(PATH, f"{melt}.npy"))

# Compile fraction folded expressions.
comp_frac_folded_dict = {}
for construct in constructs:
    frac_folded_string = frac_folded_dict[construct + "_frac_folded"]
    comp_frac_folded = compile(
        frac_folded_string, "{}_comp_ff".format(construct), "eval"
    )
    comp_frac_folded_dict[construct + "_comp_ff"] = comp_frac_folded

# CREATE INITIAL GUESSES
# First, thermodynamic parameters.  These are Global.
init_guesses = lmfit.Parameters()
init_guesses.add("dGN", value=6)
init_guesses.add("dGR", value=5)
init_guesses.add("dGC", value=6)
init_guesses.add("dGinter", value=-12)
init_guesses.add("mi", value=-1.0)

# Next, baseline parameters.  These are local.
for melt in melts:
    init_guesses.add("af_{}".format(melt), value=0.02)
    init_guesses.add("bf_{}".format(melt), value=1)
    init_guesses.add("au_{}".format(melt), value=0.0)
    init_guesses.add("bu_{}".format(melt), value=0.0)

# Transfers init_guesses to params for fitting, but init_guesses are maintained.
params = init_guesses


def fitting_function(params, denat, frac_folded, melt):
    af = params["af_{}".format(melt)].value
    bf = params["bf_{}".format(melt)].value
    au = params["au_{}".format(melt)].value
    bu = params["bu_{}".format(melt)].value
    dGN = params["dGN"].value
    dGR = params["dGR"].value
    dGC = params["dGC"].value
    dGinter = params["dGinter"].value
    mi = params["mi"].value
    return ((af * denat) + bf) * frac_folded + (
        ((au * denat) + bu) * (1 - frac_folded)
    )


# Objective function creates an array of residuals to be used by lmfit minimize.
def objective(params):
    resid_dict = {}
    dGN = params["dGN"].value
    dGR = params["dGR"].value
    dGC = params["dGC"].value
    dGinter = params["dGinter"].value
    mi = params["mi"].value
    for melt in melts:
        denat = melt_data_dict[melt][:, 0]  # A numpy array of type str
        norm_sig = melt_data_dict[melt][:, 1]  # A numpy array of type str
        denat = denat.astype(float)  # A numpy array of type float
        norm_sig = norm_sig.astype(float)  # A numpy array of type float
        string_to_eval = comp_frac_folded_dict[melt[:-2] + "_comp_ff"]
        frac_folded = eval(string_to_eval)
        af = params["af_{}".format(melt)].value
        bf = params["bf_{}".format(melt)].value
        au = params["au_{}".format(melt)].value
        bu = params["bu_{}".format(melt)].value
        resid = norm_sig - fitting_function(params, denat, frac_folded, melt)
        resid_dict[melt + "_resid"] = resid
    residuals = np.concatenate(list(resid_dict.values()))
    return residuals


# Fit with lmfit
result = lmfit.minimize(objective, init_guesses)
fit_resid = result.residual

# Print out features of the data, the fit, and optimized param values
print(("There are a total of {} data sets.".format(num_melts)))
print(("There are {} observations.".format(result.ndata)))
print(("There are {} fitted parameters.".format(result.nvarys)))
print(("There are {} degrees of freedom. \n".format(result.nfree)))
print(
    ("The sum of squared residuals (SSR) is: {0:7.4f}".format(result.chisqr))
)
print(("The reduced SSR (SSR/DOF): {0:8.6f} \n".format(result.redchi)))

dGN = result.params["dGN"].value
dGR = result.params["dGR"].value
dGC = result.params["dGC"].value
dGinter = result.params["dGinter"].value
mi = result.params["mi"].value

print("Optimized parameter values:")
print(("dGN = {0:8.4f}".format(result.params["dGN"].value)))
print(("dGR = {0:8.4f}".format(result.params["dGR"].value)))
print(("dGC ={0:8.4f}".format(result.params["dGC"].value)))
print(("dGinter ={0:8.4f}".format(result.params["dGinter"].value)))
print(("mi ={0:8.4f}".format(result.params["mi"].value)))

print("\nWriting best fit parameter and baseline files")

# Compile a list of optimized Ising params and write to file.
fitted_ising_params = [
    ["dGN", result.params["dGN"].value],
    ["dGR", result.params["dGR"].value],
    ["dGC", result.params["dGC"].value],
    ["dGinter", result.params["dGinter"].value],
    ["mi", result.params["mi"].value],
    ["Chi**2", result.chisqr],
    ["RedChi", result.redchi],
]

with open(
    os.path.join(PATH, f"{proj_name}_fitted_Ising_params.csv"), "w"
) as n:
    writer = csv.writer(n, delimiter=",")
    writer.writerows(fitted_ising_params)
n.close()

# Compile a list of optimized baseline params and write to file.
fitted_base_params = []
for melt in melts:
    af = result.params["af_%s" % (melt)].value
    bf = result.params["bf_%s" % (melt)].value
    au = result.params["au_%s" % (melt)].value
    bu = result.params["bu_%s" % (melt)].value
    fitted_base_params.append([melt, af, bf, au, bu])
with open(
    os.path.join(PATH, f"{proj_name}_fitted_baseline_params.csv"), "w"
) as m:
    writer = csv.writer(m, delimiter=",")
    writer.writerows(fitted_base_params)
m.close()

stop = time.time()
runtime = stop - start
print(("\nThe elapsed time was " + str(runtime) + " sec"))

print("\nPlotting results...\n")

# The function "baseline_adj" gives an adjusted y value based on fitted baseline
# parameters (fraction folded).
def baseline_adj(y, x, params, construct):
    af = result.params["af_{}".format(construct)].value
    bf = result.params["bf_{}".format(construct)].value
    au = result.params["au_{}".format(construct)].value
    bu = result.params["bu_{}".format(construct)].value
    return (y - (bu + (au * x))) / ((bf + (af * x)) - (bu + (au * x)))


# Defining global best-fit parameters
dGN = result.params["dGN"].value
dGR = result.params["dGR"].value
dGC = result.params["dGC"].value
dGinter = result.params["dGinter"].value
mi = result.params["mi"].value

# The function fit_model used for plotting best-fit lines and for adding
# residuals to best-fit lines in bootstrapping.  Normalized, not frac folded.
def fit_model(params, x, melt):
    denat = x
    af = result.params["af_{}".format(melt)].value
    bf = result.params["bf_{}".format(melt)].value
    au = result.params["au_{}".format(melt)].value
    bu = result.params["bu_{}".format(melt)].value
    dGN = params["dGN"].value
    dGR = params["dGR"].value
    dGC = params["dGC"].value
    dGinter = params["dGinter"].value
    mi = params["mi"].value
    frac_folded = eval(comp_frac_folded_dict[melt[:-2] + "_comp_ff"])
    return ((af * denat) + bf) * frac_folded + (
        ((au * denat) + bu) * (1 - frac_folded)
    )


# Finding the maximum denaturant value out of all the melts to
# set x axis bound
denat_maxer = np.zeros(0)
for melt in melts:
    denat_maxer = np.concatenate((denat_maxer, melt_data_dict[melt][:, 0]))
denat_maxer_list = denat_maxer.tolist()
denat_max = float(max(denat_maxer_list))
denat_bound = np.around(denat_max, 1) + 0.2

# Denaturant values to use when evaluating fits.  Determines how smooth the
# fitted curve will be, based on the third value (300) in the argument below.
# I might keep using this for fraction_foldeed, but for nomralized baseline
# use a local set of points for each melt, so as not to extrapolate the
# bselines too far.
denat_fit = np.linspace(0, denat_bound, 300)

# defining a dictionary using the first melt of each construct (construct_1)
# Move this to the plotting part, and why not do this for all constructs?
construct1_data_dict = {}
for construct in constructs:
    construct1_data_dict[construct] = np.load(
        os.path.join(PATH, f"{construct}_1.npy")
    )

# The four dictionaries below define lower and upper denaturant limnits to be
# used for plotting normalized curves, so crazy-long baseline extrapolations
# are not shown.  Do both for melts and construct 1.   These are then used
# to create 300-point synthetic baselines in the fifth and sixth dictionaries.
melt_lower_denat_dict = {}
for melt in melts:
    melt_lower_denat_dict[melt] = (
        round(float(min(melt_data_dict[melt][:, 0]))) - 0.2
    )

melt_upper_denat_dict = {}
for melt in melts:
    melt_upper_denat_dict[melt] = (
        round(float(max(melt_data_dict[melt][:, 0]))) + 0.2
    )

construct1_lower_denat_dict = {}
for construct in constructs:
    construct1_lower_denat_dict[construct] = (
        round(float(min(construct1_data_dict[construct][:, 0]))) - 0.2
    )

construct1_upper_denat_dict = {}
for construct in constructs:
    construct1_upper_denat_dict[construct] = (
        round(float(max(construct1_data_dict[construct][:, 0]))) + 0.2
    )

melt_denat_synthetic_dict = {}
for melt in melts:
    melt_denat_synthetic_dict[melt] = np.linspace(
        melt_lower_denat_dict[melt], melt_upper_denat_dict[melt], 300
    )

construct1_denat_synthetic_dict = {}
for construct in constructs:
    construct1_denat_synthetic_dict[construct] = np.linspace(
        construct1_lower_denat_dict[construct],
        construct1_upper_denat_dict[construct],
        300,
    )

""" Global Plot Aesthetics"""
# Defining how the plots are colored
num_melt_colors = num_melts
num_construct_colors = num_constructs
coloration = plt.get_cmap("hsv")


# Defining title font
title_font = {
    "family": "arial",
    "color": "black",
    "weight": "normal",
    "size": 16,
}

# Defining label font
label_font = {
    "family": "arial",
    "color": "black",
    "weight": "normal",
    "size": 14,
}

"""First Plot: Fraction Folded by Melt"""
# extracting the melt data and creating plot lines for each melt
colorset = 0  # counter to control color of curves and points
for melt in melts:
    colorset = colorset + 1
    denat = melt_data_dict[melt][:, 0]  # A numpy array of type str
    norm_sig = melt_data_dict[melt][:, 1]  # A numpy array of type str
    denat = denat.astype(float)  # A numpy array of type float
    norm_sig = norm_sig.astype(float)  # A numpy array of type float
    y_adj = baseline_adj(norm_sig, denat, result.params, melt)
    y_fit = fit_model(result.params, denat_fit, melt)
    y_fit_adj = baseline_adj(y_fit, denat_fit, result.params, melt)
    plt.plot(
        denat,
        y_adj,
        "o",
        color=coloration(colorset / num_melt_colors),
        label=melt[:-2] + " melt " + melt[-1],
    )
    plt.plot(
        denat_fit, y_fit_adj, "-", color=coloration(colorset / num_melt_colors)
    )

# set axis limits
axes = plt.gca()
axes.set_xlim([-0.1, denat_bound])
axes.set_ylim([-0.1, 1.1])
axes.set_aspect(5.5)

# lot aesthetics and labels
plt.legend(loc="center", bbox_to_anchor=(1.25, 0.5), fontsize=8)
plt.title("Fraction Folded by Melt", fontdict=title_font)
plt.xlabel("Denaturant (Molar)", fontdict=label_font)
plt.ylabel("Fraction Folded", fontdict=label_font)

# saving plot in individual doc
plt.savefig(
    os.path.join(PATH, f"{proj_name}_plot_frac_folded_by_melt.png"),
    dpi=500,
    bbox_inches="tight",
)

# show plot in iPython window and then close
plt.show()
plt.close()
plt.clf

"""Second Plot: Normalized Signal by Melt"""
colorset = 0
for melt in melts:
    colorset = colorset + 1
    denat = melt_data_dict[melt][:, 0]  # A numpy array of type str
    norm_sig = melt_data_dict[melt][:, 1]  # A numpy array of type str
    denat = denat.astype(float)  # A numpy array of type float
    norm_sig = norm_sig.astype(float)  # A numpy array of type float
    y_fit = fit_model(result.params, melt_denat_synthetic_dict[melt], melt)
    plt.plot(
        denat,
        norm_sig,
        "o",
        color=coloration(colorset / num_melt_colors),
        label=melt[:-2] + " melt " + melt[-1],
    )
    plt.plot(
        melt_denat_synthetic_dict[melt],
        y_fit,
        "-",
        color=coloration(colorset / num_melt_colors),
    )

# set axis limits
axes = plt.gca()
axes.set_xlim([-0.1, denat_bound])
axes.set_ylim([-0.1, 1.1])
axes.set_aspect(5.5)

# plot aesthetics and labels
plt.legend(loc="center", bbox_to_anchor=(1.25, 0.5), fontsize=8)
plt.title("Normalized Signal by Melt", fontdict=title_font)
plt.xlabel("Denaturant (Molar)", fontdict=label_font)
plt.ylabel("Normalized Signal", fontdict=label_font)

# saving plot in individual doc
plt.savefig(
    os.path.join(PATH, f"{proj_name}_plot_normalized_by_melt.png"),
    dpi=500,
    bbox_inches="tight",
)

# show plot in iPython window and then close
plt.show()
plt.close()
plt.clf

"""Third Plot: Fraction Folded by Construct"""
colorset = 0
for construct in constructs:
    colorset = colorset + 1
    denat = construct1_data_dict[construct][:, 0]  # A numpy array of type str
    denat_line = construct1_data_dict[construct][
        :, 0
    ]  # A numpy array of type str
    norm_sig = construct1_data_dict[construct][
        :, 1
    ]  # A numpy array of type str
    denat = denat.astype(float)  # A numpy array of type float
    denat_line = denat_line.astype(float)  # A numpy array of type float
    norm_sig = norm_sig.astype(float)  # A numpy array of type float
    y_adj = baseline_adj(norm_sig, denat_line, result.params, construct + "_1")
    y_fit = fit_model(result.params, denat_fit, construct + "_1")
    y_fit_adj = baseline_adj(y_fit, denat_fit, result.params, construct + "_1")
    plt.plot(
        denat,
        y_adj,
        "o",
        color=coloration(colorset / num_construct_colors),
        label=construct,
    )
    plt.plot(
        denat_fit,
        y_fit_adj,
        "-",
        color=coloration(colorset / num_construct_colors),
    )

# set axis limits
axes = plt.gca()
axes.set_xlim([-0.1, denat_bound])
axes.set_ylim([-0.1, 1.1])
axes.set_aspect(5.5)

# plot aesthetics and labels
plt.legend(loc="center", bbox_to_anchor=(1.15, 0.5), fontsize=8)
plt.title("Fraction Folded by Construct", fontdict=title_font)
plt.xlabel("Denaturant (Molar)", fontdict=label_font)
plt.ylabel("Fraction Folded", fontdict=label_font)

# saving plot in individual doc
plt.savefig(
    os.path.join(PATH, f"{proj_name}_plot_frac_folded_by_construct.png"),
    dpi=500,
    bbox_inches="tight",
)

# show plot in iPython window and then close
plt.show()
plt.close()
plt.clf

"""Fourth Plot: Normalized Signal by Construct"""
colorset = 0
for construct in constructs:
    colorset = colorset + 1
    denat = construct1_data_dict[construct][:, 0]  # A numpy array of type str
    norm_sig = construct1_data_dict[construct][
        :, 1
    ]  # A numpy array of type str
    denat = denat.astype(float)  # A numpy array of type float
    norm_sig = norm_sig.astype(float)  # A numpy array of type float
    y_fit = fit_model(
        result.params,
        construct1_denat_synthetic_dict[construct],
        construct + "_1",
    )
    plt.plot(
        denat,
        norm_sig,
        "o",
        color=coloration(colorset / num_construct_colors),
        label=construct,
    )
    plt.plot(
        construct1_denat_synthetic_dict[construct],
        y_fit,
        "-",
        color=coloration(colorset / num_construct_colors),
    )

# set axis limits
axes = plt.gca()
axes.set_xlim([-0.1, denat_bound])
axes.set_ylim([-0.1, 1.1])
axes.set_aspect(5.5)

# plot aesthetics and labels
plt.legend(loc="center", bbox_to_anchor=(1.15, 0.5), fontsize=8)
plt.title("Normalized Signal by Construct", fontdict=title_font)
plt.xlabel("Denaturant (Molar)", fontdict=label_font)
plt.ylabel("Normalized Signal", fontdict=label_font)

# saving plot in individual doc
plt.savefig(
    os.path.join(PATH, f"{proj_name}_plot_normalized_by_construct.png"),
    dpi=500,
    bbox_inches="tight",
)

# show plot in iPython window and then close
plt.show()
plt.close()
plt.clf


"""BOOTSTRAP ANALYSIS"""
# Create list to store bootstrap iterations of values and define column titles
bs_param_values = []
bs_param_values.append(
    [
        "Bootstrap Iter",
        "dGN",
        "dGR",
        "dGC",
        "dGinter",
        "mi",
        "redchi**2",
        "bestchi**2",
    ]
)
# total number of bootstrap iterations
bs_iter_tot = eval(input("How many bootstrap iterations? "))

if bs_iter_tot == 0:
    sys.exit()

# bs_iter_tot = 10
bs_iter_count = 0  # Iteration counter
fit_resid_index = len(fit_resid) - 1

y_fitted_dict = {}
# Dictionary of 'true' normalized y values from fit at each denaturant value.
for melt in melts:
    denat = melt_data_dict[melt][:, 0]  # A numpy array of type str
    denat = denat.astype(float)  # A numpy array of type float
    y_fitted_dict[melt] = np.array(fit_model(result.params, denat, melt))

# Arrays to store bs fitted param values
dGN_vals = []
dGR_vals = []
dGC_vals = []
dGinter_vals = []
mi_vals = []

# Add residuals chosen at random (with replacement) to expected
# y values. Note-residuals are combined ACROSS melts.
for j in range(0, bs_iter_tot):
    rand_resid_dict = (
        {}
    )  # Clears the random data for each bootsterap iteration
    bs_iter_count = bs_iter_count + 1
    print(
        "Bootstrap iteration {0} out of {1}".format(bs_iter_count, bs_iter_tot)
    )

    for melt in melts:
        rand_resid = []
        denat = melt_data_dict[melt][:, 0]  # A numpy array of type str
        denat = denat.astype(float)  # A numpy array of type float

        for x in range(0, len(denat)):  # Creastes a list of random residuals
            rand_int = np.random.randint(0, fit_resid_index)
            rand_resid.append(fit_resid[rand_int])

        rand_resid_dict[melt] = np.array(rand_resid)
        y_bootstrap = y_fitted_dict[melt] + rand_resid_dict[melt]
        z_max, z_min = y_bootstrap.max(), y_bootstrap.min()
        melt_data_dict[melt][:, 1] = (y_bootstrap - z_min) / (z_max - z_min)

    bs_result = lmfit.minimize(objective, init_guesses)
    bs_chisqr = bs_result.chisqr
    bs_red_chisqr = bs_result.redchi

    dGN = bs_result.params["dGN"].value
    dGR = bs_result.params["dGR"].value
    dGC = bs_result.params["dGC"].value
    dGinter = bs_result.params["dGinter"].value
    mi = bs_result.params["mi"].value

    # Store each value in a list for plotting
    dGN_vals.append(dGN)
    dGR_vals.append(dGR)
    dGC_vals.append(dGC)
    dGinter_vals.append(dGinter)
    mi_vals.append(mi)

    # Append bootstrapped global parameter values for ouput to a file
    bs_param_values.append(
        [bs_iter_count, dGN, dGR, dGC, dGinter, mi, bs_red_chisqr, bs_chisqr]
    )

with open(os.path.join(PATH, f"{proj_name}_bootstrap_params.csv"), "w") as n:
    writer = csv.writer(n, delimiter=",")
    writer.writerows(bs_param_values)
n.close()

## BOOTSTRAP STATISTICS
bs_param_values_fullarray = np.array(bs_param_values)
bs_param_values_array = bs_param_values_fullarray[1:, 1:-2].astype(
    np.float
)  # End at -2 since last two columns
# are chi square statistics

bs_param_names = bs_param_values_fullarray[0][1:-2]

statistics = [
    "mean",
    "median",
    "stdev",
    "2.5% CI",
    "16.6% CI",
    "83.7% CI",
    "97.5% CI",
]

bs_statistics_df = pd.DataFrame(columns=statistics)

i = 0
for param in bs_param_names:
    bs_statistics = []
    bs_statistics.append(np.mean(bs_param_values_array[:, i]))
    bs_statistics.append(np.median(bs_param_values_array[:, i]))
    bs_statistics.append(np.std(bs_param_values_array[:, i]))
    bs_statistics.append(np.percentile(bs_param_values_array[:, i], 2.5))
    bs_statistics.append(np.percentile(bs_param_values_array[:, i], 16.7))
    bs_statistics.append(np.percentile(bs_param_values_array[:, i], 83.3))
    bs_statistics.append(np.percentile(bs_param_values_array[:, i], 97.5))
    bs_statistics_df.loc[param] = bs_statistics
    i = i + 1

bs_statistics_df.to_csv(os.path.join(PATH, f"{proj_name}_bootstrap_stats.csv"))

# Compute and write out correlation coefficient matrix
corr_coef_matrix = np.corrcoef(bs_param_values_array, rowvar=False)
corr_coef_df = pd.DataFrame(
    corr_coef_matrix, columns=bs_param_names, index=bs_param_names
)
corr_coef_df.to_csv(
    os.path.join(PATH, f"{proj_name}_bootstrap_corr_coefs.csv")
)


########### BOOTSTRAP PLOTTING

# Specify the names of parameters to be compared to see correlation.
corr_params = ["dGN", "dGR", "dGC", "dGinter", "mi"]

# These are a second set of parameter names that follow in the same order
# as in corr_params.  They are formatted using TeX-style names so that Deltas
# and subscripts will be plotted.  The would not be good key names for dictionaries
corr_param_labels = [
    "$\Delta$G$_N$",
    "$\Delta$G$_R$",
    "$\Delta$G$_C$",
    "$\Delta$G$_{i, i-1}$",
    "m$_i$",
]

num_corr_params = len(corr_params)
gridsize = num_corr_params  # Determines the size of the plot grid.

# Dictionary of fitted parameter values.
corr_params_dict = {
    "dGN": dGN_vals,
    "dGR": dGR_vals,
    "dGC": dGC_vals,
    "dGinter": dGinter_vals,
    "mi": mi_vals,
}

# PDF that stores a grid of the correlation plots
with PdfPages(os.path.join(PATH, f"{proj_name}_Corr_Plots.pdf")) as pdf:
    fig, axs = plt.subplots(ncols=gridsize, nrows=gridsize, figsize=(12, 12))

    # Turns off axes on lower triangle
    axs[1, 0].axis("off")
    axs[2, 0].axis("off")
    axs[2, 1].axis("off")
    axs[3, 0].axis("off")
    axs[3, 1].axis("off")
    axs[3, 2].axis("off")
    axs[4, 0].axis("off")
    axs[4, 1].axis("off")
    axs[4, 2].axis("off")
    axs[4, 3].axis("off")

    # Defines the position of the y paramater from the array of params
    hist_param_counter = 0
    while hist_param_counter < num_corr_params:
        hist_param_label = corr_param_labels[hist_param_counter]
        hist_param = corr_params[hist_param_counter]
        # Start fixing labels here
        # plt.xticks(fontsize=8)
        # axs[hist_param_counter, hist_param_counter].tick_params(fontsize=8)
        # axs[hist_param_counter, hist_param_counter].yticks(fontsize=8)
        axs[hist_param_counter, hist_param_counter].hist(
            corr_params_dict[hist_param]
        )
        axs[hist_param_counter, hist_param_counter].set_xlabel(
            hist_param_label, fontsize=14, labelpad=5
        )
        hist_param_counter = hist_param_counter + 1

    # This part generates the correlation plots
    y_param_counter = 0
    while y_param_counter < num_corr_params - 1:
        # Pulls the parameter name for the y-axis label (with TeX formatting)
        yparam_label = corr_param_labels[y_param_counter]
        # Pulls the parameter name to be plotted on the y-axis
        yparam = corr_params[y_param_counter]

        # Defines the position of the x paramater from the array of params.
        # The + 1 offest avoids correlating a parameter with itself.
        x_param_counter = y_param_counter + 1

        while x_param_counter < num_corr_params:
            # pulls the parameter name for the x-axis label (with TeX formatting)
            xparam_label = corr_param_labels[x_param_counter]
            # Pulls the parameter name to be plotted on the x-axis
            xparam = corr_params[x_param_counter]

            x_vals = corr_params_dict[xparam]
            y_vals = corr_params_dict[yparam]

            # plt.xticks(fontsize=8)
            # plt.yticks(fontsize=8)
            # plotting scatters with axes.  +1 shifts a plot to the right from main diagonal
            axs[y_param_counter, x_param_counter].plot(x_vals, y_vals, ".")

            # The if statement below turns off numbers on axes if not the right column and
            # not the main diagonal.
            if x_param_counter < num_corr_params - 1:
                axs[y_param_counter, x_param_counter].set_xticklabels([])
                axs[y_param_counter, x_param_counter].set_yticklabels([])

            if y_param_counter == 0:  # Puts labels above axes on top row
                axs[y_param_counter, x_param_counter].xaxis.set_label_position(
                    "top"
                )
                axs[y_param_counter, x_param_counter].set_xlabel(
                    xparam_label, labelpad=10, fontsize=14
                )
                axs[y_param_counter, x_param_counter].xaxis.tick_top()
                if (
                    x_param_counter < num_corr_params - 1
                ):  # Avoids eliminating y-scale from upper right corner
                    axs[y_param_counter, x_param_counter].set_yticklabels([])

            if (
                x_param_counter == num_corr_params - 1
            ):  #  Puts labels right of right column
                axs[y_param_counter, x_param_counter].yaxis.set_label_position(
                    "right"
                )
                axs[y_param_counter, x_param_counter].set_ylabel(
                    yparam_label, rotation=0, labelpad=30, fontsize=14
                )
                axs[y_param_counter, x_param_counter].set_xticklabels([])
                axs[y_param_counter, x_param_counter].yaxis.tick_right()

            # Determin correlation coefficient and display under subplot title
            # Note, there is no code that displays this value at the moment.
            # corr_coef = np.around(np.corrcoef(x_vals, y_vals), 3)

            # min and max values of the x param
            x_min = min(x_vals)
            x_max = max(x_vals)

            # fitting a straight line to the correlation scatterplot
            fit_array = np.polyfit(x_vals, y_vals, 1)
            fit_deg1_coef = fit_array[0]
            fit_deg0_coef = fit_array[1]
            fit_x_vals = np.linspace(x_min, x_max, 10)
            fit_y_vals = fit_deg1_coef * fit_x_vals + fit_deg0_coef

            # plotting correlation line fits
            axs[y_param_counter, x_param_counter].plot(fit_x_vals, fit_y_vals)
            plt.subplots_adjust(wspace=0, hspace=0)

            x_param_counter = x_param_counter + 1
        y_param_counter = y_param_counter + 1

    pdf.savefig(bbox_inches="tight")
