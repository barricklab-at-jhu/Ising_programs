'''
Generates a partition function and fraction folded expressions for an NRXC
capped heteropolymer repeat array.

Most recent revision: 05/21/2020, Doug Barrick

'''

from __future__ import division
import sympy as sp
import numpy as np
import json
import time

path = '/Users/dougbarrick/OneDrive - Johns Hopkins/Manuscripts/Ising_program/\
Scripts/scripts_mutations/T4V_two_mi_values_2020_05_10/T4V_mR_mX_spyder/'

proj_name = 'T4V_NRC_2mi'

start = time.time()

print('\nGenerating partition functions and fraction folded expressions.  This may take a minute...')

# Parameters for partition function calculation.  Note these are sympy symbols.
RT = sp.Symbol('RT')
dGN = sp.Symbol('dGN')
dGR = sp.Symbol('dGR')
dGX = sp.Symbol('dGX')
dGC = sp.Symbol('dGC')
dGRR = sp.Symbol('dGRR')
dGXX = sp.Symbol('dGXX')
dGXR = sp.Symbol('dGXR')
dGRX = sp.Symbol('dGRX')
mR = sp.Symbol('mR')
mX = sp.Symbol('mX')
denat = sp.Symbol('denat')
KN = sp.Symbol('KN')
KR = sp.Symbol('KR')
KX = sp.Symbol('KX')
KC = sp.Symbol('KC')
TRR = sp.Symbol('TRR')
TXX = sp.Symbol('TXX')
TXR = sp.Symbol('TXR')
TRX = sp.Symbol('TRX')

np.exp = sp.Function('np.exp')

with open('{0}{1}_constructs.txt'.format(path, proj_name)) as cons:
    constructs = json.load(cons)

#define weight matricies and end vectors to be used to calculate partition functions
begin = sp.Matrix([[0,1]])
woN = sp.Matrix([[KN,1],[KN,1]]) # Leave off coupling term.  No zeroth repeat
woR = sp.Matrix([[KR,1],[KR,1]]) # Leave off coupling term.  No zeroth repeat
woX = sp.Matrix([[KX,1],[KX,1]]) # Leave off coupling term.  No zeroth repeat

wnR = sp.Matrix([[KR*TRR,1],[KR,1]]) # Treat coupling same as an rR interface, as usual
wrR = sp.Matrix([[KR*TRR,1],[KR,1]]) # Same as wnR above
wxR = sp.Matrix([[KR*TXR,1],[KR,1]])

wnX = sp.Matrix([[KX*TRX,1],[KX,1]]) # Different coupling than 
wrX = sp.Matrix([[KX*TRX,1],[KX,1]]) # Different coupling than wnX above
wxX = sp.Matrix([[KX*TXX,1],[KX,1]])

wrC = sp.Matrix([[KC*TRR,1],[KC,1]]) # Treat copuling same as an rR interface, as usual
wxC = sp.Matrix([[KC*TXR,1],[KC,1]])

end = sp.Matrix([[1],[1]])

# Build a dictionary of these matrices
w_dict = {'oN': woN, 
         'oR': woR,
         'oX': woX,
         'nR': wnR,
         'rR': wrR,
         'xR': wxR,
         'nX': wnX,
         'rX': wrX,
         'xX': wxX,
         'rC': wrC,
         'xC': wxC}

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
    
    repeat_list = construct.split('_')
    repeat_list.insert(0, 'o')
    
    pairs_list = []
    i = 1
    while i < len(repeat_list):
        pair = repeat_list[i-1].lower() + repeat_list[i] # Need to convert prev rept to lower case
        pairs_list.append(pair)
        i = i + 1
    
    q_dict[construct + '_q'] = begin
    
    for pair in pairs_list:
        q_dict[construct + '_q'] = q_dict[construct + '_q'] * w_dict[pair]
    
    q_dict[construct + '_q']  = q_dict[construct + '_q'] * end
    
    # Next two lines convert from sp.Matrix to np.array to something else.
    # Not sure the logic here, but it works.
    q_dict[construct + '_q'] = np.array(q_dict[construct + '_q']) 
    q_dict[construct + '_q'] = q_dict[construct + '_q'].item(0)
        
    # Partial derivs wrt KN dictionary.
    dqdKN_dict[construct + '_dqdKN'] \
        = sp.diff(q_dict[construct + '_q'], KN)

    # Partial derivs wrt KR dictionary.
    dqdKR_dict[construct + '_dqdKR'] \
        = sp.diff(q_dict[construct + '_q'], KR)
    
    # Partial derivs wrt KX dictionary.
    dqdKX_dict[construct + '_dqdKX'] \
        = sp.diff(q_dict[construct + '_q'], KX)

    # Partial derivs wrt KC dictionary.
    dqdKC_dict[construct + '_dqdKC'] \
        = sp.diff(q_dict[construct + '_q'], KC)
    
    # Fraction folded dictionary.  
    frac_folded_dict[construct + '_frac_folded'] \
        = (KN/(q_dict[construct + '_q']) * dqdKN_dict[construct + '_dqdKN'] \
         + KR/(q_dict[construct + '_q']) * dqdKR_dict[construct + '_dqdKR'] \
         + KX/(q_dict[construct + '_q']) * dqdKX_dict[construct + '_dqdKX'] \
         + KC/(q_dict[construct + '_q']) * dqdKC_dict[construct + '_dqdKC']) \
        / (len(pairs_list))

# The loop below replaces K's and W's the fraction folded terms in the 
# dictionary with DGs, ms, and denaturant concentrations.  The simplify line
# is really important for making compact expressions for fraction folded.
# This simplification greatly speeds up fitting.  The last line
# converts from a sympy object to a string, to allow for json dump.

for construct in frac_folded_dict:
    frac_folded_dict[construct] = frac_folded_dict[construct].subs({
    KN:(np.exp(-((dGN + (mR*denat))/RT))), 
    KR:(np.exp(-((dGR + (mR*denat))/RT))), 
    KX:(np.exp(-((dGX + (mX*denat))/RT))),
    KC:(np.exp(-((dGC + (mR*denat))/RT))),  
    TRR:(np.exp(-dGRR/RT)), 
    TXX:(np.exp(-dGXX/RT)), 
    TXR:(np.exp(-dGXR/RT)),
    TRX:(np.exp(-dGRX/RT)) }) 
    frac_folded_dict[construct] = sp.simplify(frac_folded_dict[construct])
    frac_folded_dict[construct] = str(frac_folded_dict[construct])                                

with open('{0}{1}_frac_folded_dict.txt'.format(path, proj_name), 'wb') as f:
    json.dump(frac_folded_dict, f)
         
stop = time.time()
runtime = stop - start
print('\nThe elapsed time was ' + str(runtime) + ' sec') 

print('\n')
print('''NOTE! You must quit the kernel before running the Ising fitter
      to clear the namespace. Otherwise you will get an error because
      numpy.exp has been reassigned.''')
      
