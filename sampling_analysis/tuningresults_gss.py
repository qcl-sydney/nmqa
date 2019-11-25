import sys
sys.path.append('../paduaq')
from pdpoints import dims_padua_set

SIMULATIONSDICT = {}
GRIDSIZE=25 # number of data qubits if padua qubits are used. else, total number of regularly arranged sensing qubits
MULTIPLIER=5

COARSEGRID=16
FINEGRID=81
REG4 = 4
REG9 = 9
REG36= 36
REMOVE_DUPLICATES=25

for padua_order in ["no_padua", "reg4", "reg9", "regcoarse", "reg36", "regfine", 1, 2, 3, 4, 5, 10] : # Padua order for cheb2chev function
    
    SIMULATIONSDICT[padua_order] = {}
    
    if padua_order == "no_padua":
        SIMULATIONSDICT[padua_order]["max_iterations"] = GRIDSIZE * MULTIPLIER
        SIMULATIONSDICT[padua_order]["num_of_nodes"] = GRIDSIZE
    
    if padua_order == "reg4":
        SIMULATIONSDICT[padua_order]["max_iterations"] = REG4 * MULTIPLIER
        SIMULATIONSDICT[padua_order]["num_of_nodes"] = GRIDSIZE + COARSEGRID

    if padua_order == "reg9":
        SIMULATIONSDICT[padua_order]["max_iterations"] = REG9 * MULTIPLIER
        SIMULATIONSDICT[padua_order]["num_of_nodes"] = GRIDSIZE + COARSEGRID
            
    if padua_order == "regcoarse":
        SIMULATIONSDICT[padua_order]["max_iterations"] = COARSEGRID * MULTIPLIER
        SIMULATIONSDICT[padua_order]["num_of_nodes"] = GRIDSIZE + COARSEGRID
    
    if padua_order == "reg36":
        SIMULATIONSDICT[padua_order]["max_iterations"] = REG36 * MULTIPLIER
        SIMULATIONSDICT[padua_order]["num_of_nodes"] = GRIDSIZE + COARSEGRID

    if padua_order == "regfine":
        SIMULATIONSDICT[padua_order]["max_iterations"] = (FINEGRID - REMOVE_DUPLICATES) * MULTIPLIER # RED0
        SIMULATIONSDICT[padua_order]["num_of_nodes"] = FINEGRID # 56 unique sensor qubits + 25 data qubits; after overlapping sensor qubits are removed

    if not isinstance(padua_order, str):
        SIMULATIONSDICT[padua_order]["max_iterations"] = dims_padua_set(padua_order) * MULTIPLIER
        SIMULATIONSDICT[padua_order]["num_of_nodes"] = dims_padua_set(padua_order) + GRIDSIZE
    
    SIMULATIONSDICT[padua_order]["linear"] = False
    SIMULATIONSDICT[padua_order]["repts"] = 50
    SIMULATIONSDICT[padua_order]["functype"]= 'cheb2fun'
    
    for idx_expandtype in ["Uniform", "TruncGauss"]:
    
        SIMULATIONSDICT[padua_order][idx_expandtype] = {}
        SIMULATIONSDICT[padua_order][idx_expandtype]["optimal"] = {}
        SIMULATIONSDICT[padua_order][idx_expandtype]["optimal"]["idx_1"] = None
        SIMULATIONSDICT[padua_order][idx_expandtype]["optimal"]["idx_2"] = None   
        SIMULATIONSDICT[padua_order][idx_expandtype]["optimal"]["idx_3"] = None  
    
        SIMULATIONSDICT[padua_order][idx_expandtype]["zerolambda"] = {}
        SIMULATIONSDICT[padua_order][idx_expandtype]["zerolambda"]["idx_1"] = None
        SIMULATIONSDICT[padua_order][idx_expandtype]["zerolambda"]["idx_2"] = 0
        SIMULATIONSDICT[padua_order][idx_expandtype]["zerolambda"]["idx_3"] = None   

    SIMULATIONSDICT[padua_order]["GRAPHMIN"] = None # Figure plotting
    SIMULATIONSDICT[padua_order]["GRAPHMAX"] = None #1 # Figure plotting
    SIMULATIONSDICT[padua_order]["R_MAX"] = None # Figure plotting
    SIMULATIONSDICT[padua_order]["F_MAX"] = None # Figure plotting
    


###########################################################
# NUMERICAL OPTIMISATION RESULTS
###########################################################




