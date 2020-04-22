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

for padua_order in ["no_padua", "reg4", "reg9", "regcoarse", "reg36", "regfine", 1, 2, 3, 4, 5, 10] : 
    
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

SIMULATIONSDICT["no_padua"]["Opt_Beta_Expn"] = "TruncGauss"
# Optimal trunc termination error: 6.99192672094 Params: (23, 18, 4) DONE
# Optimal uniform termination error: 7.45751852236 Params: (24, 13, 0)


SIMULATIONSDICT["no_padua"]["Uniform"]["optimal"]["idx_1"] = 24
SIMULATIONSDICT["no_padua"]["Uniform"]["optimal"]["idx_2"] = 13
SIMULATIONSDICT["no_padua"]["Uniform"]["optimal"]["idx_3"] = 0
SIMULATIONSDICT["no_padua"]["Uniform"]["zerolambda"]["idx_1"] = 24
SIMULATIONSDICT["no_padua"]["Uniform"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT["no_padua"]["Uniform"]["zerolambda"]["idx_3"] = None

SIMULATIONSDICT["no_padua"]["TruncGauss"]["optimal"]["idx_1"] = 23
SIMULATIONSDICT["no_padua"]["TruncGauss"]["optimal"]["idx_2"] = 18
SIMULATIONSDICT["no_padua"]["TruncGauss"]["optimal"]["idx_3"] = 4
SIMULATIONSDICT["no_padua"]["TruncGauss"]["zerolambda"]["idx_1"] = 24
SIMULATIONSDICT["no_padua"]["TruncGauss"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT["no_padua"]["TruncGauss"]["zerolambda"]["idx_3"] = None


SIMULATIONSDICT["reg4"]["Opt_Beta_Expn"] = "TruncGauss"
# Optimal trunc termination error: 19.0186320524 Params: (12, 7, 0) - DONE
# Optimal uniform termination error: 19.3677669125 Params: (4, 18, 0)


SIMULATIONSDICT["reg4"]["Uniform"]["optimal"]["idx_1"] = 4
SIMULATIONSDICT["reg4"]["Uniform"]["optimal"]["idx_2"] = 18
SIMULATIONSDICT["reg4"]["Uniform"]["optimal"]["idx_3"] = 0
SIMULATIONSDICT["reg4"]["Uniform"]["zerolambda"]["idx_1"] = 4
SIMULATIONSDICT["reg4"]["Uniform"]["zerolambda"]["idx_2"] = 0
SIMULATIONSDICT["reg4"]["Uniform"]["zerolambda"]["idx_3"] = None

SIMULATIONSDICT["reg4"]["TruncGauss"]["optimal"]["idx_1"] = 12
SIMULATIONSDICT["reg4"]["TruncGauss"]["optimal"]["idx_2"] = 7
SIMULATIONSDICT["reg4"]["TruncGauss"]["optimal"]["idx_3"] = 0
SIMULATIONSDICT["reg4"]["TruncGauss"]["zerolambda"]["idx_1"] = 12
SIMULATIONSDICT["reg4"]["TruncGauss"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT["reg4"]["TruncGauss"]["zerolambda"]["idx_3"] = None

SIMULATIONSDICT["reg9"]["Opt_Beta_Expn"] = "TruncGauss"
# Optimal trunc termination error: 15.6759202334 Params: (28, 23, 0) - DONE
# Optimal uniform termination error: 15.7886268841 Params: (6, 18, 0)

SIMULATIONSDICT["reg9"]["Uniform"]["optimal"]["idx_1"] = 6
SIMULATIONSDICT["reg9"]["Uniform"]["optimal"]["idx_2"] = 18
SIMULATIONSDICT["reg9"]["Uniform"]["optimal"]["idx_3"] = 0
SIMULATIONSDICT["reg9"]["Uniform"]["zerolambda"]["idx_1"] = 6
SIMULATIONSDICT["reg9"]["Uniform"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT["reg9"]["Uniform"]["zerolambda"]["idx_3"] = None

SIMULATIONSDICT["reg9"]["TruncGauss"]["optimal"]["idx_1"] = 28
SIMULATIONSDICT["reg9"]["TruncGauss"]["optimal"]["idx_2"] = 23
SIMULATIONSDICT["reg9"]["TruncGauss"]["optimal"]["idx_3"] = 0
SIMULATIONSDICT["reg9"]["TruncGauss"]["zerolambda"]["idx_1"] = 28
SIMULATIONSDICT["reg9"]["TruncGauss"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT["reg9"]["TruncGauss"]["zerolambda"]["idx_3"] = None


SIMULATIONSDICT["regcoarse"]["Opt_Beta_Expn"] = "TruncGauss"
# Optimal trunc termination error: 15.5288221825 Params: (2, 18, 3) DONE
# Optimal uniform termination error: 16.4057729368 Params: (10, 13, 2) 

SIMULATIONSDICT["regcoarse"]["Uniform"]["optimal"]["idx_1"] = 10
SIMULATIONSDICT["regcoarse"]["Uniform"]["optimal"]["idx_2"] = 13
SIMULATIONSDICT["regcoarse"]["Uniform"]["optimal"]["idx_3"] = 2
SIMULATIONSDICT["regcoarse"]["Uniform"]["zerolambda"]["idx_1"] = 10
SIMULATIONSDICT["regcoarse"]["Uniform"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT["regcoarse"]["Uniform"]["zerolambda"]["idx_3"] = None

SIMULATIONSDICT["regcoarse"]["TruncGauss"]["optimal"]["idx_1"] = 2
SIMULATIONSDICT["regcoarse"]["TruncGauss"]["optimal"]["idx_2"] = 18
SIMULATIONSDICT["regcoarse"]["TruncGauss"]["optimal"]["idx_3"] = 3
SIMULATIONSDICT["regcoarse"]["TruncGauss"]["zerolambda"]["idx_1"] = 2
SIMULATIONSDICT["regcoarse"]["TruncGauss"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT["regcoarse"]["TruncGauss"]["zerolambda"]["idx_3"] = None


SIMULATIONSDICT["reg36"]["Opt_Beta_Expn"] = "Uniform"
# Optimal trunc termination error: 21.8931761221 Params: (3, 12, 0) - DONE
# Optimal uniform termination error: 21.6740412099 Params: (29, 13, 0)

SIMULATIONSDICT["reg36"]["Uniform"]["optimal"]["idx_1"] = 29
SIMULATIONSDICT["reg36"]["Uniform"]["optimal"]["idx_2"] = 13
SIMULATIONSDICT["reg36"]["Uniform"]["optimal"]["idx_3"] = 0
SIMULATIONSDICT["reg36"]["Uniform"]["zerolambda"]["idx_1"] = 29
SIMULATIONSDICT["reg36"]["Uniform"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT["reg36"]["Uniform"]["zerolambda"]["idx_3"] = None

SIMULATIONSDICT["reg36"]["TruncGauss"]["optimal"]["idx_1"] = 3
SIMULATIONSDICT["reg36"]["TruncGauss"]["optimal"]["idx_2"] = 12
SIMULATIONSDICT["reg36"]["TruncGauss"]["optimal"]["idx_3"] = 0
SIMULATIONSDICT["reg36"]["TruncGauss"]["zerolambda"]["idx_1"] = 3
SIMULATIONSDICT["reg36"]["TruncGauss"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT["reg36"]["TruncGauss"]["zerolambda"]["idx_3"] = None


SIMULATIONSDICT["regfine"]["Opt_Beta_Expn"] = "TruncGauss"
# Optimal trunc termination error: 25.7545311417 Params: (8, 18, 3) - DONE
# Optimal uniform termination error: 26.9271985266 Params: (29, 12, 0)

SIMULATIONSDICT["regfine"]["Uniform"]["optimal"]["idx_1"] = 29
SIMULATIONSDICT["regfine"]["Uniform"]["optimal"]["idx_2"] = 12
SIMULATIONSDICT["regfine"]["Uniform"]["optimal"]["idx_3"] = 0
SIMULATIONSDICT["regfine"]["Uniform"]["zerolambda"]["idx_1"] = 29
SIMULATIONSDICT["regfine"]["Uniform"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT["regfine"]["Uniform"]["zerolambda"]["idx_3"] = None

SIMULATIONSDICT["regfine"]["TruncGauss"]["optimal"]["idx_1"] = 8
SIMULATIONSDICT["regfine"]["TruncGauss"]["optimal"]["idx_2"] = 18
SIMULATIONSDICT["regfine"]["TruncGauss"]["optimal"]["idx_3"] = 3
SIMULATIONSDICT["regfine"]["TruncGauss"]["zerolambda"]["idx_1"] = 8
SIMULATIONSDICT["regfine"]["TruncGauss"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT["regfine"]["TruncGauss"]["zerolambda"]["idx_3"] = None


SIMULATIONSDICT[1]["Opt_Beta_Expn"] = "TruncGauss"
# Optimal trunc termination error: 20.7438254644 Params: (10, 13, 4) - DONE
# Optimal uniform termination error: 20.9343160986 Params: (29, 18, 4)

SIMULATIONSDICT[1]["Uniform"]["optimal"]["idx_1"] = 29
SIMULATIONSDICT[1]["Uniform"]["optimal"]["idx_2"] = 18
SIMULATIONSDICT[1]["Uniform"]["optimal"]["idx_3"] = 4
SIMULATIONSDICT[1]["Uniform"]["zerolambda"]["idx_1"] = 29
SIMULATIONSDICT[1]["Uniform"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT[1]["Uniform"]["zerolambda"]["idx_3"] = None

SIMULATIONSDICT[1]["TruncGauss"]["optimal"]["idx_1"] = 10
SIMULATIONSDICT[1]["TruncGauss"]["optimal"]["idx_2"] = 13
SIMULATIONSDICT[1]["TruncGauss"]["optimal"]["idx_3"] = 4
SIMULATIONSDICT[1]["TruncGauss"]["zerolambda"]["idx_1"] = 10
SIMULATIONSDICT[1]["TruncGauss"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT[1]["TruncGauss"]["zerolambda"]["idx_3"] = None


SIMULATIONSDICT[2]["Opt_Beta_Expn"] = "Uniform"
# Optimal trunc termination error: 15.8254356917 Params: (1, 20, 1) - DONE
# Optimal uniform termination error: 15.6766254203 Params: (13, 13, 0)

SIMULATIONSDICT[2]["Uniform"]["optimal"]["idx_1"] = 13
SIMULATIONSDICT[2]["Uniform"]["optimal"]["idx_2"] = 13
SIMULATIONSDICT[2]["Uniform"]["optimal"]["idx_3"] = 0
SIMULATIONSDICT[2]["Uniform"]["zerolambda"]["idx_1"] = 13
SIMULATIONSDICT[2]["Uniform"]["zerolambda"]["idx_2"] = 0
SIMULATIONSDICT[2]["Uniform"]["zerolambda"]["idx_3"] = None

SIMULATIONSDICT[2]["TruncGauss"]["optimal"]["idx_1"] = 1
SIMULATIONSDICT[2]["TruncGauss"]["optimal"]["idx_2"] = 20
SIMULATIONSDICT[2]["TruncGauss"]["optimal"]["idx_3"] = 1
SIMULATIONSDICT[2]["TruncGauss"]["zerolambda"]["idx_1"] = 1
SIMULATIONSDICT[2]["TruncGauss"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT[2]["TruncGauss"]["zerolambda"]["idx_3"] = None


SIMULATIONSDICT[3]["Opt_Beta_Expn"] = "TruncGauss"
# Optimal trunc termination error: 15.4765167883 Params: (7, 18, 4) - DONE
# Optimal uniform termination error: 15.6477518697 Params: (6, 7, 0)

SIMULATIONSDICT[3]["Uniform"]["optimal"]["idx_1"] = 6
SIMULATIONSDICT[3]["Uniform"]["optimal"]["idx_2"] = 7
SIMULATIONSDICT[3]["Uniform"]["optimal"]["idx_3"] = 0
SIMULATIONSDICT[3]["Uniform"]["zerolambda"]["idx_1"] = 6
SIMULATIONSDICT[3]["Uniform"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT[3]["Uniform"]["zerolambda"]["idx_3"] = None

SIMULATIONSDICT[3]["TruncGauss"]["optimal"]["idx_1"] = 7
SIMULATIONSDICT[3]["TruncGauss"]["optimal"]["idx_2"] = 18
SIMULATIONSDICT[3]["TruncGauss"]["optimal"]["idx_3"] = 4
SIMULATIONSDICT[3]["TruncGauss"]["zerolambda"]["idx_1"] = 7
SIMULATIONSDICT[3]["TruncGauss"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT[3]["TruncGauss"]["zerolambda"]["idx_3"] = None


SIMULATIONSDICT[4]["Opt_Beta_Expn"] = "TruncGauss"
# Optimal trunc termination error: 16.0805297721 Params: (29, 13, 3) - DONE
# Optimal uniform termination error: 16.0981723624 Params: (24, 8, 2)


SIMULATIONSDICT[4]["Uniform"]["optimal"]["idx_1"] = 24
SIMULATIONSDICT[4]["Uniform"]["optimal"]["idx_2"] = 8
SIMULATIONSDICT[4]["Uniform"]["optimal"]["idx_3"] = 2
SIMULATIONSDICT[4]["Uniform"]["zerolambda"]["idx_1"] = 24
SIMULATIONSDICT[4]["Uniform"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT[4]["Uniform"]["zerolambda"]["idx_3"] = None

SIMULATIONSDICT[4]["TruncGauss"]["optimal"]["idx_1"] = 29
SIMULATIONSDICT[4]["TruncGauss"]["optimal"]["idx_2"] = 13
SIMULATIONSDICT[4]["TruncGauss"]["optimal"]["idx_3"] = 3
SIMULATIONSDICT[4]["TruncGauss"]["zerolambda"]["idx_1"] = 29
SIMULATIONSDICT[4]["TruncGauss"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT[4]["TruncGauss"]["zerolambda"]["idx_3"] = None


SIMULATIONSDICT[5]["Opt_Beta_Expn"] = "Uniform"
# Optimal trunc termination error: 17.6846170351 Params: (4, 18, 4) - DONE
# Optimal uniform termination error: 17.3214186092 Params: (29, 18, 1)


SIMULATIONSDICT[5]["Uniform"]["optimal"]["idx_1"] = 29
SIMULATIONSDICT[5]["Uniform"]["optimal"]["idx_2"] = 18
SIMULATIONSDICT[5]["Uniform"]["optimal"]["idx_3"] = 1
SIMULATIONSDICT[5]["Uniform"]["zerolambda"]["idx_1"] = 29
SIMULATIONSDICT[5]["Uniform"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT[5]["Uniform"]["zerolambda"]["idx_3"] = None

SIMULATIONSDICT[5]["TruncGauss"]["optimal"]["idx_1"] = 4
SIMULATIONSDICT[5]["TruncGauss"]["optimal"]["idx_2"] = 18
SIMULATIONSDICT[5]["TruncGauss"]["optimal"]["idx_3"] = 4
SIMULATIONSDICT[5]["TruncGauss"]["zerolambda"]["idx_1"] = 4
SIMULATIONSDICT[5]["TruncGauss"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT[5]["TruncGauss"]["zerolambda"]["idx_3"] = None


SIMULATIONSDICT[10]["Opt_Beta_Expn"] = "Uniform"
# Optimal trunc termination error: 29.7479929433 Params: (3, 18, 0) - DONE
# Optimal uniform termination error: 28.5415876627 Params: (29, 18, 0)


SIMULATIONSDICT[10]["Uniform"]["optimal"]["idx_1"] = 29
SIMULATIONSDICT[10]["Uniform"]["optimal"]["idx_2"] = 18
SIMULATIONSDICT[10]["Uniform"]["optimal"]["idx_3"] = 0
SIMULATIONSDICT[10]["Uniform"]["zerolambda"]["idx_1"] = 29
SIMULATIONSDICT[10]["Uniform"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT[10]["Uniform"]["zerolambda"]["idx_3"] = None

SIMULATIONSDICT[10]["TruncGauss"]["optimal"]["idx_1"] = 3
SIMULATIONSDICT[10]["TruncGauss"]["optimal"]["idx_2"] = 18
SIMULATIONSDICT[10]["TruncGauss"]["optimal"]["idx_3"] = 0
SIMULATIONSDICT[10]["TruncGauss"]["zerolambda"]["idx_1"] = 3
SIMULATIONSDICT[10]["TruncGauss"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT[10]["TruncGauss"]["zerolambda"]["idx_3"] = None



