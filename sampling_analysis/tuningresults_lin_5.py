import sys
sys.path.append('../paduaq')
from pdpoints import dims_padua_set

SIMULATIONSDICT = {}
GRIDSIZE=25 # number of data qubits if padua qubits are used. else, total number of regularly arranged sensing qubits
MULTIPLIER=5

for padua_order in ["no_padua", 2, 3, 4, 5, 10] : # Padua order for linear function
    
    SIMULATIONSDICT[padua_order] = {}
    
    if padua_order == "no_padua":
        SIMULATIONSDICT[padua_order]["max_iterations"] = GRIDSIZE * MULTIPLIER
        SIMULATIONSDICT[padua_order]["num_of_nodes"] = GRIDSIZE
    
    if padua_order != "no_padua":
        SIMULATIONSDICT[padua_order]["max_iterations"] = dims_padua_set(padua_order) * MULTIPLIER
        SIMULATIONSDICT[padua_order]["num_of_nodes"] = dims_padua_set(padua_order) + GRIDSIZE
    
    SIMULATIONSDICT[padua_order]["linear"] = False # WTF is this
    SIMULATIONSDICT[padua_order]["repts"] = 50
    SIMULATIONSDICT[padua_order]["functype"]= 'lin'
    
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
# Optimal trunc termination error: 7.74578403656 Params: (8, 13, 4)
# Optimal uniform termination error: 8.39182734811 Params: (22, 18, 0)

SIMULATIONSDICT["no_padua"]["Uniform"]["optimal"]["idx_1"] = 22
SIMULATIONSDICT["no_padua"]["Uniform"]["optimal"]["idx_2"] = 18
SIMULATIONSDICT["no_padua"]["Uniform"]["optimal"]["idx_3"] = 0
SIMULATIONSDICT["no_padua"]["Uniform"]["zerolambda"]["idx_1"] = 22
SIMULATIONSDICT["no_padua"]["Uniform"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT["no_padua"]["Uniform"]["zerolambda"]["idx_3"] = None

SIMULATIONSDICT["no_padua"]["TruncGauss"]["optimal"]["idx_1"] = 8
SIMULATIONSDICT["no_padua"]["TruncGauss"]["optimal"]["idx_2"] = 13
SIMULATIONSDICT["no_padua"]["TruncGauss"]["optimal"]["idx_3"] = 4
SIMULATIONSDICT["no_padua"]["TruncGauss"]["zerolambda"]["idx_1"] = 8
SIMULATIONSDICT["no_padua"]["TruncGauss"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT["no_padua"]["TruncGauss"]["zerolambda"]["idx_3"] = None


SIMULATIONSDICT[2]["Opt_Beta_Expn"] = "Uniform"

# Optimal trunc termination error: 18.586002274 Params: (21, 12, 0)
# Optimal uniform termination error: 18.4668149842 Params: (23, 20, 0)

# *_Rlin_padua_ord_2_rand_21_12_*.npz done
# *_Rlin_padua_ord_2_rand_21_0_*.npz done
# *_Rlin_padua_ord_2_rand_23_20_*.npz done
# *_Rlin_padua_ord_2_rand_23_0_*.npz done

SIMULATIONSDICT[2]["Uniform"]["optimal"]["idx_1"] = 23
SIMULATIONSDICT[2]["Uniform"]["optimal"]["idx_2"] = 20
SIMULATIONSDICT[2]["Uniform"]["optimal"]["idx_3"] = 0
SIMULATIONSDICT[2]["Uniform"]["zerolambda"]["idx_1"] = 23
SIMULATIONSDICT[2]["Uniform"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT[2]["Uniform"]["zerolambda"]["idx_3"] = None

SIMULATIONSDICT[2]["TruncGauss"]["optimal"]["idx_1"] = 21
SIMULATIONSDICT[2]["TruncGauss"]["optimal"]["idx_2"] = 12
SIMULATIONSDICT[2]["TruncGauss"]["optimal"]["idx_3"] = 0
SIMULATIONSDICT[2]["TruncGauss"]["zerolambda"]["idx_1"] = 21
SIMULATIONSDICT[2]["TruncGauss"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT[2]["TruncGauss"]["zerolambda"]["idx_3"] = None


SIMULATIONSDICT[3]["Opt_Beta_Expn"] = "Uniform"
# Optimal trunc termination error: 18.3636849816 Params: (4, 7, 0)
# Optimal uniform termination error: 18.3882597628 Params: (28, 13, 0)
# *_Rlin_padua_ord_3_rand_4_7_*.npz done
# *_Rlin_padua_ord_3_rand_4_0_*.npz done
# *_Rlin_padua_ord_3_rand_28_13_*.npz done
# *_Rlin_padua_ord_3_rand_28_0_*.npz done

SIMULATIONSDICT[3]["Uniform"]["optimal"]["idx_1"] = 4
SIMULATIONSDICT[3]["Uniform"]["optimal"]["idx_2"] = 7
SIMULATIONSDICT[3]["Uniform"]["optimal"]["idx_3"] = 0
SIMULATIONSDICT[3]["Uniform"]["zerolambda"]["idx_1"] = 4
SIMULATIONSDICT[3]["Uniform"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT[3]["Uniform"]["zerolambda"]["idx_3"] = None

SIMULATIONSDICT[3]["TruncGauss"]["optimal"]["idx_1"] = 28
SIMULATIONSDICT[3]["TruncGauss"]["optimal"]["idx_2"] = 13
SIMULATIONSDICT[3]["TruncGauss"]["optimal"]["idx_3"] = 0
SIMULATIONSDICT[3]["TruncGauss"]["zerolambda"]["idx_1"] = 28
SIMULATIONSDICT[3]["TruncGauss"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT[3]["TruncGauss"]["zerolambda"]["idx_3"] = None


SIMULATIONSDICT[4]["Opt_Beta_Expn"] = "Uniform"
# Optimal trunc termination error: 19.6976266344 Params: (22, 18, 0)
# Optimal uniform termination error: 18.9699583396 Params: (5, 13, 0)

# *_Rlin_padua_ord_4_rand_22_18_*.npz done
# *_Rlin_padua_ord_4_rand_22_0_*.npz done
# *_Rlin_padua_ord_4_rand_5_13_*.npz done
# *_Rlin_padua_ord_4_rand_5_0_*.npz done

SIMULATIONSDICT[4]["Uniform"]["optimal"]["idx_1"] = 5
SIMULATIONSDICT[4]["Uniform"]["optimal"]["idx_2"] = 13
SIMULATIONSDICT[4]["Uniform"]["optimal"]["idx_3"] = 0
SIMULATIONSDICT[4]["Uniform"]["zerolambda"]["idx_1"] = 5
SIMULATIONSDICT[4]["Uniform"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT[4]["Uniform"]["zerolambda"]["idx_3"] = None

SIMULATIONSDICT[4]["TruncGauss"]["optimal"]["idx_1"] = 22
SIMULATIONSDICT[4]["TruncGauss"]["optimal"]["idx_2"] = 18
SIMULATIONSDICT[4]["TruncGauss"]["optimal"]["idx_3"] = 0 
SIMULATIONSDICT[4]["TruncGauss"]["zerolambda"]["idx_1"] = 22
SIMULATIONSDICT[4]["TruncGauss"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT[4]["TruncGauss"]["zerolambda"]["idx_3"] = None

SIMULATIONSDICT[5]["Opt_Beta_Expn"] = "Uniform"
# Optimal trunc termination error: 21.3448381468 Params: (28, 18, 1)
# Optimal uniform termination error: 20.9187338625 Params: (22, 23, 0)
# *_Rlin_padua_ord_5_rand_22_23_*.npz done
# *_Rlin_padua_ord_5_rand_22_0_*.npz done
# *_Rlin_padua_ord_5_rand_28_18_*.npz done
# *_Rlin_padua_ord_5_rand_28_0_*.npz done

SIMULATIONSDICT[5]["Uniform"]["optimal"]["idx_1"] = 22 
SIMULATIONSDICT[5]["Uniform"]["optimal"]["idx_2"] = 23
SIMULATIONSDICT[5]["Uniform"]["optimal"]["idx_3"] = 0
SIMULATIONSDICT[5]["Uniform"]["zerolambda"]["idx_1"] = 22
SIMULATIONSDICT[5]["Uniform"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT[5]["Uniform"]["zerolambda"]["idx_3"] = None

SIMULATIONSDICT[5]["TruncGauss"]["optimal"]["idx_1"] = 28
SIMULATIONSDICT[5]["TruncGauss"]["optimal"]["idx_2"] = 18
SIMULATIONSDICT[5]["TruncGauss"]["optimal"]["idx_3"] = 1
SIMULATIONSDICT[5]["TruncGauss"]["zerolambda"]["idx_1"] = 28
SIMULATIONSDICT[5]["TruncGauss"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT[5]["TruncGauss"]["zerolambda"]["idx_3"] = None


SIMULATIONSDICT[10]["Opt_Beta_Expn"] = "NAME"

SIMULATIONSDICT[10]["Uniform"]["optimal"]["idx_1"] = 0 
SIMULATIONSDICT[10]["Uniform"]["optimal"]["idx_2"] = 0
SIMULATIONSDICT[10]["Uniform"]["optimal"]["idx_3"] = 0
SIMULATIONSDICT[10]["Uniform"]["zerolambda"]["idx_1"] = 0
SIMULATIONSDICT[10]["Uniform"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT[10]["Uniform"]["zerolambda"]["idx_3"] = None

SIMULATIONSDICT[10]["TruncGauss"]["optimal"]["idx_1"] = 0
SIMULATIONSDICT[10]["TruncGauss"]["optimal"]["idx_2"] = 0
SIMULATIONSDICT[10]["TruncGauss"]["optimal"]["idx_3"] = 1
SIMULATIONSDICT[10]["TruncGauss"]["zerolambda"]["idx_1"] = 0
SIMULATIONSDICT[10]["TruncGauss"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT[10]["TruncGauss"]["zerolambda"]["idx_3"] = None
