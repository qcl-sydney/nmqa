import sys
sys.path.append('../paduaq')
from pdpoints import dims_padua_set

SIMULATIONSDICT = {}

for padua_order in [2, 3, 4, 5] : # Padua order for linear function
    
    SIMULATIONSDICT[padua_order] = {}
    
    SIMULATIONSDICT[padua_order]["num_of_nodes"] = dims_padua_set(padua_order) + 25
    SIMULATIONSDICT[padua_order]["linear"] = False
    SIMULATIONSDICT[padua_order]["repts"] = 50
    SIMULATIONSDICT[padua_order]["max_iterations"] = dims_padua_set(padua_order) * 5
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
SIMULATIONSDICT[2]["Opt_Beta_Expn"] = "Uniform"

# Optimal trunc termination error: 18.586002274 Params: (21, 12, 0)
# Optimal uniform termination error: 18.4668149842 Params: (23, 20, 0)

# *_Rlin_padua_ord_2_rand_21_12_*.npz 
# *_Rlin_padua_ord_2_rand_21_0_*.npz
# *_Rlin_padua_ord_2_rand_23_20_*.npz 
# *_Rlin_padua_ord_2_rand_23_0_*.npz 

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
# *_Rlin_padua_ord_3_rand_4_0_*.npz
# *_Rlin_padua_ord_3_rand_28_13_*.npz 
# *_Rlin_padua_ord_3_rand_28_0_*.npz 

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
# not done

SIMULATIONSDICT[4]["Uniform"]["optimal"]["idx_1"] = 9
SIMULATIONSDICT[4]["Uniform"]["optimal"]["idx_2"] = 12
SIMULATIONSDICT[4]["Uniform"]["optimal"]["idx_3"] = 2
SIMULATIONSDICT[4]["Uniform"]["zerolambda"]["idx_1"] = 9
SIMULATIONSDICT[4]["Uniform"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT[4]["Uniform"]["zerolambda"]["idx_3"] = None

SIMULATIONSDICT[4]["TruncGauss"]["optimal"]["idx_1"] = 29
SIMULATIONSDICT[4]["TruncGauss"]["optimal"]["idx_2"] = 8
SIMULATIONSDICT[4]["TruncGauss"]["optimal"]["idx_3"] = 0 
SIMULATIONSDICT[4]["TruncGauss"]["zerolambda"]["idx_1"] = 29
SIMULATIONSDICT[4]["TruncGauss"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT[4]["TruncGauss"]["zerolambda"]["idx_3"] = None

SIMULATIONSDICT[5]["Opt_Beta_Expn"] = "Uniform"
# not done

SIMULATIONSDICT[5]["Uniform"]["optimal"]["idx_1"] = 17 
SIMULATIONSDICT[5]["Uniform"]["optimal"]["idx_2"] = 18
SIMULATIONSDICT[5]["Uniform"]["optimal"]["idx_3"] = 2
SIMULATIONSDICT[5]["Uniform"]["zerolambda"]["idx_1"] = 17
SIMULATIONSDICT[5]["Uniform"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT[5]["Uniform"]["zerolambda"]["idx_3"] = None

SIMULATIONSDICT[5]["TruncGauss"]["optimal"]["idx_1"] = 3
SIMULATIONSDICT[5]["TruncGauss"]["optimal"]["idx_2"] = 28
SIMULATIONSDICT[5]["TruncGauss"]["optimal"]["idx_3"] = 0 
SIMULATIONSDICT[5]["TruncGauss"]["zerolambda"]["idx_1"] = 3
SIMULATIONSDICT[5]["TruncGauss"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT[5]["TruncGauss"]["zerolambda"]["idx_3"] = None
