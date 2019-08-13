import sys
sys.path.append('../paduaq')
from pdpoints import dims_padua_set

SIMULATIONSDICT = {}

for padua_order in [2, 3, 4, 5] : # Padua order for cheb2chev function
    
    SIMULATIONSDICT[padua_order] = {}
    
    SIMULATIONSDICT[padua_order]["num_of_nodes"] = dims_padua_set(padua_order) + 25
    SIMULATIONSDICT[padua_order]["linear"] = False
    SIMULATIONSDICT[padua_order]["repts"] = 50
    SIMULATIONSDICT[padua_order]["max_iterations"] = dims_padua_set(padua_order) * 5s
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
SIMULATIONSDICT[2]["Opt_Beta_Expn"] = "TruncGauss"

# Optimal trunc termination error: 17.675734373 Params: (2, 13, 4)
# Optimal uniform termination error: 17.7089964369 Params: (16, 7, 0)

# *_Rcheb2fun_padua_ord_2_rand_2_13_*.npz 
# *_Rcheb2fun_padua_ord_2_rand_2_0_*.npz
# *_Rcheb2fun_padua_ord_2_rand_16_7_*.npz 
# *_Rcheb2fun_padua_ord_2_rand_16_0_*.npz 

SIMULATIONSDICT[2]["Uniform"]["optimal"]["idx_1"] = 16
SIMULATIONSDICT[2]["Uniform"]["optimal"]["idx_2"] = 7
SIMULATIONSDICT[2]["Uniform"]["optimal"]["idx_3"] = 0
SIMULATIONSDICT[2]["Uniform"]["zerolambda"]["idx_1"] = 16
SIMULATIONSDICT[2]["Uniform"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT[2]["Uniform"]["zerolambda"]["idx_3"] = None

SIMULATIONSDICT[2]["TruncGauss"]["optimal"]["idx_1"] = 2
SIMULATIONSDICT[2]["TruncGauss"]["optimal"]["idx_2"] = 13
SIMULATIONSDICT[2]["TruncGauss"]["optimal"]["idx_3"] = 4
SIMULATIONSDICT[2]["TruncGauss"]["zerolambda"]["idx_1"] = 2
SIMULATIONSDICT[2]["TruncGauss"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT[2]["TruncGauss"]["zerolambda"]["idx_3"] = None


SIMULATIONSDICT[3]["Opt_Beta_Expn"] = "Uniform"
# not done

SIMULATIONSDICT[3]["Uniform"]["optimal"]["idx_1"] = 21
SIMULATIONSDICT[3]["Uniform"]["optimal"]["idx_2"] = 20
SIMULATIONSDICT[3]["Uniform"]["optimal"]["idx_3"] = 2
SIMULATIONSDICT[3]["Uniform"]["zerolambda"]["idx_1"] = 21
SIMULATIONSDICT[3]["Uniform"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT[3]["Uniform"]["zerolambda"]["idx_3"] = None

SIMULATIONSDICT[3]["TruncGauss"]["optimal"]["idx_1"] = 26
SIMULATIONSDICT[3]["TruncGauss"]["optimal"]["idx_2"] = 18
SIMULATIONSDICT[3]["TruncGauss"]["optimal"]["idx_3"] = 3
SIMULATIONSDICT[3]["TruncGauss"]["zerolambda"]["idx_1"] = 26
SIMULATIONSDICT[3]["TruncGauss"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT[3]["TruncGauss"]["zerolambda"]["idx_3"] = None


SIMULATIONSDICT[4]["Opt_Beta_Expn"] = "TruncGauss"
# not done

SIMULATIONSDICT[4]["Uniform"]["optimal"]["idx_1"] = 17
SIMULATIONSDICT[4]["Uniform"]["optimal"]["idx_2"] = 8
SIMULATIONSDICT[4]["Uniform"]["optimal"]["idx_3"] = 2
SIMULATIONSDICT[4]["Uniform"]["zerolambda"]["idx_1"] = 17
SIMULATIONSDICT[4]["Uniform"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT[4]["Uniform"]["zerolambda"]["idx_3"] = None

SIMULATIONSDICT[4]["TruncGauss"]["optimal"]["idx_1"] = 24
SIMULATIONSDICT[4]["TruncGauss"]["optimal"]["idx_2"] = 15
SIMULATIONSDICT[4]["TruncGauss"]["optimal"]["idx_3"] = 0 
SIMULATIONSDICT[4]["TruncGauss"]["zerolambda"]["idx_1"] = 24
SIMULATIONSDICT[4]["TruncGauss"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT[4]["TruncGauss"]["zerolambda"]["idx_3"] = None


SIMULATIONSDICT[5]["Opt_Beta_Expn"] = "Uniform"
# not done

SIMULATIONSDICT[5]["Uniform"]["optimal"]["idx_1"] = 25 
SIMULATIONSDICT[5]["Uniform"]["optimal"]["idx_2"] = 8
SIMULATIONSDICT[5]["Uniform"]["optimal"]["idx_3"] = 3
SIMULATIONSDICT[5]["Uniform"]["zerolambda"]["idx_1"] = 25
SIMULATIONSDICT[5]["Uniform"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT[5]["Uniform"]["zerolambda"]["idx_3"] = None

SIMULATIONSDICT[5]["TruncGauss"]["optimal"]["idx_1"] = 3
SIMULATIONSDICT[5]["TruncGauss"]["optimal"]["idx_2"] = 12
SIMULATIONSDICT[5]["TruncGauss"]["optimal"]["idx_3"] = 0 
SIMULATIONSDICT[5]["TruncGauss"]["zerolambda"]["idx_1"] = 3
SIMULATIONSDICT[5]["TruncGauss"]["zerolambda"]["idx_2"] = 0 
SIMULATIONSDICT[5]["TruncGauss"]["zerolambda"]["idx_3"] = None
