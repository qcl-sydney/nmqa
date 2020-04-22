import sys
import os
import copy
import traceback
import numpy as np
# import matplotlib.pyplot as plt

sys.path.append('../qslam')
sys.path.append('../paduaq')
from pdinter_MM import pd_interpolant, calc_padua_cgl
from true_functions import true_function, generate_data_qubits_coords

########################
# Find qslam modules
########################

from singlerun import SingleRunAnalysis
from riskanalysis import EngineeredTruth
from qslamdesignparams import GLOBALDICT


########################
# Taking in bash parameters
########################
padua_order = int(sys.argv[1]) # Padua order. 1, 2, 3, 4, 5,...
idx_job_array = int(sys.argv[2]) # job array starts from 1. Used for tuning
idx_functype = int(sys.argv[3])
data_qubit_num = 25
data_qubit_flag ='uniform'

TUNING_MULTIPLIER=5

########################
# Set True Field
########################

if idx_functype ==0:
    true_function_type = 'cheb2fun'
    from tuningresults_nonpoly import SIMULATIONSDICT
    
if idx_functype ==1:
    true_function_type = 'lin'
    from tuningresults_linear import SIMULATIONSDICT

if idx_functype ==2:
    true_function_type = 'franke'
    from tuningresults_franke import SIMULATIONSDICT
    
if idx_functype ==3:
    true_function_type = 'franke_2'
    from tuningresults_franke_2 import SIMULATIONSDICT

if idx_functype ==4:
    true_function_type = 'gss'
    from tuningresults_gss import SIMULATIONSDICT
       
########################
# Save to path 
########################

path = '/scratch/QCL_RG/qslam_padua_paper/' # on Artemis './data/'

########################
# Generate Sensing Qubits
########################

# Sensor-qubits in Padua formation

if padua_order > 0:
    sensing_qubits = calc_padua_cgl(padua_order)[0]
    
# No data-qubits, basic analysis

if padua_order == -1: 
    sensing_qubits = generate_data_qubits_coords(data_qubit_num,
                                          flag=data_qubit_flag)

# Sensor-qubits in regular (non-Padua) formation

if padua_order == -2:
    FINEGRID = 81
    sensing_qubits = generate_data_qubits_coords(FINEGRID, flag=data_qubit_flag)

if padua_order == -3:
    COARSEGRID = 16
    sensing_qubits = generate_data_qubits_coords(COARSEGRID, flag=data_qubit_flag)
    
    # Re-position grid inside square region
    sensing_qubits = list(np.asarray(sensing_qubits) * 0.75)

if padua_order == -4:
    REG4 = 4
    sensing_qubits = generate_data_qubits_coords(REG4, flag=data_qubit_flag)
    sensing_qubits = list(np.asarray(sensing_qubits) * 0.75)

if padua_order == -5:
    REG9 = 9
    sensing_qubits = generate_data_qubits_coords(REG9, flag=data_qubit_flag)
    sensing_qubits = list(np.asarray(sensing_qubits) * 0.75)

if padua_order == -6:
    REG36 = 36
    sensing_qubits = generate_data_qubits_coords(REG36, flag=data_qubit_flag)

########################
# Generate Data Qubits
########################

if padua_order > 0:
    data_qubits = generate_data_qubits_coords(data_qubit_num, flag=data_qubit_flag)

    GLOBALDICT["DATA_QUBITS"] = np.arange(len(sensing_qubits),  len(sensing_qubits) + data_qubit_num, dtype='int')
    GLOBALDICT["INTERPOLATE_FLAG"] = padua_order
    prefix = true_function_type +'_padua_ord_'+str(padua_order)+'_'

if padua_order == -1: 
    GLOBALDICT["DATA_QUBITS"] = []
    GLOBALDICT["INTERPOLATE_FLAG"] = None
    prefix = true_function_type +'_no_padua_'

if padua_order == -2:
    data_qubits = generate_data_qubits_coords(data_qubit_num, flag=data_qubit_flag)
    # remove duplicate sensors:
    sensing_qubits = list(set(sensing_qubits) - set(data_qubits))

    # update dictionary params:
    GLOBALDICT["DATA_QUBITS"] = np.arange(len(sensing_qubits), len(sensing_qubits) + data_qubit_num, dtype='int')
    GLOBALDICT["INTERPOLATE_FLAG"] = 'Rbf' # change this to 'Rbf'
    prefix = true_function_type +'_regfine_'
    
    # reset key for SIMULATIONSDICT
    padua_order = "regfine"

if padua_order <= -3:
    data_qubits = generate_data_qubits_coords(data_qubit_num, flag=data_qubit_flag)
    GLOBALDICT["DATA_QUBITS"] = np.arange(len(sensing_qubits),  len(sensing_qubits) + data_qubit_num, dtype='int')
    GLOBALDICT["INTERPOLATE_FLAG"] = 'Rbf' 
    
if padua_order == -3:
    prefix = true_function_type +'_regcoarse_'
    padua_order = "regcoarse"
    
if padua_order == -4:
    prefix = true_function_type +'_reg4_'
    padua_order = "reg4"
    
if padua_order == -5:
    prefix = true_function_type +'_reg9_'
    padua_order = "reg9"

if padua_order == -6:
    prefix = true_function_type +'_reg36_'
    padua_order = "reg36"
    
########################
# Set true map and qubit grid
########################
TRUTHKWARGS = {}
TRUTHKWARGS["truthtype"] = "UseFunction"
TRUTHKWARGS["true_function"] = true_function
TRUTHKWARGS["true_function_type"] = true_function_type

TRUTHKWARGS["all_qubit_locations"] = sensing_qubits 
if len(GLOBALDICT["DATA_QUBITS"]) > 0:
    TRUTHKWARGS["all_qubit_locations"] = sensing_qubits + data_qubits

num_of_nodes = len(TRUTHKWARGS["all_qubit_locations"])
true_map_ =  EngineeredTruth(num_of_nodes, TRUTHKWARGS).get_map()

GLOBALDICT["GRIDDICT"] = {}
for idx_posy in range(num_of_nodes):
    
    point = TRUTHKWARGS["all_qubit_locations"][idx_posy]

    if idx_posy < 10 :

        GLOBALDICT["GRIDDICT"]["QUBIT_0" + str(idx_posy)] = (point[0], point[1])

    if idx_posy >= 10 :

        GLOBALDICT["GRIDDICT"]["QUBIT_" + str(idx_posy)] =  (point[0], point[1])

########################
# Set Defaults
########################

GLOBALDICT["MODELDESIGN"]["MSMTS_PER_NODE"] = 1
GLOBALDICT["MODELDESIGN"]["MULTIPLER_R_MAX"] = 4.

repts = 50
particleconfigs = [ [3,2], [9,6], [15,10], [21,14], [30, 20]]


lambda_paris_2 = np.load('lambda_pairs_2.npz')
random_variances = np.load('random_variances.npz')

IDX1_SHP = len(random_variances['g1var'])
IDX2_SHP = 30 # len(lambda_paris_2['lambda_1']) not doing all 250
IDX3_SHP = len(particleconfigs)

####################################
# Tuning Script (Data Generation)
####################################


idx_1, idx_2 = np.unravel_index(idx_job_array - 1 , (IDX1_SHP, IDX2_SHP) )

GLOBALDICT["NOISEPARAMS"]["SIGMOID_APPROX_ERROR"]["SIGMA"] = random_variances['g2var'][idx_1]
GLOBALDICT["NOISEPARAMS"]["QUANTISATION_UNCERTY"]["SIGMA"] = random_variances['g1var'][idx_1]
GLOBALDICT["MODELDESIGN"]["LAMBDA_1"] = lambda_paris_2['lambda_1'][idx_2]
GLOBALDICT["MODELDESIGN"]["LAMBDA_2"] = lambda_paris_2['lambda_2'][idx_2]

fname_likelihood = 'rand_'+str(idx_1)+'_'+str(idx_2)+'_'

max_iterations = len(sensing_qubits) * TUNING_MULTIPLIER
GLOBALDICT["MODELDESIGN"]["MAX_NUM_ITERATIONS"] = max_iterations

for idx_3 in range(IDX3_SHP):

    GLOBALDICT["MODELDESIGN"]["P_ALPHA"] = particleconfigs[idx_3][0]
    GLOBALDICT["MODELDESIGN"]["P_BETA"] = particleconfigs[idx_3][1]

    SAMPLE_GLOBAL_MODEL = copy.deepcopy(GLOBALDICT)

    uniform_r_expt = SingleRunAnalysis(SAMPLE_GLOBAL_MODEL, true_map_, repts, beta_expansion_mode=False, beta_skew_adjust=False)
    uniform_r_expt.run_analysis(path+'Uni_R'+prefix+fname_likelihood+str(idx_3))

    trunc_r_expt = SingleRunAnalysis(SAMPLE_GLOBAL_MODEL, true_map_, repts, beta_expansion_mode=True, beta_skew_adjust=False)
    trunc_r_expt.run_analysis(path+'Trunc_R'+prefix+fname_likelihood+str(idx_3))

