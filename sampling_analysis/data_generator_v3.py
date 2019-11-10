import sys
import os
import copy
import traceback



########################
# Find qslam modules
########################
sys.path.append('../qslam/')

from qslamdesignparams import GLOBALDICT
from riskanalysis import CreateQslamExpt as riskqslam
from riskanalysis import CreateNaiveExpt as risknaive
from riskanalysis import EngineeredTruth
from visualiserisk import *

sys.path.append('../qslam')
sys.path.append('../paduaq')
from pdinter_MM import pd_interpolant, calc_padua_cgl
from true_functions import true_function, generate_data_qubits_coords

sys.path.append('./')


########################
# Taking in bash parameters
########################
padua_order = int(sys.argv[1])  # Padua order 1, 2, 3, 4, 5 for lin; 3, 4, 5 for cheb2fun

idx_functype = int(sys.argv[2])
if idx_functype ==0:
    true_function_type = 'cheb2fun'
    from tuningresults_cheb2fun_5 import SIMULATIONSDICT
if idx_functype ==1:
    true_function_type = 'lin'
    from tuningresults_lin_5 import SIMULATIONSDICT
 
data_qubit_num = 25
data_qubit_flag ='uniform'

########################
# Save to path 
########################

path = '/scratch/QCL_RG/qslam_padua_paper/' # on Artemis './data/'

########################
# Generate Sensor Qubits
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

########################
# Generate Data Qubits
########################

if padua_order > 0:
    data_qubits = generate_data_qubits_coords(data_qubit_num, flag=data_qubit_flag)

    GLOBALDICT["DATA_QUBITS"] = np.arange(len(sensing_qubits),  len(sensing_qubits) + data_qubit_num, dtype='int')
    GLOBALDICT["INTERPOLATE_FLAG"] = padua_order
    prefix = true_function_type +'_padua_ord_'+str(padua_order)+'_'

if padua_order == -1: 
    GLOBALDICT["DATA_QUBITS"] = None
    GLOBALDICT["INTERPOLATE_FLAG"] = None
    prefix = true_function_type +'_no_padua_'
    padua_order = "no_padua"

if padua_order == -2:
    data_qubits = generate_data_qubits_coords(data_qubit_num, flag=data_qubit_flag)
    
    # remove duplicate sensors:
    sensing_qubits = list(set(sensing_qubits) - set(data_qubits))
    
    # update dictionary params:
    GLOBALDICT["DATA_QUBITS"] = np.arange(len(sensing_qubits), len(sensing_qubits) + data_qubit_num, dtype='int')
    GLOBALDICT["INTERPOLATE_FLAG"] = 'linear'
    prefix = true_function_type +'_regfine_'
    
    # reset key for SIMULATIONSDICT
    padua_order = "regfine"

if padua_order == -3:
    data_qubits = generate_data_qubits_coords(data_qubit_num, flag=data_qubit_flag)
    
    # update dictionary params:
    GLOBALDICT["DATA_QUBITS"] = np.arange(len(sensing_qubits),  len(sensing_qubits) + data_qubit_num, dtype='int')
    GLOBALDICT["INTERPOLATE_FLAG"] = 'linear'
    prefix = true_function_type +'_regcoarse_'
    
    # reset key for SIMULATIONSDICT
    padua_order = "regcoarse"
    
    
########################
# Set hardware and true map
########################
TRUTHKWARGS = {}
TRUTHKWARGS["truthtype"] = "UseFunction"
TRUTHKWARGS["true_function"] = true_function
TRUTHKWARGS["true_function_type"] = true_function_type

TRUTHKWARGS["all_qubit_locations"] = sensing_qubits 
if GLOBALDICT["DATA_QUBITS"] is not None:
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
# Set Simulation Params
########################

GLOBALDICT["MODELDESIGN"]["MSMTS_PER_NODE"] = 1
GLOBALDICT["MODELDESIGN"]["MULTIPLER_R_MAX"] = 4.
repts = 50
particleconfigs = [ [3,2], [9,6], [15,10], [21,14], [30, 20]]
lambda_paris_2 = np.load('lambda_pairs_2.npz')
random_variances = np.load('random_variances.npz')

########################
# Set Loop Parameters
########################

Multiples = [1, 3, 5, 7, 9, 10, 15, 20, 50]

########################
# Run Script
######## ################
  

opt_method = "Uniform" # SIMULATIONSDICT[padua_order]["Opt_Beta_Expn"] 

idx_1 = SIMULATIONSDICT[padua_order][opt_method]["optimal"]["idx_1"]
idx_2 = SIMULATIONSDICT[padua_order][opt_method]["optimal"]["idx_2"]
idx_3 = SIMULATIONSDICT[padua_order][opt_method]["optimal"]["idx_3"]

GLOBALDICT["NOISEPARAMS"]["SIGMOID_APPROX_ERROR"]["SIGMA"] = random_variances['g2var'][idx_1]
GLOBALDICT["NOISEPARAMS"]["QUANTISATION_UNCERTY"]["SIGMA"] = random_variances['g1var'][idx_1]
GLOBALDICT["MODELDESIGN"]["LAMBDA_1"] = lambda_paris_2['lambda_1'][idx_2]
GLOBALDICT["MODELDESIGN"]["LAMBDA_2"] = lambda_paris_2['lambda_2'][idx_2]
GLOBALDICT["MODELDESIGN"]["P_ALPHA"] = particleconfigs[idx_3][0]
GLOBALDICT["MODELDESIGN"]["P_BETA"] = particleconfigs[idx_3][1]

fname_likelihood = 'optidx_'+str(idx_1)+'_'+str(idx_2)+'_'+str(idx_3)

for idx_msmt_var in Multiples:

    unique_id = path + prefix + fname_likelihood + '_m_' + str(idx_msmt_var)
    
    GLOBALDICT["MODELDESIGN"]["MAX_NUM_ITERATIONS"] = len(sensing_qubits) * idx_msmt_var
    GLOBALDICT["MODELDESIGN"]["ID"] = unique_id

    qslam_br = 0.
    naive_br = 0.
    qslamdata = 0.
    naivedata = 0.

    try:
        qslam_br = riskqslam(copy.deepcopy(TRUTHKWARGS), copy.deepcopy(GLOBALDICT))
        naive_br = risknaive(copy.deepcopy(TRUTHKWARGS), copy.deepcopy(GLOBALDICT))
        qslam_br.naive_implementation(randomise='OFF')
        naive_br.naive_implementation()
        

    except:
        print "Index: %s was not completed..." %(idx_msmt_var)
        print "Error information:"
        print "Type", sys.exc_info()[0]
        print "Value", sys.exc_info()[1]
        print "Traceback", traceback.format_exc()
