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
from tuningresults import SIMULATIONSDICT

########################
# Taking in bash parameters
########################
padua_order = int(sys.argv[1]) + 2 # Padua order 3, 4, 5
true_function_type = 'cheb2fun'
data_qubit_num = 25
data_qubit_flag ='uniform'

########################
# Save to path 
########################

path = '/scratch/QCL_RG/qslam_padua_paper/' # on Artemis './data/'

########################
# Generate Padua Qubits
########################

sensing_qubits = calc_padua_cgl(padua_order)[0]

########################
# Generate Data Qubits
########################

data_qubits = generate_data_qubits_coords(data_qubit_num,
                                          flag=data_qubit_flag)

GLOBALDICT["DATA_QUBITS"] = np.arange(len(sensing_qubits),  len(sensing_qubits) + data_qubit_num, dtype='int')
GLOBALDICT["INTERPOLATE_FLAG"] = padua_order

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
prefix = '_padua_ord_'+str(padua_order)+'_'
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
