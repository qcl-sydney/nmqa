import sys
import os
import copy
import traceback
import numpy as np
import matplotlib.pyplot as plt

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
data_qubit_num = 25
data_qubit_flag ='uniform'

########################
# Save to path 
########################

path = '/scratch/QCL_RG/qslam_padua_paper/' # on Artemis './data/'
true_function_type = 'cheb2fun'

########################
# Generate Padua Points
########################

total_qubits = calc_padua_cgl(padua_order)[0]

########################
# Set Data Points
########################

GLOBALDICT["GRIDDICT"]["DATA_QUBITS"] = generate_data_qubits_coords(data_qubit_num,
                                                                    flag=data_qubit_flag)

########################
# Set hardware and true map
########################

if GLOBALDICT["DATA_QUBITS"] is not None:
    total_qubits = total_qubits + GLOBALDICT["DATA_QUBITS"]
num_of_nodes = len(total_qubits)

true_map = []
GLOBALDICT["GRIDDICT"] = {}
for idx_posy in range(num_of_nodes):
    
    point = total_qubits[idx_posy]
    true_map.append(true_function(point[0], point[1], d=true_function_type))

    if idx_posy < 10 :

        GLOBALDICT["GRIDDICT"]["QUBIT_0" + str(idx_posy)] = (point[0], point[1])

    if idx_posy >= 10 :

        GLOBALDICT["GRIDDICT"]["QUBIT_" + str(idx_posy)] =  (point[0], point[1])
        
true_map_ = np.asarray(true_map)

########################
# Set Defaults
########################

GLOBALDICT["MODELDESIGN"]["MSMTS_PER_NODE"] = 1
GLOBALDICT["MODELDESIGN"]["MULTIPLER_R_MAX"] = 4.

repts = 0
particleconfigs = [ [3,2], [9,6], [15,10], [21,14], [30, 20]]

prefix = '_padua_ord_'+str(padua_order)+'_'
lambda_paris_2 = np.load('lambda_pairs_2.npz')
random_variances = np.load('random_variances.npz')

IDX1_SHP = len(random_variances['g1var'])
IDX2_SHP = 30 # len(lambda_paris_2['lambda_1']) not doing all 250
IDX3_SHP = len(particleconfigs)

########################
# Run Script
########################

datatype=['Uni_R', 'Trunc_R']

error_scaling_matrix_0 = np.zeros((len(particleconfigs), max_iterations))
error_scaling_matrix_1 = np.zeros((len(particleconfigs), max_iterations))
uniform_tuning_params= 100*np.ones((IDX1_SHP, IDX2_SHP, IDX3_SHP))
truncg_tuning_params = 100*np.ones((IDX1_SHP, IDX2_SHP, IDX3_SHP))

for idx_1 in range(IDX1_SHP):
    for idx_2 in range(IDX2_SHP):
        for idx_pconfig in range(IDX3_SHP):
            
            try:

                test_data_0 = np.load(path+datatype[0]+prefix+'rand_'+str(idx_1)+'_'+str(idx_2)+'_'+str(idx_pconfig)+'.npz')
                error_scaling_matrix_0[idx_pconfig, :] = np.sum(np.mean(test_data_0['absolute_errors_matrix'], axis=0), axis=1)

                uniform_tuning_params[idx_1, idx_2, idx_pconfig] = error_scaling_matrix_0[idx_pconfig, max_iterations-1]

                test_data_1 = np.load(path+datatype[1]+prefix+'rand_'+str(idx_1)+'_'+str(idx_2)+'_'+str(idx_pconfig)+'.npz')
                error_scaling_matrix_1[idx_pconfig, :] = np.sum(np.mean(test_data_1['absolute_errors_matrix'], axis=0), axis=1)

                truncg_tuning_params[idx_1, idx_2, idx_pconfig] = error_scaling_matrix_1[idx_pconfig, max_iterations-1]
                
            except:
                continue
                
ordered_error_vals_T = np.argsort(truncg_tuning_params.flatten())
min_error_params_T = [np.unravel_index(ordered_error_vals_T[idx], truncg_tuning_params.shape) for idx in range(len(ordered_error_vals_T)) ] 
ordered_error_vals_U = np.argsort(uniform_tuning_params.flatten())
min_error_params_U = [np.unravel_index(ordered_error_vals_U[idx], uniform_tuning_params.shape) for idx in range(len(ordered_error_vals_U)) ] 

np.savez(path+prefix+'opt_params_list', min_error_params_T=min_error_params_T, min_error_params_U=min_error_params_U)


##########################
# Print to Standard Output
##########################

opt_trunc = np.unravel_index(np.argmin(truncg_tuning_params), truncg_tuning_params.shape)
opt_uniform = np.unravel_index(np.argmin(uniform_tuning_params), uniform_tuning_params.shape)
print "Optimal trunc termination error:", truncg_tuning_params[opt_trunc], "Params:", opt_trunc
print "Optimal uniform termination error:", uniform_tuning_params[opt_uniform], "Params:", opt_uniform
print "--------"
print "--------"
print "--------"
print 'First 10 minimal parameters for Uniform sampling'
print min_error_params_U[0:10]
print "--------"
print "--------"
print "--------"
print 'First 10 minimal parameters for Trunc Gaussian sampling'
print min_error_params_T[0:10]
print "--------"
print "--------"
print "-- end ---"

