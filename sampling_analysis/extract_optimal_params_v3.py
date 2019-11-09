import sys
import os
import copy
import traceback
import numpy as np

sys.path.append('../paduaq')
from pdpoints import dims_padua_set

########################
# Taking in bash parameters
########################
padua_order = int(sys.argv[1]) # Padua order. 1, 2, 3, 4, 5,...

idx_functype = int(sys.argv[2])
if idx_functype ==0:
    true_function_type = 'cheb2fun'
    MULTIPLIER=5
if idx_functype ==1:
    true_function_type = 'lin'
    MULTIPLIER=5
    
data_qubit_num = 25

# For each type of grid, data is analysed when the num of msmts per sensor = MULTIPLIER

if padua_order >0:
    max_iterations = int(dims_padua_set(padua_order) * MULTIPLIER) 
    prefix = true_function_type+'_padua_ord_'+str(padua_order)+'_'

if padua_order == -1: 
    ALLGRID = data_qubit_num
    max_iterations = ALLGRID * MULTIPLIER 
    prefix = true_function_type +'_no_padua_'
    
if padua_order == 'REG_COARSE':
    COARSEGRID = 16
    data_qubits = COARSEGRID * MULTIPLIER 
    prefix = true_function_type +'_regcoarse_'
    
if padua_order == 'REG_FINE':
    FINEGRID = 81
    data_qubits = FINEGRID * MULTIPLIER 
    prefix = true_function_type +'_regfine_'
    
########################
# Save to path 
########################

path = '/scratch/QCL_RG/qslam_padua_paper/' # on Artemis './data/'

########################
# Set Defaults
########################

particleconfigs = [ [3,2], [9,6], [15,10], [21,14], [30, 20]]

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
                print "Type", sys.exc_info()[0]
                print "Value", sys.exc_info()[1]
                print "Traceback", traceback.format_exc()
                
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

