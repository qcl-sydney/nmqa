import sys
import os
import copy
import traceback
import numpy as np

########################
# Find qslam modules
########################
sys.path.append('../qslam')

from qslamdesignparams import GLOBALDICT
from exptrisk import EmpiricalRisk
from visualiserisk import DataCube, cm2inch
from experimentaldata import DataKeys

########################
# Taking in bash parameters /
########################
datakey = int(sys.argv[1]) # data key
num_qubits = DataKeys[str(datakey)]['parameters']['N']
idx_job_array = int(sys.argv[2]) # job array starts from 1

########################
# Save to path 
########################
savetopath = '/scratch/QCL_RG/qslamdatapaper_v3/' # on Artemis
# savetopath = './' # local

########################
# Set 1D Hardware  to Ion Trap
########################

change_gridconfig = True

# assume equi-distant linear array

if change_gridconfig is True:

    GLOBALDICT["GRIDDICT"] = {}
    
    for idx_posy in range(num_qubits):
        if idx_posy < 10 :
            GLOBALDICT["GRIDDICT"]["QUBIT_0" + str(idx_posy)] = (0.0, float(idx_posy))
        if idx_posy >= 10 :
            GLOBALDICT["GRIDDICT"]["QUBIT_" + str(idx_posy)] = (0.0, float(idx_posy))

########################
# Set Defaults
########################

GLOBALDICT["MODELDESIGN"]["MSMTS_PER_NODE"] = 1
GLOBALDICT["MODELDESIGN"]["P_ALPHA"] = 15
GLOBALDICT["MODELDESIGN"]["P_BETA"] = 10


########################
# Set All Loop Parameters
########################

msmt_per_qubit_scan = [1]
meta_max_iter_scan = [ 5, 10, 15, 20, 25, 50, 75, 100, 125, 250]
lambda_paris_2 = np.load('lambda_pairs_2.npz')
random_variances = np.load('random_variances.npz')

IDX1_SHP = len(random_variances['g1var'])
IDX2_SHP = 30 # len(lambda_paris_2['lambda_1']) not doing all 250
              
########################
# Run Script
########################

idx_1, idx_2 = np.unravel_index(idx_job_array - 1 , (IDX1_SHP, IDX2_SHP) )

GLOBALDICT["NOISEPARAMS"]["SIGMOID_APPROX_ERROR"]["SIGMA"] = random_variances['g2var'][idx_1]
GLOBALDICT["NOISEPARAMS"]["QUANTISATION_UNCERTY"]["SIGMA"] = random_variances['g1var'][idx_1]
GLOBALDICT["MODELDESIGN"]["LAMBDA_1"] = lambda_paris_2['lambda_1'][idx_2]
GLOBALDICT["MODELDESIGN"]["LAMBDA_2"] = lambda_paris_2['lambda_2'][idx_2]

param_index = str(idx_1)+'_'+str(idx_2)


meta_ssim_pairs_1 = []
meta_empr_pairs_1 = [] 
ssim_qslam = []
err_qslam = []
err_naive = []
ssim_naive =[]
    
for idx_msmt_iter in range(len(meta_max_iter_scan)):

    GLOBALDICT["MODELDESIGN"]["MAX_NUM_ITERATIONS"] = meta_max_iter_scan[idx_msmt_iter]
    
    expt = EmpiricalRisk(GLOBALDICT, datakey)
    err, ssim = expt.calculate_risk(number_of_trials=50)
    
    ssim_qslam.append(ssim[0])
    err_qslam.append(err[0])
    ssim_naive.append(ssim[1])
    err_naive.append(err[1])

meta_ssim_pairs_1.append([ssim_qslam, ssim_naive])
meta_empr_pairs_1.append([err_qslam, err_naive])

np.savez('2019_Jun_qslam_exptdata_'+str(datakey)+'_param_'+param_index, 
         meta_ssim_pairs=meta_ssim_pairs_1, 
         meta_empr_pairs=meta_empr_pairs_1)

