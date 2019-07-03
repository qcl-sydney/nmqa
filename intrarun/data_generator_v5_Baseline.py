import sys
import os
import copy
import traceback
import numpy as np

########################
# Find qslam modules
########################
sys.path.append('../qslam')

from singlerun import SingleRunAnalysis
from riskanalysis import EngineeredTruth
from qslamdesignparams import GLOBALDICT

########################
# Taking in bash parameters
########################
idx_prefix = int(sys.argv[1])
# truth flag type (five options)
# 0 - Linear 1D field, 25 qubit array
# 1 - Square 2D field, 25 qubit array
# 2 - Gaussian 2D array, 25 qubit array
# 3 - Square 2D field, 9 qubit array
# 4 - Square 2D field, 16 qubit array
max_iterations = int(sys.argv[2])

if max_iterations == 75:
    OPTIMAL_LAMBDA = {}
    OPTIMAL_LAMBDA["max_iter"] = 75
    OPTIMAL_LAMBDA[0]={}
    OPTIMAL_LAMBDA[0]["LAMBDA_1"] = 0.92804888
    OPTIMAL_LAMBDA[0]["LAMBDA_2"] = 0.68006484
    OPTIMAL_LAMBDA[1]={}
    OPTIMAL_LAMBDA[1]["LAMBDA_1"] = 0.92804888
    OPTIMAL_LAMBDA[1]["LAMBDA_2"] = 0.68006484 
    OPTIMAL_LAMBDA[2]={}
    OPTIMAL_LAMBDA[2]["LAMBDA_1"] = 0.81146779
    OPTIMAL_LAMBDA[2]["LAMBDA_2"] = 0.62820286 

########################
# Truth Parameters
########################

# Choose defaults to match floor case (heights didn't work)
TRUTHKWARGS = {}

BARRIER_FLOOR = 0.25*np.pi
BARRIER_HEIGHT = 0.75*np.pi
FLOOR_RATIO = 0.4 # matches floor case

TRUTHKWARGS["OneStepdheight"] = {"low": BARRIER_FLOOR, 
                                 "high": BARRIER_HEIGHT} 
TRUTHKWARGS["OneStepdfloorarea"] = FLOOR_RATIO 


########################
# Save to path 
########################

path = './data_v3/'

########################
# Set true field
########################

if idx_prefix == 0:
    change_gridconfig = True # 1D
    TRUTHFLAG = None # use TRUTHKWARGS
    TRUTHKWARGS["truthtype"] = 'OneStepd' 

if idx_prefix == 1: 
    change_gridconfig = False # 2D
    TRUTHFLAG = None # use TRUTHKWARGS
    TRUTHKWARGS["truthtype"] = 'OneStepq' 

if idx_prefix == 2:
    change_gridconfig = False # 2D
    TRUTHFLAG = None # use TRUTHKWARGS
    TRUTHKWARGS["truthtype"] = 'Gaussian' 
    
if idx_prefix > 2: # Covers all truth flag options except 2.
    change_gridconfig = False # 2D
    TRUTHFLAG = None # use TRUTHKWARGS
    TRUTHKWARGS["truthtype"] = 'OneStepq' 

########################
# Set 1D Hardware if req
########################

num_of_nodes = len(GLOBALDICT["GRIDDICT"])

if change_gridconfig is True:
    
    GLOBALDICT["GRIDDICT"] = {}
    
    for idx_posy in range(num_of_nodes):
        
        if idx_posy < 10 :
            
            GLOBALDICT["GRIDDICT"]["QUBIT_0" + str(idx_posy)] = (0.0, float(idx_posy))
            
        if idx_posy >= 10 :
            
            GLOBALDICT["GRIDDICT"]["QUBIT_" + str(idx_posy)] = (0.0, float(idx_posy))

########################
# Change 2D Grid Size if req
########################

if idx_prefix == 3:
    from qslamdesignparams import GRIDDICT_9
    GLOBALDICT["GRIDDICT"] = copy.deepcopy(GRIDDICT_9)
    num_of_nodes = len(GLOBALDICT["GRIDDICT"])

if idx_prefix == 4:
    from qslamdesignparams import GRIDDICT_16
    GLOBALDICT["GRIDDICT"] = copy.deepcopy(GRIDDICT_16)
    num_of_nodes = len(GLOBALDICT["GRIDDICT"])

########################
# Set Baseline Defaults
########################

change_MSMTS_PER_NODE = 1
change_SIGMOID_APPROX_ERROR = 10.0**(-6)
change_QUANTISATION_UNCERTY = 10.0**(-4)
change_P_ALPHA = 15 
change_P_BETA = 10 


GLOBALDICT["MODELDESIGN"]["MAX_NUM_ITERATIONS"] = OPTIMAL_LAMBDA["max_iter"]
GLOBALDICT["MODELDESIGN"]["MSMTS_PER_NODE"] = change_MSMTS_PER_NODE
GLOBALDICT["NOISEPARAMS"]["SIGMOID_APPROX_ERROR"]["SIGMA"] = change_SIGMOID_APPROX_ERROR
GLOBALDICT["NOISEPARAMS"]["QUANTISATION_UNCERTY"]["SIGMA"] = change_QUANTISATION_UNCERTY
GLOBALDICT["MODELDESIGN"]["P_ALPHA"] = change_P_ALPHA
GLOBALDICT["MODELDESIGN"]["P_BETA"] = change_P_BETA
GLOBALDICT["MODELDESIGN"]["LAMBDA_1"] = OPTIMAL_LAMBDA[idx_prefix]["LAMBDA_1"]
GLOBALDICT["MODELDESIGN"]["LAMBDA_2"] = OPTIMAL_LAMBDA[idx_prefix]["LAMBDA_2"]

repts = 50
prefix = '_idx_prefix_'+str(idx_prefix)+'_'

########################
# Run Script
########################

truemap_generator = EngineeredTruth(num_of_nodes, TRUTHKWARGS)
true_map_ = truemap_generator.get_map()

fname_likelihood = 'baseline_'

uniform_r_expt = SingleRunAnalysis(GLOBALDICT, true_map_, repts, beta_expansion_mode=False, beta_skew_adjust=False)
uniform_r_expt.run_analysis(path+'Uni_R'+prefix+fname_likelihood)

trunc_r_expt = SingleRunAnalysis(GLOBALDICT, true_map_, repts, beta_expansion_mode=True, beta_skew_adjust=False)
trunc_r_expt.run_analysis(path+'Trunc_R'+prefix+fname_likelihood)

