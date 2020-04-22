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

# RANDOM POLYNOMIAL FUNCTIONALITY
from randpoly import NEARESTPADUA

sys.path.append('./')


########################
# Taking in bash parameters
########################
padua_order = int(sys.argv[1])
idx_functype = int(sys.argv[2])

########################
# Set True Field
########################

if idx_functype ==5:
    true_function_type = 'randpoly'

if idx_functype ==6:
    true_function_type = 'randpolymax'    

data_qubit_num = 25
data_qubit_flag ='uniform'

########################
# Save to path 
########################
path = './data/'


if true_function_type != 'randpolymax':

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

    if padua_order == -7:
        # specific grid to address random polynomial functionality
        sensing_qubits = generate_data_qubits_coords(25, flag=data_qubit_flag)
        sensing_qubits = list(np.asarray(sensing_qubits)*0.75)

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
        padua_order = "no_padua"

    if padua_order == -2:
        data_qubits = generate_data_qubits_coords(data_qubit_num, flag=data_qubit_flag)
        # remove duplicate sensors:
        sensing_qubits = list(set(sensing_qubits) - set(data_qubits))

        # update dictionary params:
        GLOBALDICT["DATA_QUBITS"] = np.arange(len(sensing_qubits), len(sensing_qubits) + data_qubit_num, dtype='int')
        GLOBALDICT["INTERPOLATE_FLAG"] = 'Rbf' 
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

    if padua_order == -7:
        prefix = true_function_type +'_reg25_'
        padua_order = "reg25"


if true_function_type == 'randpolymax':
    
    REGGRID=81
    
    data_qubits = generate_data_qubits_coords(data_qubit_num, flag=data_qubit_flag)
    
    if REGGRID == 16:
        sensing_qubits = generate_data_qubits_coords(REGGRID, flag=data_qubit_flag)
        sensing_qubits = list(np.asarray(sensing_qubits) * 0.75)
    
    if REGGRID == 81:
        sensing_qubits = generate_data_qubits_coords(REGGRID, flag=data_qubit_flag)        
        sensing_qubits = list(set(sensing_qubits) - set(data_qubits))

    # update dictionary params:
    GLOBALDICT["DATA_QUBITS"] = np.arange(len(sensing_qubits), len(sensing_qubits) + data_qubit_num, dtype='int')
    GLOBALDICT["INTERPOLATE_FLAG"] = 'Rbf' 
    prefix = true_function_type +'_reg_'+str(REGGRID)+'_'+'ord_'+str(padua_order)+'_' 
    
    
    
########################
# Set hardware and true map
########################
TRUTHKWARGS = {}
TRUTHKWARGS["truthtype"] = "UseFunction"
TRUTHKWARGS["true_function"] = true_function
TRUTHKWARGS["true_function_type"] = 'randpoly'

# PADUA PROJECT: RANDOM POLYNOMIALS FUNCTIONALITY
if true_function_type == 'randpoly': 
    TRUTHKWARGS["randpoly"]={}
    if padua_order > 0:
        TRUTHKWARGS["randpoly"]["n"]= padua_order 
    if isinstance(padua_order, str):
        TRUTHKWARGS["randpoly"]["n"] = NEARESTPADUA[str(padua_order)]    
    TRUTHKWARGS["randpoly"]["trial"]= None  # New polynomial each time .get_map() is called


if true_function_type == 'randpolymax':
    TRUTHKWARGS["randpoly"]={}
    TRUTHKWARGS["randpoly"]["n"]= padua_order  
    TRUTHKWARGS["randpoly"]["trial"]= None  # New polynomial each time .get_map() is called


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
# Set Simulation Params
########################

GLOBALDICT["MODELDESIGN"]["MSMTS_PER_NODE"] = 1
GLOBALDICT["MODELDESIGN"]["MULTIPLER_R_MAX"] = 4.
repts = 50

########################
# Set Loop Parameters
########################

Multiples = [1, 3, 5, 7, 9, 10, 15, 20, 50]

########################
# Run Script
######## ################

# RANDOM POLYNOMIAL FUNCTIONALITY
fname_likelihood = 'interpol_only' 

for idx_msmt_var in Multiples:

    unique_id = path + prefix + fname_likelihood + '_m_' + str(idx_msmt_var)
    GLOBALDICT["MODELDESIGN"]["MAX_NUM_ITERATIONS"] = len(sensing_qubits) * idx_msmt_var
    GLOBALDICT["MODELDESIGN"]["ID"] = unique_id
    naive_br = 0.
    naivedata = 0.
    
    try:
        naive_br = risknaive(copy.deepcopy(TRUTHKWARGS), copy.deepcopy(GLOBALDICT))
        naive_br.naive_implementation()
        

    except:
        print "Index: %s was not completed..." %(idx_msmt_var)
        print "Error information:"
        print "Type", sys.exc_info()[0]
        print "Value", sys.exc_info()[1]
        print "Traceback", traceback.format_exc()
