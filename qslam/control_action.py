'''
Created on Thu Apr 20 19:20:43 2017
@author: riddhisw

.. module:: control_action

    :synopsis: Implements controller specifying measurement locations for qslam.

    Module Level Functions:
    ----------------------
        control_lengthscale_uncertainty : Return location(s) for next set of
            single qubit measurement(s) selecting qubit locations with highest state
            estimate uncertainty.

        control_user_input : Return location(s) for next set oF single qubit
            measurement(s) as defined by the user.

.. moduleauthor:: Riddhi Gupta <riddhi.sw@gmail.com>
'''

import numpy as np

# == VERSION OCT 2019 ==========================================================
# Fix: data qubits are never measured. 

def control_lengthscale_uncertainty2(listofcontrolparameters,
                                    next_control_neighbourhood,
                                    number_of_diff_nodes=1,
                                    list_of_dataqubits=None,
                                    dtype=[('Node', int), ('ControlParam', float)]
                                   ):

    ''' Return location(s) for next set of single qubit measurement(s)
        selecting qubit locations with highest state estimate uncertainty.

    Parameters:
    ----------
        listofcontrolparameters (`float64`| numpy array):
            Avg. uncertainity metric for correlation lengthscales at each node.
            Dims: Grid.number_of_nodes

        next_control_neighbourhood (`int`| list) :
            List of indices for qubit locations, denoted a control region.
            The highest uncertainity location will be chosen from this set
            and be selected for the next physical measurement.

        number_of_diff_nodes (`int`| scalar | optional):
            Number of single qubit measurements at different locations
            that can be simultaneously performed on the hardware grid.

        list_of_dataqubits (`int`| list) :
            List of indices for qubit locations for data qubits that
            cannot be used for sensing measurements.
            Defaults to None.

        dtype ( List of tuples | optional) :
            Specifies how control parameters and control neighbourhoods are read,
            and handled by this function. Should be a hidden local variable.

    Returns:
    -------
        controls_list (`float64` | numpy array):
            Location(s) for performnng the next single qubit measurement.
            Dims: number_of_diff_nodes
    '''
    
    # structured_array labels all sensor and data qubit locations
    labelled_params = [iterate for iterate in enumerate(listofcontrolparameters)]
    structured_array = np.asarray(labelled_params, dtype=dtype) 
    grid_size = len(labelled_params)
    
    # filter out qubits not in scope (e.g. data qubits)
    mask = generate_control_mask(grid_size, next_control_neighbourhood, list_of_dataqubits)
    structured_array = structured_array[mask]
    
    # get location(s) of next physical measurement
    controls_list = extract_max_var_region(structured_array, number_of_diff_nodes)
    
    return controls_list

def extract_max_var_region(structured_array, number_of_diff_nodes):
    '''Helper function to control_lengthscale_uncertainty2.
    
    Return the set of qubit indicies on which physical measurements are scheduled.
    
    structured_array ('structured numpy array' | dtype=[('Node', int), ('ControlParam', float)):
        
        Structured array with uncertainty metric given by 'ControlParam'
            and qubit location given by 'None'. 
        number_of_diff_nodes number of Node indices with highest variance metrics
            are selected. Selection of qubits is uniformly randomly sampled 
            if variance metrics are the same for a number of qubits greater than 
            number_of_diff_nodes. 
    
    number_of_diff_nodes (`int`| scalar ):
        Number of single qubit measurements at different locations
        that can be simultaneously performed on the hardware grid.
    '''
    
    # sort the region from max to min uncertainty.
    sorted_array = np.sort(structured_array, order=['ControlParam', 'Node'])[::-1]

    # Zeroth term of sorted_array is the node corresponding to maximal ucertainity.
    max_val = sorted_array['ControlParam'][0]

    # Now find multiple instances of the maximum value if they exist.
    # If the counter > 1, then multiple maxima exist.
    counter = 1
    for controlparamval in sorted_array['ControlParam'][1:]:
        if controlparamval == max_val:
            counter += 1
        elif controlparamval != max_val:
            break
    # Return controls based on highest uncertainty.
    # If equi-certain options exist, choose between these options at random.
    if number_of_diff_nodes <= counter:
        # Returns random choices among equicertain nodes.
        chosen_node_indices = np.random.randint(low=0, high=counter, size=number_of_diff_nodes)
        return sorted_array['Node'][chosen_node_indices]
 
    elif number_of_diff_nodes > counter:
        # Return nodes in descending order of uncertainty
        return sorted_array['Node'][0 : number_of_diff_nodes]


def generate_control_mask(grid_size, next_control_neighbourhood, list_of_dataqubits):
    ''' Helper function to control_lengthscale_uncertainty2.
    
    Return boolean mask for hiding qubit location indicies not in scope for the controller.
    
    Notes:
        Qubits are filtered by location indices. 
        Physical coordinates for each qubit with different location indices are not 
            checked for degeneracy (e.g. two labels refer to same qubit). 
        Mask is assigned True if sensor qubit, and False if data qubit. 
        If a qubit is (mistakenly) assigned both as a sensor and as a data qubit, 
            then mask will be False.
        
    grid_size : Total number of sensor qubits and data qubits.
    
    next_control_neighbourhood: List of qubit indicies (spanning a subset of the grid or the whole grid) 
        from which the controller decides the next physical measurement. 
        If empty, all sensor qubits on the grid are considered in scope. 
        Mask elements are set to True.
        
    list_of_dataqubits: List of data qubit indicies (spanning a subset of the grid) 
        which cannot be physically measured. 
        Mask elements are set to False.
    
    '''

    if len(next_control_neighbourhood) > 0:
        mask = np.zeros(grid_size, dtype=bool) # Mask hides all qubits on grid
        mask[next_control_neighbourhood] = True # Shows qubits only within next_control_neighbourhood.
        if list_of_dataqubits is not None:
            mask[list_of_dataqubits] = False # Hide data qubits

    if len(next_control_neighbourhood) == 0 or np.sum(mask) == 0 :     
        mask = np.ones(grid_size, dtype=bool) # Mask shows all qubits on grid
        if list_of_dataqubits is not None:
            mask[list_of_dataqubits] = False # Hide data qubits
            
        if np.sum(mask) == 0:
            print("No sensor qubits specified; raise RunTime Error")
            raise RuntimeError
            
    return mask

# == OLD VERSION ===============================================================

def control_lengthscale_uncertainty(listofcontrolparameters,
                                    next_control_neighbourhood,
                                    number_of_diff_nodes=1,
                                    list_of_dataqubits=None,
                                    dtype=[('Node', int), ('ControlParam', float)]
                                   ):

    ''' Return location(s) for next set of single qubit measurement(s)
        selecting qubit locations with highest state estimate uncertainty.

    Parameters:
    ----------
        listofcontrolparameters (`float64`| numpy array):
            Avg. uncertainity metric for correlation lengthscales at each node.
            Dims: Grid.number_of_nodes

        next_control_neighbourhood (`int`| list) :
            List of indices for qubit locations, denoted a control region.
            The highest uncertainity location will be chosen from this set
            and be selected for the next physical measurement.

        number_of_diff_nodes (`int`| scalar | optional):
            Number of single qubit measurements at different locations
            that can be simultaneously performed on the hardware grid.

        list_of_dataqubits (`int`| list) :
            List of indices for qubit locations for data qubits that
            cannot be used for sensing measurements.
            Defaults to None.

        dtype ( List of tuples | optional) :
            Specifies how control parameters and control neighbourhoods are read,
            and handled by this function. Should be a hidden local variable.

    Returns:
    -------
        controls_list (`float64` | numpy array):
            Location(s) for performnng the next single qubit measurement.
            Dims: number_of_diff_nodes
    '''

    # COMMENT: store control parameters in a structureed numpy array.
    # Then, mask control parameters which are corresspond to be out of the
    # control neighbourhood. If the control is empty, include all qubits
    # in the analysis for choosing the next measurement.

    labelled_params = [iterate for iterate in enumerate(listofcontrolparameters)]
    structured_array = np.asarray(labelled_params, dtype=dtype) 
    
    if len(next_control_neighbourhood) > 0:
        mask = np.zeros(len(labelled_params), dtype=bool) # Mask for hiding all values.
        mask[next_control_neighbourhood] = True # Show nodes only from next_control_neighbourhood.
        
        if list_of_dataqubits is not None:
            mask[list_of_dataqubits] = False # Show nodes only from next_control_neighbourhood.

    if len(next_control_neighbourhood) == 0 or np.sum(mask) == 0 : 
        mask = np.ones(len(labelled_params), dtype=bool) # Mask for showing all values.
        # print "Control List empty; randomly chosen qubit on grid."

    structured_array = structured_array[mask] # No mask for lists - make array and bring back.

    # Now, sort the region from max to min uncertainty.
    sorted_array = np.sort(structured_array, order=['ControlParam', 'Node'])[::-1]

    # Zeroth term of sorted_array is the node corresponding to maximal ucertainity.
    max_val = sorted_array['ControlParam'][0]

    # Now find multiple instances of the maximum value if they exist.
    # If the counter > 1, then multiple maxima exist.
    counter = 1
    for controlparamval in sorted_array['ControlParam'][1:]:
        if controlparamval == max_val:
            counter += 1
        elif controlparamval != max_val:
            break

    # Return controls based on highest uncertainty.
    # If equi-certain options exist, choose between these options at random.
    if number_of_diff_nodes < counter:
        # Returns random choices among equicertain nodes.
        chosen_node_indices = np.random.randint(low=0, high=counter, size=number_of_diff_nodes)
        return sorted_array['Node'][chosen_node_indices]

    elif number_of_diff_nodes >= counter:
        # Return nodes in descending order of uncertainty
        return sorted_array['Node'][0 : number_of_diff_nodes]


def control_user_input(input_controls, next_control_neighbourhood):
    ''' Return location(s) for next set of single qubit measurement(s). No control
    protocol specified. '''

    # do nothing

    return input_controls


PROTOCOL = {"userinput" : control_user_input,
            "lenvar": control_lengthscale_uncertainty2
           }


def controller(listofcontrolparameters, next_control_neighbourhood, 
               controltype='lenvar',
               list_of_dataqubits=None,
               number_of_nodes=1):
    ''' Return location(s) for next set of single qubit measurement(s) based on
        a selected control protocol.

        Parameters:
        ----------
        listofcontrolparameters (`float64`| numpy array-like):
            List of input information for the control protocol specified by controltype.

        controltype (`str`| optional):
            Specifies control protocol to be used:
                'lenvar' : Calls protocol control_lengthscale_uncertainty.
                'userinput' : Returns listofcontrolparameters (no control).
            Defaults to 'lenvar'.

        number_of_nodes (`int`| scalar | optional):
            Number of single qubit measurements that can be simultaneously
            performend on the hardware grid.
            Defaults to 1.

        list_of_dataqubits (`int`| list) :
            List of indices for qubit locations for data qubits that
            cannot be used for sensing measurements.
            Defaults to None.

        Returns:
        -------
        controls_list (`float64` | numpy array):
            Location(s) for performnng the next single qubit measurement.
            Dims: number_of_nodes
    '''

    if len(listofcontrolparameters) > 0:

        controls_list = PROTOCOL[controltype](listofcontrolparameters, next_control_neighbourhood, number_of_nodes)
        return controls_list[0:number_of_nodes]

    elif len(listofcontrolparameters) == 0:

        print "Input control parameters empty."
        print "Runtime Error raised"
        raise RuntimeError
