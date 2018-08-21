'''
Created on Thu Apr 20 19:20:43 2017
@author: riddhisw

.. module:: hardware

    :synopsis: Incorporates a collection of physical qubits on hardware,
        tracking both measurements and state information for each qubit.

    Module Level Classes:
    ----------------------
        Node : Maintains spatial coordinates, state estimates and measurement data
            for a single qubit location.
        Grid : Maintains spatial coordinates, state estimates and measurement data
            for an arbirary arrangement of qubits.

.. moduleauthor:: Riddhi Gupta <riddhi.sw@gmail.com>
'''
import numpy as np
from model_design import INITIALDICT

LAMBDA = INITIALDICT["LAMBDA"]
PARTICLE_STATE = ["x_state", "y_state", "f_state", "r_state"]


###############################################################################
# CHIP STRUCTURE
###############################################################################
class Node(object):
    ''' Maintains spatial coordinates, state estimates and measurement data
        for a single qubit location.

    Attributes:
    ----------
        r_state (`float64` | scalar ):
            State variable for the correlation lengthscale about this location.
        x_state (`float64` | scalar ):
            State variable for x-coordinate of location of the qubit.
        y_state (`float64` | scalar ):
            State variable for y-coordinate of location of the qubit.
        counter_tau (`int` | scalar ):
            Counter for the number of times a physical measurement is
                made at this location.
        counter_beta (`int` | scalar ):
            Counter for the number of times a quasi measurement is
                made at this location.
        lambda_factor (`int` | scalar ):
            Hyper-parameter governing the rate at which quasi-measuretments
                are forgotten as physical measurements accumulate.

    Properties:
    ----------
        f_state (`float64` | scalar | non-public attribute):
            State variable for the dephasing noise field about this location
                caculated based on physical measurements and inferred data.
                Initial Value: 0.0
                Setter Function: None.
        r_state_variance (`float64` | scalar ):
            Derived quantity from state variables tracking the level of uncertainty
                in qslam's knowledge of correlation lengths at this location.
                Initial Value: 10**4
                Setter Function: r_state_variance(var_metric)
        physcmsmtsum (`float64` | scalar ):
            Sum of binary physical measurement outcomes.
                Initial Value: 0.0
                Setter Function: physcmsmtsum(next_phys_msmt)
        quasimsmtsum (`float64` | scalar ):
            Sum of binary quasi measurement outcomes.
                Initial Value: 0.0
                Setter Function: quasimsmtsum(next_quasi_msmt)

    Methods:
    -------
        sample_prob_from_msmts:
            Returns a sample probability for observing the qubit up-state based
                on physical and quasi-measurement attributes, mediated by lambda
                forgetting factor.

    Static Methods:
    --------------
        born_rule:
            Return Born probability given qubit phase for a Ramsey experiment
                under dephasing noise.
        inverse_born:
            Return qubit phase estimate given Born probability estimate for a
                Ramsey experiment under dephasing noise.

    '''

    def __init__(self):

        self._f_state = np.random.uniform(low=0.0, high=np.pi) # cant set _f_state
        self.r_state = 0.0
        self.__r_state_variance = 10**4 # COMMENT: default value. Very large to make threshold.
        self.x_state = 0.0
        self.y_state = 0.0
        self.counter_tau = 0
        self.counter_beta = 0
        self.__physcmsmtsum = 0.0
        self.__quasimsmtsum = 0.0
        self.lambda_factor = LAMBDA

    @property
    def physcmsmtsum(self):
        '''Return sum of binary physical measurement outcomes.'''
        return self.__physcmsmtsum
    @physcmsmtsum.setter
    def physcmsmtsum(self, next_phys_msmt):
        '''Update Node attributes with a new physical measurement.'''
        self.__physcmsmtsum += next_phys_msmt
        self.counter_tau += 1

    @property
    def quasimsmtsum(self):
        '''Return sum of binary quasi measurement outcomes.'''
        return self.__quasimsmtsum
    @quasimsmtsum.setter
    def quasimsmtsum(self, next_quasi_msmt):
        '''Update Node attributes with a new quasi-measurement.'''
        self.__quasimsmtsum += next_quasi_msmt
        self.counter_beta += 1


    ############################################################################
    # PLACEHOLDER kew based control for SLAM
    ############################################################################
    @property
    def r_state_variance(self):
        ''' Return level of uncertainty in qslam's knowledge of correlation lengths
            at this location. '''
        return self.__r_state_variance
    @r_state_variance.setter
    def r_state_variance(self, skew_metric):
        ''' Update r_state_variance with new value.'''
        self.__r_state_variance = skew_metric


    ############################################################################
    @property
    def f_state(self): # no .setter function
        '''Return state variable for the dephasing noise field about this location.'''
        prob_sample = self.sample_prob_from_msmts()
        if prob_sample is None:
            # self._f_state = np.random.uniform(low=0.0, high=np.pi)
            # print "Map value returned is default; prob_sample NONE in f_state", self._f_state
            return self._f_state
        if prob_sample >= 0.0 and prob_sample <= 1.0:
            self._f_state = Node.inverse_born(prob_sample)
            return self._f_state
        # print "INVALID prob_sample value encountered in calling f_state", prob_sample
        print prob_sample
        raise RuntimeError

    def sample_prob_from_msmts(self): # TODO Data Association
        ''' Returns a sample probability for observing the qubit up-state based
            on physical and quasi-measurements.'''

        forgetting_factor = self.lambda_factor**self.counter_tau

        w_q = None
        prob_q = 0.0
        if self.counter_beta != 0:
            prob_q = self.quasimsmtsum / self.counter_beta*1.0
            w_q = 0.5*forgetting_factor

        w_p = None
        prob_p = 0.0
        if self.counter_tau != 0:
            prob_p = self.physcmsmtsum / self.counter_tau*1.0
            w_p = 0.5 + 0.5*(1 - forgetting_factor)

        if w_p is None and w_q is None:
            # print "NONE returned in sample_prob_from_msmts"
            return  None
        elif w_p is not None  and w_q is None:
            w_p = 1.0
            w_q = 0.0
        elif w_p is None and w_q is not None:
            w_q = 1.0
            w_p = 0.0
        elif w_p is not None and w_q is not None:
            pass
        prob_j = w_p*prob_p + w_q*prob_q

        if prob_j > 1 or prob_j < 0:
            raise RuntimeError
        return prob_j


    @staticmethod
    def born_rule(map_val):
        '''Return Born probability given qubit phase for a Ramsey experiment
            under dephasing noise. '''
        born_prob = np.cos(map_val / 2.0)**2
        return born_prob


    @staticmethod
    def inverse_born(born_prob):
        '''Return qubit phase estimate given Born probability estimate
             for a Ramsey experiment under dephasing noise.'''
        map_val = np.arccos(2.0*born_prob  - 1.0)
        return map_val

class Grid(object):
    ''' Maintains spatial coordinates, state estimates and measurement data
    for an arbirary arrangement of qubits.

    Attributes:
    ----------
        list_of_nodes_positions (`float64` | list ):
            List of (x, y) pairs for defining spatial coordinates for qubits on hardware.
        engineeredtruemap (`dtype` | scalar ):
            Map of true (engineered) dephasing noise field for sumulation analysis.

        number_of_nodes (`int` | scalar ):
            Total number of qubits on hardware.
        nodes (`Node` class object | list ):
            List of `Node` objects, where one `Node` object exists for each qubit.
        control_sequence (`dtype` | list ):
            List of control actions (qubits to measure) for every iteration step of qslam.
            Initalised as an empty list.

    Properties:
    ----------
        state_vector (`float64` | numpy array ):
            State variable representing all state variables for all qubit nodes in Grid.
                Initial Value: 0.0
                Setter Function: state_vector(new_state_vector).
    Methods:
    -------
        measure_node : Return a 0 or 1 physical measurement for performing a
            single shot Ramsey experiment on a qubit.

    '''

    def __init__(self, list_of_nodes_positions=None, engineeredtruemap=None):
        if list_of_nodes_positions is None:
            print "No node positions specified"
            raise RuntimeError

        self.list_of_nodes_positions = list_of_nodes_positions
        self.number_of_nodes = len(self.list_of_nodes_positions)
        self.nodes = [Node() for i in range(self.number_of_nodes)]
        self.__state_vector = 0.0
        self.state_vector = np.zeros(self.number_of_nodes*len(PARTICLE_STATE))

        self.engineeredtruemap = engineeredtruemap
        self.control_sequence = []

        for item in xrange(self.number_of_nodes):
            self.nodes[item].x_state, self.nodes[item].y_state = self.list_of_nodes_positions[item]

    def get_all_nodes(self, attribute_list):
        '''Returns attribute for all nodes in a single vector, stacked by nodes, then attributes '''
        vector = [[getattr(node, attr) for node in self.nodes] for attr in attribute_list]
        return np.asarray(vector).flatten()

    def set_all_nodes(self, single_attribute, attribute_values):
        '''Sets attribute for all nodes '''
        for item in xrange(len(attribute_values)):
            setattr(self.nodes[item], single_attribute, attribute_values[item])

    @property
    def state_vector(self):
        ''' Return state variable representing all state variables for all qubit
            nodes in Grid.'''
        self.__state_vector = self.get_all_nodes(PARTICLE_STATE)
        return self.__state_vector
    @state_vector.setter
    def state_vector(self, new_state_vector):
        ''' Set state variable representing all state variables for all qubit nodes
            in Grid.'''
        for state in xrange(len(PARTICLE_STATE)):
            if state == 2: # f_state cannot be set, only calculated
                pass
            elif state != 2:
                bgn = state * self.number_of_nodes
                end = bgn + self.number_of_nodes
                self.set_all_nodes(PARTICLE_STATE[state], new_state_vector[bgn:end])

    def measure_node(self, node_j):
        '''Return a 0 or 1 physical measurement for performing a
            single shot Ramsey experiment on a qubit.'''

        # ENGINEERED NOISE MSMT
        qubit_phase = self.engineeredtruemap[node_j]
        # TODO: Perturb randomly for a real noise field.
        born_prob = Node.born_rule(qubit_phase)
        msmt = np.random.binomial(1, born_prob)

        self.control_sequence.append(node_j)

        return msmt
