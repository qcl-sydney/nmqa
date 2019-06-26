'''
Created on Thu Apr 20 19:20:43 2017
@author: riddhisw

.. module:: exptrisk

    :synopsis: Computes expected value of estimated error for map reconstruction
        using experimental measurement data.

    Module Level Classes:
    ----------------------
        EmpiricalRisk : Computes expected value of estimated map error using QSLAM
            and experimental measurement data.
        NaiveEstimatorExpt : Computes expected value of estimated map error using naive
            brute force measurement strategy and experimental measurement data.

.. moduleauthor:: Riddhi Gupta <riddhi.sw@gmail.com>
'''
import copy
import numpy as np
from visualiserisk import Metric
import qslamr as qs
from experimentaldata import RealData
from hardware import Node


class EmpiricalRisk(object):
    '''Computes expected value of estimated map error using QSLAM
       and experimental measurement data.

    Class Methods:
    --------------
    qslam_trial : Return posterior map, posterior lengthscales and control path
        for one QSLAM run.
    naive_trial : Return recontructed maps using experimental measurement data.
    calculate_risk : Return error matrices for map reconstructions for QSLAM and Naive
        measurement strategies using experimental measurement data.
    '''

    def __init__(self, GLOBALDICT, data_key):

        self.GLOBALDICT = GLOBALDICT
        self.data_key = data_key

    def qslam_trial(self, measurements_controls_=None,
                    autocontrol_="ON",
                    var_thres_=1.0):
        '''
        Return posterior map, posterior lengthscales and control path for one
        QSLAM run.

        Parameters:
        -----------

        measurements_controls_ (`float` | dims: 2):
            External input for the most recent control directive and measurement
            outcome for the qubit, denoted by locaton index, node_j.
            Default value: None.
        autocontrol_ (`dtype` | string binary - "ON", "OFF"):
            "OFF" - next qubit measured is specified as a user input via measurement_controls_
            "ON" - next qubit measured is chosen by the algorithm.
            Default value: "ON".
        var_thres [NOT USED]:
            Error variance threshold where if the variance
            of length scales is less than interqubit separation, then algorithm
            terminates.

        Returns:
        --------
        posterior_f_state (`dtype` | numpy array):
            Posterior dephasing noise field estimate at each qubit on grid in
            vectorised form.
            Dims: len(self.GLOBALDICT["GRIDDICT"])
        posterior_r_state (`dtype` | numpy array):
            Posterior correlation length estimate at each qubit on grid in
            vectorised form.
            Dims: len(self.GLOBALDICT["GRIDDICT"])
        control_sequence (`dtype` | list ):
            List of control actions (qubits to measure) for every iteration step of qslam.
            Initialised as an empty list.

        '''

        qslamobj = qs.ParticleFilter(copy.deepcopy(self.GLOBALDICT),
                                     real_data=True,
                                     real_data_key=self.data_key)

        qslamobj.qslamr(max_num_iterations=self.GLOBALDICT["MODELDESIGN"]["MAX_NUM_ITERATIONS"],
                        measurements_controls=measurements_controls_,
                        autocontrol=autocontrol_,
                        var_thres=var_thres_)

        posterior_f_state = qslamobj.QubitGrid.get_all_nodes(["f_state"])
        posterior_r_state = qslamobj.QubitGrid.get_all_nodes(["r_state"])
        control_sequence = qslamobj.QubitGrid.control_sequence

        return posterior_f_state, posterior_r_state, control_sequence

    def naive_trial(self):
        '''
        Return recontructed maps using experimental measurement data.

            posterior_map (`dtype` | numpy array ):
                Reconstructed map using brute force measurement strategy but with
                a fixed measurement budget, where a subset of experimental
                measurements are sampled from the database.
                Dims: len(self.GLOBALDICT["GRIDDICT"]).

            true_map_ (`dtype` | numpy array ):
                Empirical map estimate using all experimental
                measurements in the database (infinite measurements).
                Dims: len(self.GLOBALDICT["GRIDDICT"]).
        '''

        RealDataObject = RealData(self.data_key)

        naiveobj = NaiveEstimatorExpt(RealDataObject,
                                      msmt_per_node=self.GLOBALDICT["MODELDESIGN"]["MSMTS_PER_NODE"],
                                      numofnodes=len(self.GLOBALDICT["GRIDDICT"]),
                                      max_num_iterations=self.GLOBALDICT["MODELDESIGN"]["MAX_NUM_ITERATIONS"])

        posterior_map, true_map_ = naiveobj.get_empirical_est()
        return posterior_map, true_map_

    def calculate_risk(self, number_of_trials=50):
        '''
        Return error matrices for map reconstructions for QSLAM and Naive
        measurement strategies using experimental measurement data.

        Parameters:
        -----------
            number_of_trials (`dtype`: int scalar):
                Number of repetitions of single run experiments.

        Returns:
        --------
            ssim_array (`dtype` | numpy array ):
                Mean SSIM score (0 - ideal, 1 - max error) for map reconstruction.
                Dims: 2. Col index: 0 - QSLAM, 1 - Naive.

            empr_array (`dtype` | numpy array ):
                Mean single map RMSE score (0 - ideal) for map reconstruction.
                Dims: 2. Col index: 0 - QSLAM, 1 - Naive.
        '''

        ssim_array = np.zeros((number_of_trials, 2))
        empr_array = np.zeros((number_of_trials, 2))

        for idx_run in range(number_of_trials):

            posterior_qslam_map = self.qslam_trial()[0]
            posterior_naive_map, true_map_ = self.naive_trial()
            posterior_map_list = [posterior_qslam_map, posterior_naive_map]

            for idx in range(2):

                residuals = posterior_map_list[idx] - true_map_
                ssim_array[idx_run, idx] = Metric.score_ssim(posterior_map_list[idx],
                                                             true_map_,
                                                             Cone=0.01, Ctwo=0.01)

                empr_array[idx_run, idx] = Metric.singlemap_rmse(residuals, axis=0)

        return np.mean(empr_array, axis=0), np.mean(ssim_array, axis=0)


class NaiveEstimatorExpt(object):
    '''
    Computes expected value of estimated map error using naive
    brute force measurement strategy and experimental measurement data.

    Attributes:
    -----------
    msmt_per_node (`dtype` | scalar int-like):
        Number of measurements per qubit per iteration before information is
            exchanged with neighbours.
    numofnodes (`dtype` | scalar int-like):
        Total number of qubit locations.
    max_num_iterations (`dtype` | scalar int-like):
        Max number of iterations for a qslam algorithm.
        A single control directive corressponds to one iteration.
    expt_data_generator (`dtype` | RealData class instance):
        A RealData class object instance.
    empirical_estimate (`dtype` | numpy array):
        Reconstructed map using brute force measurement strategy but with
        a fixed measurement budget, where a subset of experimental
        measurements are sampled from the database.
        Dims: numofnodes.


    Class Methods:
    --------------
    total_msmt_budget : Return total msmt budget given algorithm configurations.
    get_empirical_truth : Return true map implied by the empirical mean of msmt data.
    get_empirical_est : Return reconstructed map from Naive measurement strategy.

    '''

    def __init__(self,
                 RealDataObject,
                 msmt_per_node=1,
                 numofnodes=25,
                 max_num_iterations=None):
        '''
        Creates NaiveEstimatorExpt class object.

        Parameters:
        -----------
        RealDataObject (`dtype` | RealData class instance):
            A RealData class object instance.

        msmt_per_node (`dtype` | scalar int-like):
            Number of measurements per qubit per iteration before information is
                exchanged with neighbours. Default value: 1.s
        numofnodes (`dtype` | scalar int-like):
            Total number of qubit locations. Default value: 25.
        max_num_iterations (`dtype` | scalar int-like):
            Max number of iterations for a qslam algorithm. Default value: None.
        '''
        self.msmt_per_node = msmt_per_node
        self.numofnodes = numofnodes
        self.max_num_iterations = max_num_iterations
        self.expt_data_generator = RealDataObject
        self.empirical_estimate = None

        self.__total_msmt_budget = self.msmt_per_node * self.max_num_iterations


    @property
    def total_msmt_budget(self):
        ''' Return total msmt budget given algorithm configurations. '''
        return self.msmt_per_node * self.max_num_iterations

    def get_empirical_truth(self):
        ''' Return true map implied by the empirical mean of msmt data. '''
        empirical_mean = self.expt_data_generator.get_empirical_mean()
        true_map_ = Node.inverse_born(empirical_mean)
        return true_map_

    def get_empirical_est(self):
        ''' Return reconstructed map from Naive measurement strategy.

            self.empirical_estimate (`dtype` | numpy array ):
                Reconstructed map using brute force measurement strategy but with
                a fixed measurement budget, where a subset of experimental
                measurements are sampled from the database.
                Dims: numofnodes.

            true_map_ (`dtype` | numpy array ):
                Empirical map estimate using all experimental
                measurements in the database (infinite measurements).
                Dims: numofnodes.
        '''

        true_map_ = self.get_empirical_truth()

        # Measure all qubits if msmt budget is more than total number of qubits
        if self.numofnodes <= self.max_num_iterations:

            if self.max_num_iterations / self.numofnodes == self.msmt_per_node:
                mask = np.ones(self.numofnodes, dtype=bool)

            if self.max_num_iterations / self.numofnodes != self.msmt_per_node:
                self.msmt_per_node = int(self.total_msmt_budget / self.numofnodes)
                mask = np.ones(self.numofnodes, dtype=bool)

        # Measure a random subset of qubits if msmt budget is less than total number of qubits
        if self.numofnodes > self.max_num_iterations:

            randomly_choose = np.random.choice(self.numofnodes,
                                               self.max_num_iterations,
                                               replace=False)
            mask = np.zeros(self.numofnodes, dtype=bool) # Hiding all qubit locations.
            mask[randomly_choose] = True # Revealing only randomly chosen qubit locations.


        node_labels = np.arange(self.numofnodes)

        # MARKER JUNE 2019 - ALTERNATIVE INITIALISATION
        # self.empirical_estimate = np.random.randint(low=0.0,
        #                                             high=np.pi,
        #                                             size=self.numofnodes)

        # Original Initialisation
        self.empirical_estimate = np.ones(self.numofnodes) * np.random.randint(low=0.0,
                                                                               high=np.pi,
                                                                               size=1)

        for idx_node in node_labels[mask]:

            noisy_single_shots = [self.expt_data_generator.get_real_data(idx_node) for idx_shot in range(self.msmt_per_node)]
            self.empirical_estimate[idx_node] = Node.inverse_born(np.mean(np.asarray(noisy_single_shots, dtype=float)))

        return self.empirical_estimate, true_map_
