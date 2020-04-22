'''
Created on Thu Apr 20 19:20:43 2017
@author: riddhisw

.. module:: riskanalysis

    :synopsis: Computes expected value of Bayes Risk for qslam.

    Module Level Classes:
    ----------------------
        EngineeredTruth: Generates true maps for Bayes Risk analysis
        CreateQslamExpt: Optimises qslam filter parameters and uses
            Bayes Risk metric for predictive power analysis.
        BayesRisk: Stores Bayes Risk map for a scenario specified
            by (testcase, variation).
        CreateNaiveExpt: Implements Naive Brute Force Measurement strategy and uses
            Bayes Risk metric for predictive power analysis.

.. moduleauthor:: Riddhi Gupta <riddhi.sw@gmail.com>
'''
import copy
import qslamr as qs
import numpy as np
import os

from scipy.interpolate import LinearNDInterpolator, Rbf

# PADUA COMPATIBILITY
import sys
sys.path.append('../paduaq/')
from pdinter_MM import pd_interpolant
from pdpoints import dims_padua_set

H_PARAM = ['LAMBDA_1', 'LAMBDA_2', 'SIGMOID_VAR', 'QUANT_VAR']
TYPE = 'multiquadric'

from hardware import Node
class EngineeredTruth(object):
    ''' Generates true maps for Bayes Risk analysis.'''
    def __init__(self,
                 numberofnodes,
                 TRUTHKWARGS):

        ''' Creates EngineeredTruth instance.
        Parameters:
        -----------
        number_of_nodes (`int` | scalar ):
            Total number of qubits on hardware.
        TRUTHKWARGS (`type`| dictionary object):
            Generates different true map configurations for map reconstruction
            simulations.

            TRUTHKWARGS["truthtype"]: Sets shape of true map.

                "Uniform" : sets the true field as uniform for a hardcoded value
                    of 0.456*np.pi. This value is arbitrary selected such that it
                    is not an extreme phase (0, pi) nor 0.5pi, but tests the
                    phase resolution of a map reconstruction procedure.

                "OneStep": sets the true field into 2 discontinuous regions
                    with a randomly chosen phase value in each region.
                    The ratio between two regions is set to 0.4.

                "OneStepd" : sets the true field into 2 discontinuous regions
                    with phase values TRUTHKWARGS["OneStepdheight"]["low"] and
                    TRUTHKWARGS["OneStepdheight"]["high"]. The number of qubit
                    nodes in low region relative to high region is given by
                    TRUTHKWARGS["OneStepdfloorarea"].

                "OneStepq" : sets the true field as 4 discontinous regions; of which
                     two regions have phase values set by TRUTHKWARGS["OneStepdheight"]["low"],
                     and two regions have phase values set by  RUTHKWARGS["OneStepdheight"]["high"].
                        e.g. For a linear array, this results in low | high | low |high
                        e.g. For a square grid with odd number of points on a edge,
                        this results in a quadrant where opposite corners have same phase,
                        size of low phase quadrants are larger than high phase quadrants.

                "Gaussian" : sets the true field as a 2D Gaussian with:
                    mu_x = 2.0 # mean in x
                    mu_y = 2.0 # mean in y
                    scl = 0.8 # scale

                "UseFunction" : sets true field using coordinate locations for each qubit

                    TRUTHKWARGS["all_qubit_locations"] : List of tuples [(x,y)] for coordinates
                    of qubit locations (sensing and data qubits)

                    TRUTHKWARGS["true_function"] : true_function() used for interpolation

                    TRUTHKWARGS["true_function_type"] : flag to specify calculation
                        inside true_function()

            TRUTHKWARGS["OneStepdheight"]["low"] : phase value of "low field",
                a scalar between 0 and np.pi, when "OneStepd" or "OneStepq" is selected.

            TRUTHKWARGS["OneStepdheight"]["high"] : phase value of "high field",
                a scalar between 0 and np.pi, when "OneStepd" or "OneStepq" is selected.

            TRUTHKWARGS["OneStepdfloorarea"] : Ratio of the number of qubits
            in low vs. high field region when truth type "OneStepdheight" is selected.



        Returns:
        --------
            true_map (`dtype` | numpy array):
                A set of phase values associated to each qubit location in an arrangement,
                where each phase take a value between [0, np.pi].
        '''

        self.type = TRUTHKWARGS["truthtype"]
        self.dims = numberofnodes
        self.TRUTHKWARGS = TRUTHKWARGS

    def get_map(self):

        if self.type == 'Uniform':

            truemap = np.ones(self.dims)*0.456*np.pi

        if self.type == 'OneStep':
            truemap = np.ones(self.dims)*np.random.uniform(low=0., high=np.pi*0.5)
            barrierdims = int(float(self.dims)* 0.4)
            truemap[barrierdims:] = np.ones(self.dims-barrierdims)*np.random.uniform(low=np.pi*0.5, high=np.pi)

        if self.type == 'OneStepd':
            lowfloor = self.TRUTHKWARGS["OneStepdheight"]["low"]
            highbarrier = self.TRUTHKWARGS["OneStepdheight"]["high"]
            floorarea = self.TRUTHKWARGS["OneStepdfloorarea"]
            truemap = np.ones(self.dims) * lowfloor
            barrierdims = int(float(self.dims)* floorarea)
            truemap[barrierdims:] = np.ones(self.dims-barrierdims) * highbarrier

        if self.type == 'OneStepq':
            tuner = self.TRUTHKWARGS["OneStepdheight"]["high"]
            truemap = np.ones(self.dims)*self.TRUTHKWARGS["OneStepdheight"]["low"]
            quarter = int((np.sqrt(self.dims)/2.0))
            grid = truemap.reshape(int(np.sqrt(self.dims)), int(np.sqrt(self.dims)))
            grid[0:quarter, 0:quarter] = tuner*np.ones_like(grid[0:quarter, 0:quarter])
            grid[quarter:, quarter:] = tuner*np.ones_like(grid[quarter:, quarter:])
            if np.sqrt(self.dims)/2.0 % 1 != 0:
                grid[quarter, :] = tuner*np.ones_like((grid[quarter, :]))
                grid[:, quarter] = tuner*np.ones_like((grid[:, quarter]))

        if self.type == 'Gaussian':

            mu_x = 2.0
            mu_y = 2.0
            scl = 0.8

            truemap = []
            sqrdims = int(np.sqrt(self.dims))
            for xidx in range(sqrdims):
                for yidx in range(sqrdims):
                    phase = 2.5*np.pi*(1.0 / (np.sqrt(2.0*np.pi*scl)))*np.exp(-((float(xidx) - mu_x)**2 + (float(yidx) - mu_y)**2)/ 2*scl)

                    if phase > np.pi:
                        phase = np.pi

                    if phase < 0.0:
                        phase = 0.0

                    truemap.append(phase)

            truemap = np.asarray(truemap)

        if self.type == "UseFunction":

            true_map = []
            true_function_type = self.TRUTHKWARGS["true_function_type"]
            
            # PADUA PROJECT: RANDOM POLYNOMIALS FUNCTIONALITY
            if true_function_type == 'randpoly':
                n = self.TRUTHKWARGS["randpoly"]["n"]
                trial = np.random.randint(low=0,high=199) 
            
            for idx_point in range(self.dims):
                point = self.TRUTHKWARGS["all_qubit_locations"][idx_point]
                
                # PADUA PROJECT: RANDOM POLYNOMIALS FUNCTIONALITY
                if true_function_type != 'randpoly':
                    true_map.append(self.TRUTHKWARGS["true_function"](point[0], point[1], d=true_function_type))
                
                if true_function_type == 'randpoly':
                    true_map.append(self.TRUTHKWARGS["true_function"](point[0], point[1], d=true_function_type,
                                                                      n=n,trial=trial))
            
            truemap = np.asarray(true_map)

        return truemap

class NaiveEstimator(object):
    '''
    Computes expected value of estimated map error using naive
    brute force measurement strategy and simulated measurements generated from an
    EngineeredTruth instance.

    Attributes:
    -----------
    msmt_per_node (`dtype` | scalar int):
        Number of measurements per qubit per iteration before information is
            exchanged with neighbours.
    numofnodes (`dtype` | scalar int):
        Total number of qubit locations.
    max_num_iterations (`dtype` | scalar int):
        Max number of iterations for a qslam algorithm.
        A single control directive corressponds to one iteration.
    truth_generator (`dtype` | RealData class instance):
        A EngineeredTruth class object instance.
    empirical_estimate (`dtype` | numpy array):
        Reconstructed map using brute force measurement strategy.
        Dims: numofnodes.


    Class Methods:
    --------------
    total_msmt_budget : Return total msmt budget given algorithm configurations.
    get_empirical_est : Return reconstructed map from Naive measurement strategy.
    '''
    def __init__(self,
                 TRUTHKWARGS,
                 msmt_per_node=1,
                 numofnodes=25,
                 data_qubits_indicies=None, # PADUA COMPATIBILITY - DATA QUBITS CANNOT BE MEASURED
                 intepolationflag=None, # PADUA COMPATIBILITY - LINEAR or PADUA INTERPOLATION
                 max_num_iterations=None):
        '''
        Creates a NaiveEstimator instance.

        Parameters:
        -----------
            TRUTHKWARGS (`type` | dictionary object):
                Dictionary object used to instantiate a EngineeredTruth instance.

            msmt_per_node (`dtype` | scalar int):
                Number of measurements per qubit per iteration before information is
                    exchanged with neighbours.
            numofnodes (`dtype` | scalar int):
                Total number of qubit locations.
            max_num_iterations (`dtype` | scalar int):
                Max number of iterations for a qslam algorithm.
                A single control directive corressponds to one iteration.
        '''

        self.msmt_per_node = msmt_per_node
        self.numofnodes = numofnodes
        self.max_num_iterations = max_num_iterations
        self.truth_generator = EngineeredTruth(self.numofnodes, TRUTHKWARGS)

        # PADUA
        self.data_qubits_indicies = data_qubits_indicies
        self.intepolationflag = intepolationflag
        self.all_qubit_locations = None
        if self.intepolationflag is not None:
            self.all_qubit_locations = self.truth_generator.TRUTHKWARGS["all_qubit_locations"]

        self.empirical_estimate = None

        self.__total_msmt_budget = self.msmt_per_node * self.max_num_iterations

    @property
    def total_msmt_budget(self):
        return self.msmt_per_node * self.max_num_iterations

    def get_empirical_est(self, addnoise=None):

        # PADUA COMPATIBILITY - TRUE MAP SUPPLIED
        # if self.true_map_supplied is None:
        phase_map = self.truth_generator.get_map()
        # if self.true_map_supplied is not None:
        #     phase_map = self.true_map_supplied

        # PADUA COMPATIBILITY - exclude data qubits
        num_of_sensing_nodes = self.numofnodes
        sensing_qubits = np.arange(num_of_sensing_nodes)
        if self.data_qubits_indicies is not None:
            sensing_qubits = np.asarray(list(set(sensing_qubits) - set(self.data_qubits_indicies)))
            num_of_sensing_nodes = len(sensing_qubits)

        mask = np.zeros(self.numofnodes, dtype=bool) # Mask for hiding all values.
        if num_of_sensing_nodes <= self.max_num_iterations:

            if self.max_num_iterations / num_of_sensing_nodes == self.msmt_per_node:
                mask[sensing_qubits] = True
            if self.max_num_iterations / num_of_sensing_nodes != self.msmt_per_node:
                self.msmt_per_node = int(self.total_msmt_budget / num_of_sensing_nodes)
                mask[sensing_qubits] = True

        if num_of_sensing_nodes > self.max_num_iterations:
            randomly_choose = np.random.choice(sensing_qubits, self.max_num_iterations, replace=False)
            mask[randomly_choose] = True

        # Final set of qubits that are measured under a Naive Approach
        final_measured_qubits = np.arange(self.numofnodes)[mask]
        self.empirical_estimate = np.ones(self.numofnodes) * np.random.random_sample(size=1) * np.pi

        for idx_node in final_measured_qubits:

            single_shots = [Node.quantiser(Node.born_rule(phase_map[idx_node])) for idx_shot in range(self.msmt_per_node)]

            if addnoise is not None:
                noisy_single_shots = addnoise["func"](single_shots, **addnoise["args"])
                self.empirical_estimate[idx_node] = Node.inverse_born(np.mean(np.asarray(noisy_single_shots, dtype=float)))

            elif addnoise is None:
                print "There is no noise process dictionary (NaiveEstimator.get_empirical_est)."
                self.empirical_estimate[idx_node] = Node.inverse_born(np.mean(np.asarray(single_shots, dtype=float)))

        # Add interpolation for unmeasured data qubits.
        if self.data_qubits_indicies is not None and self.intepolationflag is not None:

            f_data =  self.empirical_estimate[sensing_qubits]
            data_points = np.asarray([self.all_qubit_locations[idx_point] for idx_point in sensing_qubits])
            test_points = np.asarray([self.all_qubit_locations[idx_point] for idx_point in self.data_qubits_indicies])
                    
            if isinstance(self.intepolationflag, int):  # selects Padua interpolation specified by postive integer Padua order
                
                if self.intepolationflag > 0:
                    order = self.intepolationflag
                    f_interpolated_vals = np.diag(pd_interpolant(order, f_data, [test_points[:,0], test_points[:,1]]))

            if isinstance(self.intepolationflag, str): # selects interpolation techniques (non-Padua) specified by strings
                
                if self.intepolationflag == 'linear':
                    lin_interpolator = LinearNDInterpolator(data_points, f_data, fill_value=np.nan, rescale=False)
                    f_interpolated_vals = lin_interpolator.__call__(test_points)
                    
                if self.intepolationflag == 'Rbf':
                    rbf_interp = Rbf(data_points[:, 0], data_points[:, 1], f_data, function=TYPE)
                    f_interpolated_vals = rbf_interp(test_points[:, 0], test_points[:, 1])

            self.empirical_estimate[self.data_qubits_indicies] = f_interpolated_vals

        return self.empirical_estimate, phase_map

class Bayes_Risk(object):
    ''' Stores Bayes Risk map for a scenario specified by (testcase, variation)

    Attributes:
    ----------
        bayes_params (`dtype`) : Parameters to intiate Bayes Risk class:
            max_it_BR (`int`) : Number of repetitions for a Bayes Risk calculation.
            num_randparams (`int`) : Number of random (sigma, R) sample pairs.
            space_size (`int`) : Exponent parameter to set orders of magnitude
                spanned by unknown noise variance parameters.
            loss_truncation (`int`) : Pre-determined threshold for number of
                lowest input values to return in modcule function,
                get_tuned_params(), in module common.
        doparallel (`Boolean`) : Enable parallelisation of Bayes Risk calculations [DEPRECIATED].
        lowest_pred_BR_pair (`float64`) : (sigma, R) pair with min Bayes Risk in state estimation.
        lowest_fore_BR_pair (`float64`) : (sigma, R) pair with min Bayes Risk in prediction.
        means_list (`float64`) : Helper calculation for Bayes Risk.
        skip_msmts (`int`) : Number of time-steps to skip between measurements.
            To receive measurement at every time-step, set skip_msmts=1.
        did_BR_Map (`Boolean`) : When True, indicates that a BR Map has been created.
        macro_truth (`float64`) : Matrix data container for set of true noise realisations,
            generated for the purposes of calculating the Bayes Risk metric for all
            (sigma, R) random samples.
        macro_prediction_errors (`float64`) : Matrix data container for set of state estimates.
        macro_forecastng_errors (`float64`) : Matrix data containter for set of forecasts.
        macro_hyperparams (`float64`) : Matrix data containter for random
            samples of (sigma, R).
    '''

    def __init__(self, numofnodes, TRUTHKWARGS, RISKPARAMS):
        '''Initiates a Bayes_Risk class instance.
        Parameters:
        -----------
        numofnodes (`dtype` | scalar int):
                Total number of qubit locations.
        TRUTHKWARGS (`type`| dictionary object):
            Dictionary object used to instantiate an EngineeredTruth instance.
        RISKPARAMS (`type`| dictionary object):
            Dictionary object as given in qslamdesignparams.py that sets
            parameters to conduct an empirical risk analysis.

        '''

        self.savetopath = RISKPARAMS["savetopath"]
        self.max_it_BR = RISKPARAMS["max_it_BR"]
        self.max_it_qslam = RISKPARAMS["max_it_qslam"]
        self.num_randparams = RISKPARAMS["num_randparams"]
        self.space_size = RISKPARAMS["space_size"]
        self.loss_truncation = RISKPARAMS["loss_truncation"]
        self.numofnodes = numofnodes # Default

        self.truemap_generator = EngineeredTruth(self.numofnodes, TRUTHKWARGS) # PADUA 

        self.filename_br = None
        self.macro_true_fstate = None
        self.macro_predictions = None
        self.macro_residuals = None
        self.macro_hyperparams = None
        self.lowest_pred_BR_pair = None
        self.did_BR_Map = True
        self.means_list = None

        pass

class CreateQslamExpt(Bayes_Risk):
    '''Scans QSLAM filter parameters and uses
       Bayes Risk metric for predictive power analysis.

    Class Attributes:
    -----------------
        qslamobj (`type`| ParticleFilter object):
            Stores a QSLAM ParticleFilter instance.

        filename_br (`type`| string ):
            Sets filename for storing the results generated by Bayes
            Risk analysis of QSLAM.

        GLOBALDICT  (`type`| dictionary object):
            Dictionary object storing QSLAM model parameters,
            as specified in qslamdesignparams.py

    Class Methods:
    --------------
        loss : Return residual errors between
            engineered truth and algorithm posterior state.

        map_loss_trial : Return posterior f_state, posterior r_state, residuals, and
            control path for map reconstruction from one trial of QSLAM.

        rand_param : Return a randomly sampled hyper-parameter vector.

        modify_global_model : Return new model parameters for QSLAM.
            Helper function for modifying QSLAM parameters during brute
            force optimisation.

        one_bayes_trial : Return true realisations, state etimation errors and preditions
            ver max_it_BR repetitions for one choice of QSLAM model configuration.

        naive_implementation : Return Bayes Risk analysis as a saved .npz file over max_it_BR
            repetitions of true dephasing noise and simulated datasets; for
            num_randparams number of random hyperparameters.
    '''

    def __init__(self, TRUTHKWARGS, GLOBALDICT):
        '''
        Initiates a CreateQslamExpt class instance.

        Parameters:
        -----------
        TRUTHKWARGS (`type`| dictionary object):
            Dictionary object used to instantiate an EngineeredTruth instance.
        GLOBALDICT (`type`| dictionary object):
            Dictionary object as given in qslamdesignparams.py that sets
            parameters to conduct an empirical risk analysis and configure QSLAM.

        '''
        self.GLOBALDICT = GLOBALDICT
        RISKPARAMS = self.GLOBALDICT["RISKPARAMS"]
        Bayes_Risk.__init__(self,
                            len(GLOBALDICT["GRIDDICT"]),
                            TRUTHKWARGS,
                            RISKPARAMS)
        self.qslamobj = None
        self.filename_br = self.GLOBALDICT["MODELDESIGN"]["ID"] + '_BR_Map'


    def loss(self, posterior_state, true_state_):
        '''Return residual errors between
        engineered truth and algorithm posterior state.

        Parameters:
        ----------
            posterior_state (`float64` | Numpy array | dims: N) :
                Posterior state estimates in vector form.
            true_state_ (`float64` | Numpy array | dims: N) :
                One realisation of (engineered) true state in vector form.
        Returns:
        -------
            residuals_errors : difference between algorithm
                output and engineered truth map.
        '''
        residuals_errors = posterior_state - true_state_
        return residuals_errors

    def map_loss_trial(self,
                       true_map_,
                       SAMPLE_GLOBAL_MODEL,
                       measurements_controls_=None,
                       autocontrol_="ON",
                       var_thres_=1.0):

        '''Return posterior f_state, posterior r_state, residuals, and
        control path for map reconstruction from one trial of QSLAM.

        Parameters:
        -----------
            true_map_ (`dtype`| numpy array ):
                A set of phase values associated to each qubit location in an arrangement,
                where each phase take a value between [0, np.pi].

            SAMPLE_GLOBAL_MODEL (`dtype`| numpy array ):
                Dictionary object as given in qslamdesignparams.py that sets
                parameters to conduct an empirical risk analysis and configure QSLAM.

            measurements_controls_ (`dtype`| numpy array ):
                A list containing a measurement set and the control
                directive of the location of the measured qubit. Each control directive
                is a single iteration of the algorithm.

            autocontrol_ (`dtype`| string ):
                "OFF" - next qubit measured is specified as a user input via measurement_controls_
                "ON" - next qubit measured is chosen by the algorithm.
                Default value: "OFF".

            var_thres_ (`dtype`| numpy array ):
                Error variance threshold where if the variance
                of length scales is less than interqubit separation, then algorithm
                terminates.

        Returns:
        -------

            posterior_map (`type` | numpy array ):
                Posterior f_state at the end of a QSLAM filtering run.
            map_residuals (`type` | numpy array ):
                Difference between posterior f_state  and a true map at the end
                    of a QSLAM filtering run.
            posterior_corrs (`type` | numpy array ):
                Posterior r_state at the end of a QSLAM filtering run.
            controlpath (`type` | list ):
                List of qubit locations at which a physical measurement was performed.
                A single element of this list corresponds to one iteration of QSLAM;
                length of list output at the end of a QSLAM filtering run matches
                total number of algorithm iterations.

        '''

        self.qslamobj = qs.ParticleFilter(SAMPLE_GLOBAL_MODEL)
        self.qslamobj.QubitGrid.engineeredtruemap = true_map_

        self.qslamobj.qslamr(measurements_controls=measurements_controls_,
                             autocontrol=autocontrol_,
                             max_num_iterations=SAMPLE_GLOBAL_MODEL["MODELDESIGN"]["MAX_NUM_ITERATIONS"],
                             var_thres=var_thres_)

        posterior_map = self.qslamobj.QubitGrid.get_all_nodes(["f_state"])
        posterior_corrs = self.qslamobj.QubitGrid.get_all_nodes(["r_state"])

        map_residuals = self.loss(posterior_map, true_map_)
        controlpath = self.qslamobj.QubitGrid.control_sequence

        return posterior_map, map_residuals, posterior_corrs, controlpath

    def rand_param(self, SAMPLE_GLOBAL_MODEL):
        ''' Return a randomly sampled hyper-parameter vector.
        Helper function to one_bayes_trial().

        Parameters:
        ----------
            SAMPLE_GLOBAL_MODEL (`type`| dictionary object):
                Dictionary object storing QSLAM model parameters,
                as specified in qslamdesignparams.py.

        Returns:
        --------
            samples (`type`| numpy array):
                Choice of model parameters sampled from probability
                distributions for hyper-parameters, as specified in
                in SAMPLE_GLOBAL_MODEL["HYPERDICT"]:
                    samples[0]: Random variate for LAMBDA_1
                    samples[1]: Random variate for LAMBDA_2
                    samples[2]: Random variate for SIGMOID_APPROX_ERROR
                    samples[3]: Random variate for QUANTISATION_UNCERTY
        '''

        HYPERDICT = SAMPLE_GLOBAL_MODEL["HYPERDICT"]
        samples = [HYPERDICT["DIST"](space_size=self.space_size, **HYPERDICT["ARGS"][param]) for param in H_PARAM]
        return samples

    def modify_global_model(self, samples, SAMPLE_GLOBAL_MODEL):
        '''Return new model parameters for QSLAM.
        Helper function for modifying QSLAM parameters during brute
        force optimisation.

        Parameters:
        ----------
            SAMPLE_GLOBAL_MODEL (`type`| dictionary object):
                Dictionary object storing QSLAM model parameters,
                as specified in qslamdesignparams.py.
            samples (`type`| numpy array):
                Choice of model parameters in GLOBALDICT for this Bayes trial:
                    samples[0]: LAMBDA_1
                    samples[1]: LAMBDA_2
                    samples[2]: SIGMOID_APPROX_ERROR
                    samples[3]: QUANTISATION_UNCERTY
        Returns:
        --------
            SAMPLE_GLOBAL_MODEL : New model parameters for QSLAM with dictionary
                elememts LAMBDA_1, LAMBDA_2, SIGMOID_APPROX_ERROR and QUANTISATION_UNCERTY
                updated with the values in samples.
        '''

        SAMPLE_GLOBAL_MODEL["MODELDESIGN"]["LAMBDA_1"] = samples[0]
        SAMPLE_GLOBAL_MODEL["MODELDESIGN"]["LAMBDA_2"] = samples[1]
        SAMPLE_GLOBAL_MODEL["NOISEPARAMS"]["SIGMOID_APPROX_ERROR"]["SIGMA"] = samples[2]
        SAMPLE_GLOBAL_MODEL["NOISEPARAMS"]["QUANTISATION_UNCERTY"]["SIGMA"] = samples[3]

        return SAMPLE_GLOBAL_MODEL

    def one_bayes_trial(self, samples=None):
        ''' Return true realisations, state etimation errors and preditions
        over max_it_BR repetitions for one choice of QSLAM model configuration.

        Parameters:
        ----------
            samples (`type`| numpy array):
                Choice of model parameters in GLOBALDICT for this Bayes trial:
                    samples[0]: LAMBDA_1
                    samples[1]: LAMBDA_2
                    samples[2]: SIGMOID_APPROX_ERROR
                    samples[3]: QUANTISATION_UNCERTY
                If NONE, then these samples are randomly selected from:
                    rand_param() using hyper-parameter distributions set in
                    GLOBALDICT["HYPERDICT"].
                Default value: None.

        Return:
        -------
            true_maps (`type`| numpy array):
                Set of true_maps for each trial in the Bayes Risk calculation,
                for a given configuration of the QSLAM model set by GLOBALDICT
                and samples.
            predictions (`type`| numpy array):
                Set of posterior f_state for each trial in the Bayes Risk calculation,
                for a given configuration of the QSLAM model set by GLOBALDICT
                and samples.
            map_errors (`type`| numpy array):
                Set of residual map errors for each trial in the Bayes Risk calculation,
                for a given configuration of the QSLAM model set by GLOBALDICT
                and samples.
            samples (`type`| numpy array):
                Choice of model parameters for LAMBDA_1, LAMBDA_2, SIGMOID_APPROX_ERRORs
                and QUANTISATION_UNCERTY in GLOBALDICT for this Bayes trial.
            correlations (`type`| numpy array):
                Set of posterior r_state for each trial in the Bayes Risk calculation,
                for a given configuration of the QSLAM model set by GLOBALDICT
                and samples.
            paths (`type`| numpy array):
                Set of control paths for each trial in the Bayes Risk calculation,
                for a given configuration of the QSLAM model set by GLOBALDICT
                and samples.
        '''

        SAMPLE_GLOBAL_MODEL = copy.deepcopy(self.GLOBALDICT)

        if samples is None:
            samples = self.rand_param(SAMPLE_GLOBAL_MODEL)

        SAMPLE_GLOBAL_MODEL = self.modify_global_model(samples, SAMPLE_GLOBAL_MODEL)

        predictions = []
        map_errors = []
        true_maps = []
        correlations = []
        paths = []

        for ind in xrange(self.max_it_BR):

            true_map_ = self.truemap_generator.get_map()
            posterior, errors, posterior_corrs, controlpath = self.map_loss_trial(true_map_, SAMPLE_GLOBAL_MODEL)

            true_maps.append(true_map_)
            predictions.append(posterior)
            map_errors.append(errors)
            correlations.append(posterior_corrs)
            paths.append(controlpath)

        return true_maps, predictions, map_errors, samples, correlations, paths

    def naive_implementation(self, randomise='OFF'):
        ''' Return Bayes Risk analysis as a saved .npz file over max_it_BR
        repetitions of true dephasing noise and simulated datasets; for
        num_randparams number of random hyperparameters.

        Parameters:
        -----------
            randomise (`type` | string):
                "OFF" : LAMBDA_1, LAMBDA_2, SIGMOID_APPROX_ERRORs
                and QUANTISATION_UNCERTY are set by GLOBAL DICT.
                Otherwise, these parameters are reset in GLOBALDICT by sampling
                random variates using randparams().

        Returns:
        -------
            Output .npz file containing all Bayes Risk data for analysis.
        '''
        fix_hyperparams = None
        self.macro_hyperparams = []
        self.macro_predictions = []
        self.macro_residuals = []
        self.macro_true_fstate = []
        self.macro_correlations = []
        self.macro_paths = []

        # start_outer_multp = t.time()

        for ind in xrange(self.num_randparams):

            if randomise == 'OFF':
                fix_hyperparams = np.ones(4)
                fix_hyperparams[0] = self.GLOBALDICT["MODELDESIGN"]["LAMBDA_1"]
                fix_hyperparams[1] = self.GLOBALDICT["MODELDESIGN"]["LAMBDA_2"]
                fix_hyperparams[2] = self.GLOBALDICT["NOISEPARAMS"]["SIGMOID_APPROX_ERROR"]["SIGMA"]
                fix_hyperparams[3] = self.GLOBALDICT["NOISEPARAMS"]["QUANTISATION_UNCERTY"]["SIGMA"]

            full_bayes_map = self.one_bayes_trial(samples=fix_hyperparams)

            self.macro_true_fstate.append(full_bayes_map[0])
            self.macro_predictions.append(full_bayes_map[1])
            self.macro_residuals.append(full_bayes_map[2])
            self.macro_hyperparams.append(full_bayes_map[3])
            self.macro_correlations.append(full_bayes_map[4])
            self.macro_paths.append(full_bayes_map[5])

            np.savez(os.path.join(self.savetopath, self.filename_br),
                     macro_true_fstate=self.macro_true_fstate,
                     macro_predictions=self.macro_predictions,
                     macro_residuals=self.macro_residuals,
                     macro_correlations=self.macro_correlations,
                     macro_paths=self.macro_paths,
                     macro_hyperparams=self.macro_hyperparams,
                     max_it_BR=self.max_it_BR,
                     num_randparams=self.num_randparams,
                     savetopath=self.savetopath)

        self.did_BR_Map = True


class CreateNaiveExpt(Bayes_Risk):
    '''Implements Naive Brute Force Measurement strategy and uses
      Bayes Risk metric for predictive power analysis.

    Class Attributes:
    -----------------
        naiveobj (`type`| ParticleFilter object):
            Stores a NaiveEstimator instance.

        filename_br (`type`| string ):
            Sets filename for storing the results generated by Bayes
            Risk analysis of Naive brute force measurement.

        GLOBALDICT  (`type`| dictionary object):
            Dictionary object as specified in qslamdesignparams.py
            to conduct an empirical risk analysis.

      '''

    def __init__(self, TRUTHKWARGS, GLOBALDICT):
        ''' Creates CreateNaiveExpt instance.

        Parameters:
        -----------
        TRUTHKWARGS (`type`| dictionary object):
            Dictionary object used to instantiate an EngineeredTruth instance.
        GLOBALDICT (`type`| dictionary object):
            Dictionary object as specified in qslamdesignparams.py
            to conduct an empirical risk analysis.
        '''

        self.GLOBALDICT = GLOBALDICT
        RISKPARAMS = self.GLOBALDICT["RISKPARAMS"]
        Bayes_Risk.__init__(self,
                            len(self.GLOBALDICT["GRIDDICT"]),
                            TRUTHKWARGS,
                            RISKPARAMS)
        self.naiveobj = None
        self.filename_br = self.GLOBALDICT["MODELDESIGN"]["ID"] + '_NE_Map'

    def loss(self, posterior_state, true_state_):
        ''' Return residual errors between
            engineered truth and algorithm posterior state.

        Parameters:
        ----------
            posterior_state (`float64` | Numpy array | dims: N) :
                Empirical state estimated from brute force measurement
                in vector form.
            true_state_ (`float64` | Numpy array | dims: N) :
                One realisation of (engineered) true state in vector form.
        Returns:
        -------
            residuals_errors : difference between algorithm
                output and engineered true map.
        '''
        residuals_errors = posterior_state - true_state_
        return residuals_errors


    def map_loss_trial(self):

        '''Return true (enginneered) map, estimate of f_state, and residual errors
         for map reconstruction from one trial of Naive brute force measurement
         algorithm.

         Returns:
         --------
            true_map_ (`dtype`| numpy array ):
                A set of phase values associated to each qubit location in an arrangement,
                where each phase take a value between [0, np.pi].
            posterior_map (`type` | numpy array ):
                Estimate f_state at the end of a Naive brute force measurement run.
            map_residuals (`type` | numpy array ):
                Difference between estimated f_state  and a true map at the end
                    of a Naive brute force measurement run.
         '''
        KWARGS = {}
        KWARGS["msmt_per_node"] = self.GLOBALDICT["MODELDESIGN"]["MSMTS_PER_NODE"]
        KWARGS["numofnodes"] = len(self.GLOBALDICT["GRIDDICT"])
        KWARGS["max_num_iterations"] = self.GLOBALDICT["MODELDESIGN"]["MAX_NUM_ITERATIONS"]

        if self.GLOBALDICT["DATA_QUBITS"] is not None:
            KWARGS["data_qubits_indicies"] = self.GLOBALDICT["DATA_QUBITS"]

        if self.GLOBALDICT["INTERPOLATE_FLAG"] is not None:
            KWARGS["intepolationflag"] = self.GLOBALDICT["INTERPOLATE_FLAG"]

        self.naiveobj = NaiveEstimator(self.truemap_generator.TRUTHKWARGS, **KWARGS)

        posterior_map, true_map_ = self.naiveobj.get_empirical_est(addnoise=self.GLOBALDICT["ADDNOISE"])

        map_residuals = self.loss(posterior_map, true_map_)

        return true_map_, posterior_map, map_residuals


    def one_bayes_trial(self):
        ''' Return true realisations, state etimation errors and preditions
        over max_it_BR repetitions for one choice of Naive Brute Force measurement
        model configuration.

        Returns:
        --------
            true_maps (`type`| numpy array):
                Set of true_maps for each trial in the Bayes Risk calculation,
                for a given configuration of the Naive brute force measurement
                experimnent.
            predictions (`type`| numpy array):
                Set of estimated f_state for each trial in the Bayes Risk calculation,
                for a given configuration of the Naive brute force measurement
                experimnent.
            map_errors (`type`| numpy array):
                Set of residual map errors for each trial in the Bayes Risk calculation,
                for a given configuration of the Naive brute force measurement
                experimnent.


        '''

        predictions = []
        map_errors = []
        true_maps = []

        for ind in xrange(self.max_it_BR):

            true_map_, posterior, errors = self.map_loss_trial()

            true_maps.append(true_map_)
            predictions.append(posterior)
            map_errors.append(errors)

        return true_maps, predictions, map_errors

    def naive_implementation(self, num_randparams=1):
        ''' Return Bayes Risk analysis as a saved .npz file over max_it_BR
        repetitions of true dephasing noise and simulated datasets.

        Parameters:
        -----------
            num_randparams (`type`: int):
                Number of configurations of the algorithm for which one
                Bayes trial needs to conducted.  Naive brute force measurement
                has only one (default) configuration and no other free parameters.
                Default value: 1.

        Returns:
        -------
            Output .npz file containing all Bayes Risk data for analysis.
        '''
        self.num_randparams = num_randparams
        self.macro_predictions = []
        self.macro_residuals = []
        self.macro_true_fstate = []

        for ind in xrange(self.num_randparams):

            full_bayes_map = self.one_bayes_trial()

            self.macro_true_fstate.append(full_bayes_map[0])
            self.macro_predictions.append(full_bayes_map[1])
            self.macro_residuals.append(full_bayes_map[2])

            np.savez(os.path.join(self.savetopath, self.filename_br),
                     macro_true_fstate=self.macro_true_fstate,
                     macro_predictions=self.macro_predictions,
                     macro_residuals=self.macro_residuals,
                     max_it_BR=self.max_it_BR,
                     num_randparams=self.num_randparams,
                     savetopath=self.savetopath)
        self.did_BR_Map = True
