import qslamr as qs
import numpy as np

class SingleRunAnalysis2(object):
    ''' Class for analysing empirical distributions and step-wise reuslts in
        one single QSLAM run.'''
    def __init__(self,
                 SAMPLE_GLOBAL_MODEL,
                 true_map_,
                 repts=1,
                 measurements_controls_=None,
                 autocontrol_="ON",
                 var_thres_=1.0,
                 numofnodes=25,
                 beta_expansion_mode=True,
                 beta_skew_adjust=True):
        '''
        Parameters:
        -----------
            SAMPLE_GLOBAL_MODEL ('type': dictionary object):
                Dictionary object containing all QSLAM model design parameters,
                as defined in qslamdesignparams.py.

            true_map_ (`dtype` | numpy array):
                A set of phase values associated to each qubit location in an arrangement,
                where each phase take a value between [0, np.pi].

            repts  (`int`) :
                Number of repetitions for a Bayes Risk calculation.
                Default value: 1.

            measurements_controls : A list containing a measurement set and the control
                directive of the location of the measured qubit. Each control directive
                is a single iteration of the algorithm.
                Default value: None.

            autocontrol_( "OFF" / "ON" flag) :
                "OFF" - next qubit measured is specified as a user input via measurement_controls_
                "ON" - next qubit measured is chosen by the algorithm.
                Default value: "ON".

            var_thres :[NOT USED] Error variance threshold where if the variance
                of length scales is less than interqubit separation, then algorithm
                terminates.

            numofnodes (`type` | int scalar) :
                Total number of qubit locations in an arrangement.
                Default value: 25.

            beta_expansion_mode ('type': Boolean):
                If TRUE, QSLAM expands Beta particle layer using a truncated
                Gaussian distribution  at t based on posterior r_state information.
                If FALSE, QSLAM expands Beta particle layer by sampling from
                r_state prior distribution.
                Default value: TRUE

            beta_skew_adjust ('type': Boolean):
                If TRUE, QSLAM reports the mean or the mode of the r_state
                over a collection of beta particles for each Alpha parent.
                This setting can assist in mitigating influence of skewed
                beta distributions or those with numerical outliers.
                If FALSE, QSLAM reports the mean of the r_state for each Alpha
                parent.
                Default value: TRUE
        '''

        self.SAMPLE_GLOBAL_MODEL = SAMPLE_GLOBAL_MODEL
        self.true_map_ = true_map_
        self.repts = repts
        self.measurements_controls_ = measurements_controls_
        self.autocontrol_ = autocontrol_
        self.var_thres_ = var_thres_
        self.numofnodes = len(SAMPLE_GLOBAL_MODEL["GRIDDICT"])
        self.beta_expansion_mode = beta_expansion_mode
        self.beta_skew_adjust = beta_skew_adjust

    def call_qslam(self):
        ''' Return a full QSLAM object after a single run.'''

        qslamobj = 0
        qslamobj = qs.ParticleFilter(self.SAMPLE_GLOBAL_MODEL,
                                     save_run=True,
                                     beta_expansion_mode=self.beta_expansion_mode,
                                     beta_skew_adjust=self.beta_skew_adjust)

        qslamobj.QubitGrid.engineeredtruemap = self.true_map_

        qslamobj.qslamr(measurements_controls=self.measurements_controls_,
                        autocontrol=self.autocontrol_,
                        max_num_iterations=self.SAMPLE_GLOBAL_MODEL["MODELDESIGN"]["MAX_NUM_ITERATIONS"],
                        var_thres=self.var_thres_)

        return qslamobj

    def run_analysis(self, pathtofile ='./SingleRunAnalysis'):
        ''' Save a .npz file at pathtofile (with filename)containing the following
        data keys in a single run of the QSLAM algorithm:

        alpha_labels (`type` | list of int):
            Alpha particles labelled from 1, 2, ... M, where
            M = total number of alpha particles.
        beta_labels (`type` | list of int):
            Beta particles labelled from 1, 2, ... M, where
            M = total number of beta particles.
        joint_labels (`type` | list of int ):
            Joint Alpha-Beta particles labelled from 1, 2, ... M, where
            M = total number of alpha * total number of beta particles.

        The following data matrices have a common structure where array axes
            represent:
                Axis 0 : Number of iterations, t, of the algorithm
                Axis 1 : Number of qubit locations

        absolute_errors_matrix (`type` | numpy array):
            Absolute value of the difference between estimated f_state
            and a true map, with dimensions:

        rate_of_change_of_map (`type` | numpy array):
            Absolute value of the difference of the estimate f_state
            at t and t-1.

        rate_of_change_of_len (`type` | numpy array):
            Absolute value of the difference of the estimate r_state
            at t and t-1.


        predictive_weights (`type` | numpy array ):
            Weights of an empirical distribution representing the state of particles
                after a dynmical model update.

        leaf_weights (`type` | numpy array ):
            Weights of an empirical distribution representing the state of
                beta particles after expansion i.e. a beta likelihood score.

        joint_weights (`type` | numpy array ):
            Weights of an empirical distribution representing the state of particles
                as a joint distribution over both alpha and beta weights i.e.
                 a global likelhood score.

        marginalised_weights (`type` | numpy array ):
            Weights of an empirical distribution representing the state of alpha
                survivors after the beta particle distribution is marginalised.

        marginalised_states (`type` | numpy array ):
            Empirical distribution representing the state of alpha survivors
                after the beta particle distribution is marginalised.

        posterior_weights (`type` | numpy array ):
             Weights of an empirical posterior distribution representing
                the state of alpha particles after a second re-sampling step.

        posterior_states (`type` | numpy array ):
            Empirical posterior distribution representing the state of alpha particles
                after a second re-sampling step.

        '''

        max_num_iterations=self.SAMPLE_GLOBAL_MODEL["MODELDESIGN"]["MAX_NUM_ITERATIONS"]

        # --------------------------------------------------------------------
        # Make true error matrix and rate of change matrix for a single run
        # --------------------------------------------------------------------

        absolute_errors_matrix = np.zeros((max_num_iterations, self.numofnodes))
        rate_of_change_of_map = np.zeros((max_num_iterations, self.numofnodes))
        rate_of_change_of_len = np.zeros((max_num_iterations, self.numofnodes))

        qslamobj = self.call_qslam()

        for idxt in range(0, max_num_iterations):

            f_strt = self.numofnodes * 2  # fixed by QSLAM design of state vector
            f_end = self.numofnodes * 3
            r_strt = self.numofnodes * 3

            absolute_errors_matrix[idxt, :] = abs(qslamobj.save_posterior_state[idxt][f_strt:f_end] - self.true_map_)

            if idxt > 0:

                rate_of_change_of_map[idxt, :] = abs(qslamobj.save_posterior_state[idxt][f_strt:f_end] - qslamobj.save_posterior_state[idxt-1][f_strt:f_end])
                rate_of_change_of_len[idxt, :] = abs(qslamobj.save_posterior_state[idxt][r_strt:] - qslamobj.save_posterior_state[idxt-1][r_strt:])

        # --------------------------------------------------------------------
        # Make matrices for distribution of weights for the last repetition
        # --------------------------------------------------------------------


        P_ALPHA = self.SAMPLE_GLOBAL_MODEL["MODELDESIGN"]["P_ALPHA"]
        P_BETA = self.SAMPLE_GLOBAL_MODEL["MODELDESIGN"]["P_BETA"]

        predictive_weights = np.zeros((max_num_iterations, P_ALPHA))
        posterior_weights = np.zeros((max_num_iterations, P_ALPHA))
        leaf_weights = np.zeros((max_num_iterations, P_ALPHA, P_BETA))
        joint_weights = np.zeros((max_num_iterations, P_ALPHA * P_BETA))

        marginalised_weights = np.zeros((max_num_iterations, P_ALPHA))
        marginalised_states = np.zeros((max_num_iterations, P_ALPHA, self.numofnodes*4))
        posterior_states = np.zeros((max_num_iterations, P_ALPHA, self.numofnodes*4))

        alpha_labels = range(1, P_ALPHA + 1, 1)
        beta_labels = range(1, P_BETA + 1, 1)
        joint_labels = range(1, P_ALPHA*P_BETA + 1, 1)

        for idxt in range(0, max_num_iterations):

            predictive_weights[idxt, :] = np.asarray([qslamobj.save_alpha_predictive[idxt][idx].weight for idx in range(P_ALPHA)])

            start = P_ALPHA * idxt
            end = start + P_ALPHA
            leaf_weights[idxt, :, :] = np.asarray(qslamobj.save_all_beta_weights[start :  end])

            joint_weights[idxt, :] = qslamobj.save_alpha_beta_joint_weights[idxt]

            marginalised_labels = qslamobj.save_alpha_bar_2_oldnodes[idxt]

            for idx_m in range(len(marginalised_labels)):
                idx_particle = marginalised_labels[idx_m]
                marginalised_weights[idxt, idx_particle] = np.asarray(qslamobj.save_alpha_bar_2_weights[idxt])[idx_m]
                marginalised_states[idxt, idx_particle, :] = qslamobj.save_alpha_bar_2[idxt][idx_m].particle

            posterior_weights[idxt, :] = np.asarray([qslamobj.save_alpha_posterior[idxt][idx].weight for idx in range(P_ALPHA)])

            for idx in range(P_ALPHA):
                posterior_states[idxt, idx, :] = np.asarray([qslamobj.save_alpha_posterior[idxt][idx].particle])

        np.savez(pathtofile,
                 alpha_labels=alpha_labels,
                 beta_labels=beta_labels,
                 joint_labels=joint_labels,
                 absolute_errors_matrix=absolute_errors_matrix,
                 rate_of_change_of_map=rate_of_change_of_map,
                 rate_of_change_of_len=rate_of_change_of_len,
                 predictive_weights=predictive_weights,
                 leaf_weights=leaf_weights,
                 joint_weights=joint_weights,
                 marginalised_weights=marginalised_weights,
                 marginalised_states=marginalised_states,
                 posterior_weights=posterior_weights,
                 posterior_states=posterior_states)
        return
