'''
Created on Thu Apr 20 19:20:43 2017
@author: riddhisw

.. module:: qslamr

    :synopsis: Core particle filter to implement SLAM framework using local
        measurements on a 2D arrangement of qubits to reconstruct a spatially
        varying dephasing noise field.

    Module Level Classes:
    ----------------------
        ParticleFilter : Conducts particle filtering under the qslamr framework.

.. moduleauthor:: Riddhi Gupta <riddhi.sw@gmail.com>
'''

from itertools import combinations
import numpy as np
from scipy.stats import mode, truncnorm

from hardware import Grid, PARTICLE_STATE
from control_action import controller
from hardware import Node
from particlesets import AlphaParticle, ParticleSet, BetaParticle
from particleweightcalcs import ParticleLikelihoods as pl

class ParticleFilter(Grid):
    ''' Conducts particle filtering under the qslamr framework. Inherits from
        class hardware.Grid.

    Attributes:
    ----------
        QubitGrid (`Grid` class object):
            `Grid` class object representing a 2D spatial arrangement of qubits,
            tracking  posterior state estimates, physical and quasi measurements.
        dgrid (`float` | scalar):
            Max spatial pairwise separation for any two qubits on the Grid.
        d_iterq (`float` | scalar):
            Min spatial pairwise separation for any two qubits on the Grid.
        R_MAX (`float` | scalar):
            Upper bound on maximal correlation lengths physically expected in the system.
        R_MIN (`float` | scalar):
            Lower bound on maximal correlation lengths physically expected in the system.
        MODELDESIGN (Dictionary object):
            Stores a set of free parameters and initial conditions.
        pset_alpha (`int` | scalar):
            Total number of of Alpha particles (conserved during re-sampling).
        pset_beta (`int` | scalar):
            Total number of of Beta particles for each Alpha parent.
        resample_thresh (`float` | scalar):
            Effective number of particles given variance of posterior weights.
            Can be used to introduce an adaptive re-sampling scheme.
            [NOT USED : Currently resampling at each iteration].
        measurements_controls (`float` | dims: 2):
             External input for the most recent control directive and measurement
             outcome for the qubit, denoted by locaton index, node_j.
        empty_alpha_particles = [AlphaParticle() for idx in range(self.MODELDESIGN["P_ALPHA"])]
        self.AlphaSet = ParticleSet(empty_alpha_particles, **WEIGHTFUNCDICT_ALPHA)

    Class Methods:
    -------------
        qslamr : Execute core numerical QSLAM solver for qslamr module.
        PerformMeasurement : Return measurement outcome and control directive
            i.e. the location of  the qubit node that is measured.
        RunTwoLayerParticleFilter : Return next iteration of particle filtering algorithm after processing
            an input measurement outcome and control directive.
        update_qubitgrid_via_quasimsmts : Return the posterior neighbourhood associated
            with a qubit measurement at node_j and update all neighbours of node_j
            with quasi-measurement information from the posterior state estimates at node_j.
        InitializeParticles : Intialise properties of QubitGrid and AlphaSet for t = 0 time step
            set_init_alphaparticle : Set initial properties of an Alpha particle.
            Helper function to InitializeParticles().
        ReceiveMsmt : Updates global qubit hardware variables and properties of Alpha
            particles for a new measurement outcome and control directive.
        PropagateState: Propagates all Alpha particle states according to (apriori) transition
            propability distributions for propagating states between measurements.
        ComputeAlphaWeights : Return weights for the distribution defined over Alpha particles.
        ComputeJointWeights : Return posterior weights for the joint distribution
            defined over both Alpha and Beta particles.

            generate_beta_layer : Return a set of Beta particle weights for a given Alpha parent.
                Helper function to ComputeJointWeights()
            sample_radii : Return a r_state sample to generate a new Beta particle for a
                given Alpha parent.
                Helper function to generate_beta_layer().
            generate_beta_neighbourhood : Generate phase estimates over a neighbourhood
                for a candidate Beta particle and its Alpha parent. Helper function
                to generate_beta_layer().

        ExecuteParticleBranching : Return a new set of Alpha particles with uniform weights resampled
            according to the joint probability of both Alpha and Beta particles, subject to
            conserving Alpha particle number.

            collapse_beta: Return resampled Alpha parents with collapsed Beta layers. For each Alpha
                parent, store the mode and the spread of r_states in Beta-Particle layer as
                the posterior r_state at node_j. Helper function for ExecuteParticleBranching().
            effective_particle_size : Return effective particle number based on posterior weight
                variances for adapative resampling. [Currently not used].

    Static Methods:
    --------------
        resample_constant_pset_alpha (static method) : Return indicies for Alpha
            particles resampled from posterior Alpha-Beta weights, while conserving
            Alpha particle number. Helper function for ExecuteParticleBranching().
        calc_posterior_lengthscale (static method) : Return descriptive moments
            for the distribution of r_states in a Beta particle set for an Alpha parent,
            post resampling. Helper function for collapse_beta().
        compute_dist : Return Euclidean distance between two points.
            Helper function for get_distances().
        get_distances :  Return all pairwise distances between an arrangement of qubits.
            Helper function for ParticleFilter.find_max_distance() and
            ParticleFilter.find_min_distance.
        find_max_distance : Return the pair of positions with maximum pair-wise
            separation distance.
        find_min_distance : Return  the pair of positions with minimum pair-wise s
            eparation distance.
        get_alpha_node_from_treeleaf  (static method) : Return Alpha Particle index
            based on global index for flattened layers of Alpha-Beta particles.
            Helper function for collapse_beta().
        get_beta_node_from_treeleaf (static method) : Return Beta Particle index based
            on global index for flattened layers of Alpha-Beta particles.
            Helper function for collapse_beta().
        resample_from_weights (static method) : Return samples from the posterior
            distribution of states approximated by a discrete set of posterior normalised
            weights. Helper function for ParticleFilter.resample_constant_pset_alpha().
        GetHardwareConfig (static method) : Return max qubit pairwise separation
            on hardware, and lower and upper bounds on uniform distribution of r_states as a prior.
        get_subtrees (static method) : Return a list of  pairs of (start, end) for a sub-tree
            where each sub-tree represents an Alpha particle, and its element
            leaves represent suriviving resampled Beta particles.
            Helper function for ExecuteParticleBranching.

    '''

    def __init__(self,
                 GLOBALDICT,
                 save_run=False,
                 beta_expansion_mode=False,
                 beta_skew_adjust=True,
                 real_data=False,
                 real_data_key=None):
        '''
        Creates a ParticleFilter instance.

        Parameters:
        ----------
            GLOBALDICT ('type': dictionary object):
                Dictionary object containing all QSLAM model design parameters,
                as defined in qslamdesignparams.py.
            save_run ('type': Boolean):
                If TRUE, QSLAM saves intermediary particle distribution information
                for each iteration, t, as additional class attributes.
                Default value: FALSE
            beta_expansion_mode ('type': Boolean):
                If TRUE, QSLAM expands Beta particle layer using a truncated
                Gaussian distribution  at t based on posterior r_state information.
                If FALSE, QSLAM expands Beta particle layer by sampling from
                r_state prior distribution.
                Default value: FALSE
            beta_skew_adjust ('type': Boolean):
                If TRUE, QSLAM reports the mean or the mode of the r_state
                over a collection of beta particles for each Alpha parent.
                This setting can assist in mitigating influence of skewed
                beta distributions or those with numerical outliers.
                If FALSE, QSLAM reports the mean of the r_state for each Alpha
                parent.
                Default value: TRUE
            real_data ('type': Boolean):
                If TRUE, QSLAM samples an experimental measurement database to
                simulate a spatial mapping experiment.
                Default value: FALSE
            real_data_key ('type': int):
                If real_data is TRUE, real_data_key selects experimental database
                for a simulated QSLAM experiment.
                Default value: None
        '''

        self.GLOBALDICT = GLOBALDICT
        self.MODELDESIGN = self.GLOBALDICT["MODELDESIGN"]
        NOISEPARAMS = self.GLOBALDICT["NOISEPARAMS"]

        self.save_run = save_run
        self.beta_expansion_mode = beta_expansion_mode
        self.skew_adjust = beta_skew_adjust

        if self.save_run is True:

            self.save_run_variance = [] # beta particle variance at location t-1
            self.save_controls_j = [] # pairs of (qubit msmt, qubit location)
            self.save_alpha_predictive = [] # alpha particles after dynanmic update == posterior state
            self.save_alpha_bar_1 = [] # alpha particles after h1 update
            self.save_all_beta_weights = [] # long list of all beta wieghts
            self.save_alpha_beta_joint = [] # joint particles after beta expansion
            self.save_alpha_beta_joint_weights = []  # joint particles after beta expansion
            self.save_alpha_bar_2 = []  # alpha particles after beta collapse
            self.save_alpha_bar_2_oldnodes = [] # links alpha particles to alpha_bar_1
            self.save_alpha_bar_2_weights = [] # alpha weights after beta collapse
            self.save_alpha_posterior = [] # final particle set
            self.save_posterior_state = [] # final posterior state

        self.LikelihoodObj = pl(**NOISEPARAMS)
        self.measurements_controls = None

        self.PRIORDICT = self.GLOBALDICT["PRIORDICT"]

        poskeys = sorted(self.GLOBALDICT["GRIDDICT"].keys())
        posvals = [self.GLOBALDICT["GRIDDICT"][idx_key] for idx_key in poskeys]
        LAMBDA_1 = self.MODELDESIGN["LAMBDA_1"]
        SAMPLE_F = self.PRIORDICT["SAMPLE_F"]

        ADDNOISE = self.GLOBALDICT["ADDNOISE"]
        self.QubitGrid = Grid(LAMBDA_1=LAMBDA_1,
                              addnoise=ADDNOISE,
                              real_data=real_data,
                              real_data_key=real_data_key,
                              list_of_nodes_positions=posvals,
                              **SAMPLE_F)

        self.dgrid, self.diterq, self.R_MIN, self.R_MAX = self.GetHardwareConfig()

        empty_alpha_particles = [AlphaParticle() for idx in range(self.MODELDESIGN["P_ALPHA"])]
        self.AlphaSet = ParticleSet(empty_alpha_particles,
                                    **self.LikelihoodObj.WEIGHTFUNCDICT_ALPHA)

# #############################################################################
# QSLAM
# #############################################################################

    def qslamr(self,
               measurements_controls=None,
               autocontrol="OFF",
               max_num_iterations=None,
               list_of_dataqubits=None,
               var_thres=1.0):
        ''' Execute core QSLAM module.

        Parameters:
        ----------
            measurements_controls : A list containing a measurement set and the control
                directive of the location of the measured qubit. Each control directive
                is a single iteration of the algorithm.
            autocontrol_( "OFF" / "ON" flag) :
                "OFF" - next qubit measured is specified as a user input via measurement_controls_
                "ON" - next qubit measured is chosen by the algorithm.
                Default value: "OFF".

            max_num_iterations : Maximum number of iterations at which the algorithm terminates.

            list_of_dataqubits (`int`| list) :
                List of indices for qubit locations for data qubits that
                cannot be used for sensing measurements.
                Defaults to None.

            var_thres :[NOT USED] Error variance threshold where if the variance
                of length scales is less than interqubit separation, then algorithm
                terminates.
        '''

        # Set up Data Qubits (qubits which are never physically measured):
        if list_of_dataqubits is None:
            list_of_dataqubits = self.GLOBALDICT["DATA_QUBITS"]

        # Set stopping protocol
        if max_num_iterations is None:
            max_num_iterations = self.MODELDESIGN["MAX_NUM_ITERATIONS"]

        # Create particles
        self.InitializeParticles()

        # Reset Fano factor to be of order R_Max.
        r_fano = np.ones(self.QubitGrid.number_of_nodes) * self.R_MAX
        self.QubitGrid.set_all_nodes("r_state_variance", r_fano)
        run_variance = self.QubitGrid.get_all_nodes(["r_state_variance"])

        # Initialize control strategy (automated or user driven)
        if autocontrol == "OFF" and self.measurements_controls is None:
            print "Auto-controller is off and no measurement control protocol is specified "
            raise RuntimeError

        if measurements_controls is not None:
            self.measurements_controls = measurements_controls
            PROTOCOL_ON = True

        if autocontrol == "ON":
            self.measurements_controls = [(0.0, 0.0)]
            PROTOCOL_ON = True

        max_iter_condition = max(len(self.measurements_controls), max_num_iterations)
        next_control_neighbourhood = range(0, self.QubitGrid.number_of_nodes)

        # Run QSLAM
        protocol_counter = 0
        while PROTOCOL_ON is True:

            msmt_control_pair = self.PerformMeasurement(autocontrol,
                                                        run_variance,
                                                        protocol_counter,
                                                        next_control_neighbourhood=next_control_neighbourhood,
                                                        list_of_dataqubits=list_of_dataqubits)

            next_control_neighbourhood = self.RunTwoLayerParticleFilter(msmt_control_pair)

            if self.save_run is True:
                self.save_run_variance.append(run_variance)
                self.save_controls_j.append(msmt_control_pair)

            if protocol_counter == max_iter_condition - 1:
                # print "PROTOCOL - SAFE END - Max number of measurements taken"
                PROTOCOL_ON = False

            run_variance = self.QubitGrid.get_all_nodes(["r_state_variance"])

            protocol_counter += 1


# #############################################################################
# CONFIGURE QUBIT GRID PARAMETERS
# #############################################################################

    def GetHardwareConfig(self):
        '''Return max qubit pairwise separation  on hardware, and lower and upper
        bounds on uniform distribution of r_states as a prior.

        Returns:
        -------
            dgrid (`float` | scalar):
                Max spatial pairwise separation for any two qubits on the Grid.
            d_iterq (`float` | scalar):
                Min spatial pairwise separation for any two qubits on the Grid.
            R_MAX (`float` | scalar):
                Upper bound on maximal correlation lengths physically expected
                in the system.
            R_MIN (`float` | scalar):
                Lower bound on maximal correlation lengths physically expected
                in the system.
        '''
        list_of_positions = self.QubitGrid.list_of_nodes_positions

        d_grid, _ = ParticleFilter.find_max_distance(list_of_positions)
        diterq, _ = ParticleFilter.find_min_distance(list_of_positions)

        R_MAX = d_grid * self.MODELDESIGN["MULTIPLER_R_MAX"]
        R_MIN = diterq * self.MODELDESIGN["MULTIPLER_R_MIN"]

        return d_grid, diterq, R_MIN, R_MAX

# #############################################################################
# PERFORM SIMULATED MEASUREMENTS
# #############################################################################

    def PerformMeasurement(self,
                           autocontrol,
                           listofcontrolparameters,
                           protocol_counter,
                           list_of_dataqubits=None,
                           next_control_neighbourhood=None):

        ''' Return measurement outcome and control directive i.e. the location of
            the qubit node that is measured.

        Parameters:
        ----------
            autocontrol_( "OFF" / "ON" flag) :
                "OFF" - next qubit measured is specified as a user input via measurement_controls_
                "ON" - next qubit measured is chosen by the algorithm.
            listofcontrolparameters (`float64`| numpy array):
                Avg. uncertainity metric for correlation lengthscales at each node.
                Dims: self.QubitGrid.number_of_nodes
            protocol_counter (`dtype`: int scalar):
                Tracks the number of iterations for QSLAM algorithm.
            next_control_neighbourhood (`int`| list) :
                List of indices for a qubit within a control region.
                Default value: NONE.
            list_of_dataqubits (`int`| list) :
                List of indices for qubit locations for data qubits that
                cannot be used for sensing measurements.
                Defaults to None.

        Returns:
        -------
            msmtset_j (`type`: list):
                List of 0 or 1 single qubit measurements at qubit location indexed
                by node_j.
            node_j (`dtype`: int scalar):
                Index of a qubit location on a 2D arrangement.
        '''

        if autocontrol == "OFF":
            return self.measurements_controls[protocol_counter]

        elif autocontrol == "ON":

            node_j = controller(listofcontrolparameters,
                                next_control_neighbourhood,
                                list_of_dataqubits=list_of_dataqubits,
                                number_of_nodes=1)[0]

            # TODO: adapt for arbitrary number of simultaneous msmts at different locations, number_of_nodes > 1
            msmtset_j = [self.QubitGrid.measure_node(node_j) for idx in range(self.MODELDESIGN["MSMTS_PER_NODE"])]

            return [msmtset_j, node_j]

# #############################################################################
# TWO LAYER PARTICLE FILTERING ALGORITHM
# #############################################################################

    def RunTwoLayerParticleFilter(self, msmt_control_pair):
        ''' Return next iteration of particle filtering algorithm after processing
            an input measurement outcome and control directive.

            Parameters:
            ----------
             msmt_control_pair (`type` | list):

                msmt_control_pair[0] contains a list of 0 or 1 single qubit
                measurements performed at qubit location indexed by msmt_control_pair[1]

                msmt_control_pair[1] is the qubit location index at which
                measurements in msmt_control_pair[0] were performed.

            Returns:
            -------
                next_control_neighbourhood (`type` | list):
                List of indices for qubit locations for input into the QSLAM controller.
        '''

        msmtset_j = msmt_control_pair[0]
        control_j = msmt_control_pair[1]

        for alpha_particle in self.AlphaSet.particles:
            alpha_particle.pset_beta = self.MODELDESIGN["P_BETA"]
            alpha_particle.node_j = control_j

        ###### BEGIN MSMT LOOP / ESTIMATE LOCALLY
        for next_phys_msmt_j in msmtset_j:

            self.ReceiveMsmt(control_j, next_phys_msmt_j)
            self.PropagateState(control_j)

            if self.save_run is True:
                self.save_alpha_predictive.append(self.AlphaSet.particles)

            self.ComputeAlphaWeights()

        ###### END  MSMT LOOP / ESTIMATE LOCALLY

        if self.save_run is True:
            self.save_alpha_bar_1.append(self.AlphaSet.particles)

        ###### COMPUTE JOINT DISTRIBUTION AND WEIGHTS
        joint_weights = self.ComputeJointWeights(control_j,
                                                 **self.LikelihoodObj.WEIGHTFUNCDICT_BETA)

        ###### RESAMPLE JOINT DISTRIBUTION, COLLAPSE BETA LAYER, RESAMPLE ALPHA PARTICLES
        self.ExecuteParticleBranching(joint_weights)

        ###### COMPUTE POSTERIOR STATE
        posterior_state = self.AlphaSet.posterior_state
        self.QubitGrid.state_vector = posterior_state*1.0

        if self.save_run is True: # MARKER: intermediary dist. save func. Jun-19
            self.save_alpha_posterior.append(self.AlphaSet.particles)
            self.save_posterior_state.append(posterior_state)

        ###### SHARE DATA MESSAGES, IDENTITY NEXT CONTROL NEIGHBOURHOOD
        next_control_neighbourhood = self.update_qubitgrid_via_quasimsmts(control_j,
                                                                          posterior_state)
        return next_control_neighbourhood # neighbourhood of next control action.


#       ------------------------------------------------------------------------
#       SUPPORT FUNCTION 1: INITIALISE AND SAMPLE AT t = 0
#       ------------------------------------------------------------------------

    def InitializeParticles(self):
        '''Intialise properties of QubitGrid and AlphaSet for t = 0 time step.
        Note f_state cannot be set, merely read from QubitGrid.
        '''
        for alphaparticle in self.AlphaSet.particles:
            self.set_init_alphaparticle(alphaparticle)

        self.QubitGrid.state_vector = self.AlphaSet.posterior_state


    def set_init_alphaparticle(self, alphaparticle):
        ''' Set initial properties of an Alpha particle.
            Helper function to InitializeParticles().
        '''

        self.PRIORDICT["SAMPLE_R"]["ARGS"]["R_MIN"] = self.R_MIN
        self.PRIORDICT["SAMPLE_R"]["ARGS"]["R_MAX"] = self.R_MAX
        size = self.QubitGrid.number_of_nodes

        substate_list = []
        for idx_key in ["SAMPLE_X", "SAMPLE_Y", "SAMPLE_F", "SAMPLE_R"]:

            # f_state cannot be set, merely read from QubitGrid

            if idx_key == "SAMPLE_F":

                # # All alpha particles have the same prior ie. same f_state values.
                # samples = self.QubitGrid.get_all_nodes(["f_state"])

                # MARKER JUNE 2019 - differentiate alpha maps
                # Different alpha particle initial state
                self.PRIORDICT["SAMPLE_F"]["ARGS"]["SIZE"] = size
                samples = self.PRIORDICT["SAMPLE_F"]["FUNCTION"](**self.PRIORDICT["SAMPLE_F"]["ARGS"])

            if idx_key != "SAMPLE_F":

                self.PRIORDICT[idx_key]["ARGS"]["SIZE"] = size
                samples = self.PRIORDICT[idx_key]["FUNCTION"](**self.PRIORDICT[idx_key]["ARGS"])

                if idx_key == "SAMPLE_X":
                    samples += self.QubitGrid.get_all_nodes(["x_state"])

                if idx_key == "SAMPLE_Y":
                    samples += self.QubitGrid.get_all_nodes(["y_state"])

            substate_list.append(samples)

        alphaparticle.particle = np.asarray(substate_list).flatten()
        alphaparticle.pset_beta = self.MODELDESIGN["P_BETA"]

#       ------------------------------------------------------------------------
#       SUPPORT FUNCTION 2: RECEIVE MSMT
#       ------------------------------------------------------------------------

    def ReceiveMsmt(self, control_j, next_phys_msmt_j):
        ''' Updates global qubit hardware variables and properties of Alpha
            particles for a new measurement outcome and control directive.
        '''
        # Update tau counter before likelihood calculation
        # self.QubitGrid.nodes[control_j].physcmsmtsum = next_phys_msmt_j
        # prob_j = self.QubitGrid.nodes[control_j].sample_prob_from_msmts()

        # MARKER JUN 2019 - Update tau counter after likelihood calculation
        prob_j = self.QubitGrid.nodes[control_j].sample_prob_from_msmts()
        if prob_j is None:
            initial_state = self.QubitGrid.nodes[control_j].f_state
            prob_j = self.QubitGrid.nodes[control_j].born_rule(initial_state)

        pl.update_alpha_dictionary(next_phys_msmt_j,
                                   prob_j,
                                   **self.LikelihoodObj.LIKELIHOOD_ALPHA)

        # MARKER JUN 2019 - Update tau counter after likelihood calculation
        self.QubitGrid.nodes[control_j].physcmsmtsum = next_phys_msmt_j

#       ------------------------------------------------------------------------
#       SUPPORT FUNCTION 3: PROPAGATE (ALPHA) STATES
#       ------------------------------------------------------------------------

    def PropagateState(self, control_j):
        '''Propagates all Alpha particle states according to (apriori) transition
        propability distributions for propagating states between measurements.
        '''
        for alpha_particle in self.AlphaSet.particles:
            self.sample_from_transition_dist(alpha_particle, control_j)


    def sample_from_transition_dist(self, alpha_particle, control_j):
        ''' Helper function for PropagateState(). Under identity dynamics (time invariance)
        this function does nothing in qslamr.
        '''
        # Implements time invariant states
        # Currently implements identity (do nothing).
        # TODO: Placehodler for time-invariance


#       ------------------------------------------------------------------------
#       SUPPORT FUNCTION 4: COMPUTE ALPHA WEIGHTS; GENERATE BETA WEIGHTS
#       ------------------------------------------------------------------------
    def ComputeAlphaWeights(self):
        ''' Update weights for the distribution defined over Alpha particles.
        '''
        new_alpha_weights = self.AlphaSet.calc_weights_set() # Normalised

        self.AlphaSet.weights_set = new_alpha_weights

    def ComputeJointWeights(self, control_j, **BETADICT):
        ''' Return posterior weights for the joint distribution defined over both
        Alpha and Beta particles.
        '''

        f_state_index = 2*self.QubitGrid.number_of_nodes + control_j

        posterior_weights = []

        for idx_alpha in range(self.MODELDESIGN["P_ALPHA"]):
            alpha_particle = self.AlphaSet.particles[idx_alpha]

            alpha_particle.particle[f_state_index] = self.QubitGrid.nodes[control_j].f_state

            beta_alpha_j_weights = self.generate_beta_layer(alpha_particle, **BETADICT)

            if self.save_run is True: # MARKER: intermediary dist. save func. Jun-19
                self.save_all_beta_weights.append(beta_alpha_j_weights)

            posterior_weights.append(alpha_particle.weight*beta_alpha_j_weights)

        if self.save_run is True: # MARKER: intermediary dist. save func. Jun-19
            self.save_alpha_beta_joint.append(self.AlphaSet.particles)

        posterior_weights = np.asarray(posterior_weights).flatten()
        normalisation = np.sum(posterior_weights)

        if normalisation == 0.0:
            # print "Zero value normalisation in ComputeJointWeights()"
            # print "Resetting to uniform weights"
            normalised_posterior_weights = self.posterior_reset()

            if self.save_run is True:
                self.save_alpha_beta_joint_weights.append(normalised_posterior_weights)

            return normalised_posterior_weights

        normalised_posterior_weights = posterior_weights*(1.0/normalisation)

        if np.any(np.isnan(normalised_posterior_weights)):
            # print "Invalid Nan values encountered in ComputeJointWeights()"
            # print "Resetting to uniform weights"
            normalised_posterior_weights = self.posterior_reset()

            if self.save_run is True:
                self.save_alpha_beta_joint_weights.append(normalised_posterior_weights)

            return normalised_posterior_weights

        # print "ComputeJointWeights() has no error - yay!"

        if self.save_run is True: # MARKER: intermediary dist. save func. Jun-19
            self.save_alpha_beta_joint_weights.append(normalised_posterior_weights)

        return  normalised_posterior_weights


    def posterior_reset(self):
        '''Return a uniform distribution following a singular or invalid
        numerical posterior distribution'''

        dimension = float(self.MODELDESIGN["P_ALPHA"]) * float(self.MODELDESIGN["P_BETA"])
        normalised_posterior_weights = (1.0 / dimension)* np.ones(int(dimension))
        return normalised_posterior_weights


    def generate_beta_layer(self, alpha_particle, **BETADICT):
        ''' Return a set of Beta particle weights for a given Alpha parent. Helper
            function to ComputeJointWeights().
        '''
        len_idx = self.QubitGrid.number_of_nodes*3 + alpha_particle.node_j
        parent_alpha = alpha_particle.particle.copy()
        new_beta_state = parent_alpha * 1.0 # get property

        list_of_parent_states = []
        list_of_length_samples = []

        for idx_beta in range(alpha_particle.pset_beta):
            new_length_sample = self.sample_radii(control_j=alpha_particle.node_j) # MARKER
            new_beta_state[len_idx] = new_length_sample*1.0
            list_of_parent_states.append(new_beta_state.copy())
            list_of_length_samples.append(new_length_sample)

        alpha_particle.generate_beta_pset(list_of_parent_states,
                                          list_of_length_samples,
                                          **BETADICT)
        for beta_particle_object in alpha_particle.BetaAlphaSet_j.particles:

            self.generate_beta_neighbourhood(beta_particle_object)

        beta_alpha_j_weights = alpha_particle.BetaAlphaSet_j.calc_weights_set()

        return beta_alpha_j_weights


    def sample_radii(self, previous_length_scale=0.0, control_j=None): # MARKER
        ''' Return a r_state sample to generate a new Beta particle for a given Alpha parent.
        Helper function to generate_beta_layer().
        '''

        if self.beta_expansion_mode: # MARKER
            # New trial

            mean = self.QubitGrid.get_all_nodes(["r_state"])[control_j]
            fano = self.QubitGrid.get_all_nodes(["r_state_variance"])[control_j]
            var = fano * mean
            sd = np.sqrt(var)

            if sd == 0:
                sd = np.sqrt(self.GLOBALDICT["NOISEPARAMS"]["SIGMOID_APPROX_ERROR"]["SIGMA"]) * 10.**-6
                # set very small rel to likelhood parameter variances but not zero

            a, b = (self.R_MIN - mean) / sd, (self.R_MAX - mean) / sd

            valid_sample = False
            counter = 0

            while valid_sample is False:
                sample = truncnorm.rvs(a, b, loc=mean, scale=sd, size=1)
                counter += 1

                if sample >= 0:
                    valid_sample = True

                    if counter > 1:
                        print "Repeat samples req: ", counter

            return sample

        if not self.beta_expansion_mode:

            # Pre June 2019 method: Uniformly injected samples betwn Rmin and Rmax

            if previous_length_scale < 0:
                print "Previous length scale is less than zero:", previous_length_scale
                raise RuntimeError
            lower_bound = (previous_length_scale + self.R_MIN)*0.1 + 0.9*self.R_MIN
            sample = np.random.uniform(low=lower_bound, high=self.R_MAX)

            return sample


    def generate_beta_neighbourhood(self, BetaParticle):
        ''' Generate phase estimates over a neighbourhood for input Beta particle
            and its Alpha parent. Helper function to generate_beta_layer().
         '''

        NEIGHBOURDICT = {"prev_posterior_f_state" : self.QubitGrid.get_all_nodes(["f_state"]),
                         "prev_counter_tau_state" : self.QubitGrid.get_all_nodes(["counter_tau"]),
                         "lambda_" : self.MODELDESIGN["LAMBDA_2"],
                         "kernel_function": self.MODELDESIGN["kernel_function"]}

        BetaParticle.smear_fj_on_neighbours(**NEIGHBOURDICT)

#       ------------------------------------------------------------------------
#       SUPPORT FUNCTION 5: RESAMPLE AND BETA COLLAPSE
#       ------------------------------------------------------------------------

    def ExecuteParticleBranching(self, joint_weights):
        '''Return a new set of Alpha particles with uniform weights resampled
        according to the joint probability of both Alpha and Beta particles, subject to
        conserving Alpha particle number.
        '''

        ##### RESAMPLE JOINT DISTRIBUTION
        resampled_idx = ParticleFilter.resample_constant_pset_alpha(joint_weights,
                                                                    self.MODELDESIGN["P_ALPHA"],
                                                                    self.MODELDESIGN["P_BETA"])

        new_alpha_subtrees = self.get_subtrees(resampled_idx, self.MODELDESIGN["P_BETA"])

        ##### COLLAPSE BETA LAYER AND RESAMPLE ALPHAS
        new_alpha_list = self.collapse_beta(new_alpha_subtrees, resampled_idx)

        ##### SET POSTERIOR ALPHA DISTRIBUTION
        self.AlphaSet.particles = new_alpha_list

        ##### RESET POSTERIOR TO UNIFORM WEIGHTS
        uniformprob = (1.0/self.MODELDESIGN["P_ALPHA"])
        self.AlphaSet.weights_set = uniformprob*np.ones(self.MODELDESIGN["P_ALPHA"])


    def collapse_beta(self, subtree_list, resampled_indices):
        ''' Return resampled Alpha parents with collapsed Beta layers. For each Alpha
        parent, store the mode and the spread of r_states in Beta-Particle layer as
        the posterior r_state at node_j.  Helper function for ExecuteParticleBranching().

        Parameters:
        ----------
            subtree_list (`type`| list of pairs) :
                List of  pairs of (start, end) for a sub-tree, where each sub-tree
                represents an Alpha particle, and its element leaves represent suriviving
                resampled Beta particles.
            resampled_indices (`type`| list of int) :
                List of indices formed by resampling joint distribution of alpha-beta
                particles.

        Returns:
        -------
            final_alpha_list (`type`| list of alpha particle objects):
                List of Alpha particles that constitute the posterior particle
                distribution for QSLAM
                Dims: self.MODELDESIGN["P_ALPHA"]
        '''
        state_update = 0.
        new_alpha_particle_list = []
        leaf_count_list = []
        alpha_node_list = [] # MARKER
        normaliser = (1./float(len(subtree_list))) # Num of. surviving alpha particles
        uncertainity_at_j = 0.0

        ##### ANALYSE BETA LAYER
        # Each subtree in  subtree_list represents an uneven alpha-beta tree
        # after re-sampling of the joint distribution of alpha and beta particles.

        for subtree in subtree_list:

            leaves_of_subtree = resampled_indices[subtree[0]:subtree[1]]

            # leaf_count: number of surviving beta particles for a surviving alpha parent
            leaf_count = float(len(leaves_of_subtree))

            if leaf_count != 0: # excludees the [0,0] subtree

                # Extract r_state uncertainty measures
                alpha_node = ParticleFilter.get_alpha_node_from_treeleaf(leaves_of_subtree[0], self.MODELDESIGN["P_BETA"])
                beta_alpha_nodes = [ParticleFilter.get_beta_node_from_treeleaf(leafy, self.MODELDESIGN["P_BETA"]) for leafy in leaves_of_subtree]
                r_est_index = self.QubitGrid.number_of_nodes*3 + self.AlphaSet.particles[alpha_node].node_j
                r_est_subtree_list = []

                for node in beta_alpha_nodes:
                    beta_state = self.AlphaSet.particles[alpha_node].BetaAlphaSet_j.particles[node].particle.copy()
                    beta_lengthscale = beta_state[r_est_index]*1.0
                    if np.isnan(beta_lengthscale):
                        raise RuntimeError
                    r_est_subtree_list.append(beta_lengthscale)

                # Store r_state uncertainty measures
                r_mean, r_var, r_fano = ParticleFilter.calc_posterior_lengthscale(np.asarray(r_est_subtree_list), skew_adjust=self.skew_adjust)
                parent = self.AlphaSet.particles[alpha_node].particle.copy()*1.0
                parent[r_est_index] = r_mean # r_mean could be mode depending on skew

                if np.any(np.isnan(parent)):
                    print "A resampled parent particle has an invalid value."
                    raise RuntimeError

                uncertainity_at_j += r_fano * normaliser

                # Update Alpha particle for r_state mean
                self.AlphaSet.particles[alpha_node].particle = parent*1.0

                # Store marginalised alpha particles and their weights
                new_alpha_particle_list.append(self.AlphaSet.particles[alpha_node])
                leaf_count_list.append(leaf_count)

                if self.save_run is True:
                    alpha_node_list.append(alpha_node)


        self.QubitGrid.nodes[self.AlphaSet.particles[alpha_node].node_j].r_state_variance = uncertainity_at_j
        # NB: Alpha_node is just pulling a location

        ##### RESAMPLE ALPHA PARTICLES POST MARGINALISATION OF BETA LAYER
        number_of_new_alphas = len(new_alpha_particle_list)
        marginalise_weights = leaf_count_list / np.sum(leaf_count_list)

        if self.save_run is True:
            self.save_alpha_bar_2.append(new_alpha_particle_list)
            self.save_alpha_bar_2_weights.append(marginalise_weights)
            self.save_alpha_bar_2_oldnodes.append(alpha_node_list)

        ##### GET POSTERIOR ALPHA LAYER
        p_alpha_indices = ParticleFilter.resample_from_weights(marginalise_weights,
                                                               self.MODELDESIGN["P_ALPHA"])

        final_alpha_list = [new_alpha_particle_list[idx_] for idx_ in p_alpha_indices]

        ##### RESET BETA LAYER FOR POSTERIOR ALPHA PARTICLES
        self.null_beta_layer(final_alpha_list)

        return final_alpha_list

    def null_beta_layer(self, list_of_alpha_particles):
        ''' Strips input list of Alpha particles of its individual Beta layers.
        Helper function to self.collapse_beta()
        '''

        for alpha_idx in list_of_alpha_particles:
            alpha_idx.BetaAlphaSet_j = None

    def effective_particle_size(self, posterior_weights):
        '''Return effective particle number based on posterior weight variances
        for adapative resampling. [NOT USED].
        '''

        p_size = 1.0/ np.sum(posterior_weights**2)
        return p_size

#       ------------------------------------------------------------------------
#       SUPPORT FUNCTIONS 6: SHARING MECHANISM AND NEXT CONTROL STEP
#       ------------------------------------------------------------------------

    def update_qubitgrid_via_quasimsmts(self, control_j, posterior_state):
        ''' Return the posterior neighbourhood associated with a qubit measurement
            at node_j and update all neighbours of node_j with data messages
            from the posterior state estimates at node_j.

        Parameters:
        ----------
            control_j (`dtype`: float)
                Location index of measured qubit at iteration t.
            posterior_state (`dtype`: numpy array)
                Posterior state for QSLAM at the end of iteration t.

        Returns:
        -------

        CRITICAL -- this function should only be applied after lengthscales have
            been discovered (alpha, beta_alpha) particles carried over or are
            resampled with sufficient alpha diversity; weights are set to uniform.
        '''
        #### print "In sharing mechanism"

        posterior_radius = posterior_state[self.QubitGrid.number_of_nodes*3 + control_j]
        posterior_beta_particle = BetaParticle(control_j, posterior_state, posterior_radius)
        self.generate_beta_neighbourhood(posterior_beta_particle)

        # smear_fj_on_neighbours() is being called somehow between these print statements

        #### print "Has smear phases been called?" ## Yes it has

        for idx in range(len(posterior_beta_particle.neighbourhood_qj)):

            neighbour_q = posterior_beta_particle.neighbourhood_qj[idx]
            quasi_phase_q = posterior_beta_particle.smeared_phases_qj[idx]

#            print quasi_phase_q

            if quasi_phase_q >= 0.0 and quasi_phase_q <= np.pi:
                born_prob_q = Node.born_rule(quasi_phase_q)
                quasi_msmt = Node.quantiser(born_prob_q)
                self.QubitGrid.nodes[neighbour_q].quasimsmtsum = quasi_msmt

            elif quasi_phase_q < 0.0 or quasi_phase_q > np.pi:
                print "quasi-phase posterior at q=", neighbour_q
                print "...was invalid, q_phase ", quasi_phase_q
                print "... no quasi_msmts were added."

        return posterior_beta_particle.neighbourhood_qj

#       ------------------------------------------------------------------------
#       SUPPORT FUNCTIONS 7: STATIC METHODS
#       ------------------------------------------------------------------------

    @staticmethod
    def calc_posterior_lengthscale(r_lengthscales_array, skew_adjust=True):
        '''Return descriptive moments for the distribution of r_states in a
        Beta particle set for an Alpha parent, post resampling.
        Helper function for collapse_beta().

        Parameters:
        ----------
            r_lengthscales_array (`dtype`| numpy array):
                A 1D array of r_state lengthscales at a single qubit location.
                Typically implied by a distribution of suruving beta particles
                for an alpha parent after the joint-distribution has been resampled.

            skew_adjust (`type`| Boolean):
                If FALSE, r_mean is the mean of r_lengthscales_array.
                If TRUE, r_mean is the mean of r_lengthscales_array when the
                distribution of r_lengthscales_array is not skewed; and r_mean
                is the mode if the distribution of r_lengthscales_array is skewed left
                or right.

        Returns:
        -------
            r_mean ('dtype' | float) : Empirical mean of r_lengthscales_array
            r_var ('dtype' | float) : Empirical variance of r_lengthscales_array
            r_fano ('dtype' | float) : Empirical fano factor of r_lengthscales_array

        '''
        r_mean, r_var, r_fano = ParticleFilter.calc_skew(r_lengthscales_array, skew_adjust)

        return r_mean, r_var, r_fano

    @staticmethod
    def calc_skew(r_lengthscales_array, skew_adjust=True):
        '''Helper function for  calc_posterior_lengthscale().
        Parameters:
        ----------
            r_lengthscales_array (`dtype`| numpy array):
                A 1D array of r_state lengthscales at a single qubit location.
                Typically implied by a distribution of suruving beta particles
                for an alpha parent after the joint-distribution has been resampled.

            skew_adjust (`type`| Boolean):
                If FALSE, r_mean is the mean of r_lengthscales_array.
                If TRUE, r_mean is the mean of r_lengthscales_array when the
                distribution of r_lengthscales_array is not skewed; and r_mean
                is the mode if the distribution of r_lengthscales_array is skewed left
                or right.

        Returns:
        -------
            r_mean ('dtype' | float) : Empirical mean of r_lengthscales_array
            r_var ('dtype' | float) : Empirical variance of r_lengthscales_array
            r_fano ('dtype' | float) : Empirical fano factor of r_lengthscales_array
        '''
        totalcounts = len(r_lengthscales_array)
        mean_ = np.mean(r_lengthscales_array)
        mode_, counts = mode(r_lengthscales_array)
        median_ = np.sort(r_lengthscales_array)[int(totalcounts/2) - 1]

        variance = np.var(r_lengthscales_array)
        fano = variance / mean_

        if skew_adjust is True: # MARKER changed July 2019 to output fano

            if mean_ < mode_ and counts > 1:
                return mode_, variance, fano

            if mean_ > mode_ and counts > 1:
                return mode_, variance, fano

        return mean_, variance, fano


    @staticmethod
    def compute_dist(one_pair):
        '''Return Euclidean distance between two points.
        Helper function for ParticleFilter.get_distances().'''
        xval, yval = one_pair
        return np.sqrt((xval[0] - yval[0])**2 + (xval[1] - yval[1])**2)

    @staticmethod
    def get_distances(list_of_positions):
        '''Return all pairwise distances between an arrangement of qubits.
        Helper function for ParticleFilter.find_max_distance() and
        ParticleFilter.find_min_distance.
        '''
        distance_pairs = [a_pair for a_pair in combinations(list_of_positions, 2)]
        distances = [ParticleFilter.compute_dist(one_pair)for one_pair in distance_pairs]
        return distances, distance_pairs

    @staticmethod
    def find_max_distance(list_of_positions):
        '''Return the pair of positions with maximum pair-wise separation distance. '''
        distances, distance_pairs = ParticleFilter.get_distances(list_of_positions)
        return max(distances), distance_pairs[np.argmax(distances)]

    @staticmethod
    def find_min_distance(list_of_positions):
        '''Return  the pair of positions with minimum pair-wise separation distance. '''
        distances, distance_pairs = ParticleFilter.get_distances(list_of_positions)
        return min(distances), distance_pairs[np.argmin(distances)]

    @staticmethod
    def get_alpha_node_from_treeleaf(leaf_index, pset_beta):
        '''Return Alpha Particle index based on global index for flattened layers
        of Alpha-Beta particles. Helper function for collapse_beta().

        Parameters:
        ----------
            leaf_index  (`dtype`| float):
                Index representing an (alpha, beta) pair in the joint distribution
                after flattening both particle layers.
            pset_beta (`dtype`| int):
                Total number of Beta particles for each Alpha parent.

        Returns:
        -------
            alpha_node (`dtype`| int):
                Index representing the Alpha particle parent of the
                Beta-Alpha pair represented by leaf_index.
        '''
        alpha_node = int(leaf_index//float(pset_beta))
        return alpha_node

    @staticmethod
    def get_beta_node_from_treeleaf(leaf_index, pset_beta):
        '''Return Beta Particle index based on global index for flattened layers
        of Alpha-Beta particles. Helper function for collapse_beta().

        Parameters:
        ----------
            leaf_index  (`dtype`| float):
                Index representing an (alpha, beta) pair in the joint distribution
                after flattening both particle layers.
            pset_beta (`dtype`| int):
                Total number of Beta particles for each Alpha parent.

        Returns:
        -------
            beta_node (`dtype`| int):
                Index representing the unflattened Beta particle index for an Alpha parent,
                corressponding to the Beta-Alpha pair represented by leaf_index.
        '''

        beta_node = int(leaf_index - int(leaf_index//float(pset_beta))*pset_beta)
        return beta_node

    @staticmethod
    def resample_from_weights(posterior_weights, number_of_samples):
        ''' Return samples from the posterior distribution of states approximated
        by a discrete set of posterior normalised weights.
        Helper function for static method ParticleFilter.resample_constant_pset_alpha().

        Parameters:
        ----------
            posterior_weights (`dtype` | numpy array) :
                A set of weights that form a discrete probability
                distribution. All weights >= 0 and sum to 1.

            number_of_samples (`type` | int):
                Total number of random variates to be acquired as samples from
                the discrete probability distribution represented as posterior_weights.
        Returns:
        -------
            resampled_idx (`type` | list):
                List of random variates sampled from the discrete probability distribution,
                represented as indicies of the x_axis.
                Dims: number_of_samples

        '''
        total_particles = len(posterior_weights)
        cdf_weights = np.asarray([0] + [np.sum(posterior_weights[:idx+1]) for idx in range(total_particles)])
        pdf_uniform = np.random.uniform(low=0, high=1.0, size=number_of_samples)

        resampled_idx = []
        for u_0 in pdf_uniform:
            j = 0
            while u_0 > cdf_weights[j]:
                j += 1
                if j > total_particles:
                    j = total_particles
                    break
            resampled_idx.append(j-1)
        return resampled_idx

    @staticmethod
    def resample_constant_pset_alpha(posterior_weights, pset_alpha, pset_beta):
        ''' Return indicies for Alpha particles resampled from posterior Alpha-Beta weights, while
        conserving Alpha particle number. Helper function for ExecuteParticleBranching().

        Parameters:
        ----------
            posterior_weights (`dtype` | numpy array) :
                A set of weights that form a discrete probability
                distribution. All weights >= 0 and sum to 1.
            pset_alpha (`int` | scalar):
                Total number of of Alpha particles (conserved during re-sampling).
            pset_beta (`int` | scalar):
                Total number of of Beta particles for each Alpha parent.

        Returns:
        -------
            resampled_idx (`type` | list):
                List of random variates sampled from the discrete probability distribution,
                represented as indicies of the x_axis.
                Dims: number_of_samples (== pset_alpha)
        '''

        num_of_samples = pset_alpha
        # change to: num_of_samples = pset_alpha * pset_beta
        # # no don't change it! convergence proof simplier if total trials doesn't change
        resampled_indices = ParticleFilter.resample_from_weights(posterior_weights,
                                                                 num_of_samples)
        return resampled_indices

    @staticmethod
    def get_subtrees(resampled_indices, pset_beta):
        '''
        Return a list of  pairs of (start, end) to sub-divide a sorted list of
        resampled_idx. Each subset represents a sub-tree, constituting of a surviving
        Alpha parent, and its suriviving Beta particles (leaves) after resampling the
        joint distribution over both Alpha and Beta layers.
        Helper function for ExecuteParticleBranching.

        Parameters:
        ----------
            resampled_idx (`type` | list):
                List of random variates sampled from the discrete probability distribution,
                represented as indicies of the x_axis.

            pset_beta (`int` | scalar):
                Total number of of Beta particles for each Alpha parent.
        Returns:
        -------
            new_sub_trees (`type` | list):
                A list of  pairs of (start, end) for a sub-tree, where each sub-tree
                is a sub-set of a sorted resampled_idx list. The subtree represents an Alpha parent,
                and its elements represent suriviving resampled Beta particles.
        '''
        new_sub_trees = []

        resampled_indices.sort()
        alpha_index_0 = None
        strt_counter = 0
        end_counter = 0

        for idx in resampled_indices:

            alpha_index = ParticleFilter.get_alpha_node_from_treeleaf(idx, pset_beta)
            beta_alpha_idx = ParticleFilter.get_beta_node_from_treeleaf(idx, pset_beta)

            if alpha_index_0 == alpha_index:
                end_counter += 1

            elif alpha_index_0 != alpha_index:
                new_sub_trees.append([strt_counter, end_counter])

                alpha_index_0 = alpha_index
                strt_counter = end_counter

                end_counter += 1

        if end_counter == len(resampled_indices):
            new_sub_trees.append([strt_counter, end_counter])

        return new_sub_trees
