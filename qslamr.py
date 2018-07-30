'''
MODULE: qslamr

Returns discrete posterior for the SLAM problem for qubit control as a set of
particles (true states) and their weights (posterior distribution).

CLASS: Node, Grid

METHODS

PARAMETERS

STATIC FUNCTIONS


'''
import numpy as np
import particledict as pd
from model_design import GRIDDICT, INITIALDICT
from hardware import Grid, PARTICLE_STATE
from itertools import combinations
from particlesets import AlphaParticle, ParticleSet, BetaParticle
from particledict import WEIGHTFUNCDICT_ALPHA, WEIGHTFUNCDICT_BETA, update_alpha_dictionary
from scipy.stats import mode

class ParticleFilter(Grid):
    '''doctring
    
    MU_W : 0.0, # Qubit position noise mean (dynamics)
    SIG2_W : 1.0, # Qubit position noise variance (dynamics)
    MU_R : 0.0, # Length scale noise mean (dynamics)
    SIG2_R : 1.0, # Length scale noise variance (dynamics)
    MU_MEASR : 0.0, # Map noise mean (measurement)
    SIG2_MEASR : 1.0, # Map noise variance (measurement)
    MU_F : 0.0, # True sigmoid approximation error mean
    SIG2_F : 1.0, # True sigmoid approximation error variance
    LAMBDA : 0.9, # Forgetting factor for quasi-msmt information
    GAMMA_T : 10**8, # Re-sampling threshold
    P_ALPHA : 10, # Number of alpha particles
    P_BETA : 10, # Numer of beta particles for each alpha

    '''

    def __init__(self, list_of_nodes_positions, **INITIALDICT):

        self.QubitGrid = Grid(list_of_nodes_positions=list_of_nodes_positions)
        self.dgrid, self.R_min, self.R_max = ParticleFilter.set_uniform_prior_for_correlation(self.QubitGrid.list_of_nodes_positions)
        # self.dgrid, _ = self.find_max_distance(self.QubitGrid.list_of_nodes_positions)
        self.INITIALDICT = INITIALDICT


        self.pset_alpha = self.INITIALDICT["P_ALPHA"]
        self.pset_beta = self.INITIALDICT["P_BETA"]
        empty_alpha_particles = [AlphaParticle() for idx in range(self.pset_alpha)]
        self.AlphaSet = ParticleSet(empty_alpha_particles, **WEIGHTFUNCDICT_ALPHA)

        self.resample_thresh = INITIALDICT["GAMMA_T"]
        # self.R_max = self.dgrid * 0.5 # change this for min distance btwn qubits
        # self.R_min= self.dgrid * 1.5
        # print self.R_max


    def qslamr(self, measurements_controls):

        self.InitializeParticles()

        for t_item in measurements_controls:

            next_phys_msmt_j = t_item[0]
            control_j = t_item[1]
            self.ReceiveMsmt(control_j, next_phys_msmt_j)
            self.PropagateState(control_j)
            posterior_weights = self.ComputeWeights(control_j)
            self.ResampleParticles(posterior_weights)
            posterior_state = self.AlphaSet.posterior_state
            self.QubitGrid.state_vector = posterior_state.copy()*1.0
            self.update_qubitgrid_via_quasimsmts(control_j, posterior_state)

#       ------------------------------------------------------------------------
#       DATA ASSOCIATION FUNCTIONS
#       ------------------------------------------------------------------------

    def update_alpha_map_via_born_rule(self, alpha_particle): # TODO: Data association
        '''doc string
        this function must always come after we havw computed alpha weights
        '''

        for alpha_particle in self.AlphaSet.particles:

            born_prob = self.sample_prob_from_msmts(alpha_particle.node_j)

            if np.isnan(born_prob): 
                print "A born_prob has NAN value"
                raise RuntimeError

            if born_prob > 1:
                print "A born_prob >1 value"
                raise RuntimeError
            
            if born_prob < 0:
                print "A born_prob <0 value"
                raise RuntimeError
            
            if born_prob is None:
                print "No quasi msmts or data msmts available. Something is wrong"
                raise RuntimeError
            
            map_val = ParticleFilter.inverse_born(born_prob)
            
            if map_val > np.pi or np.isnan(map_val): 
                        print "A map_val has an invalid value ", map_val
                        print "the born probability was", born_prob
                        raise RuntimeError
            
            map_idx = self.QubitGrid.number_of_nodes * 2 + alpha_particle.node_j
            parent_particle = alpha_particle.particle.copy()*1.0 # get property
            parent_particle[map_idx] = map_val
            alpha_particle.particle = parent_particle # assign updated property


    def sample_prob_from_msmts(self, control_j): # TODO Data Association
        '''docstring'''
        
        forgetting_factor = self.INITIALDICT["LAMBDA"]**self.QubitGrid.nodes[control_j].counter_tau

        w_q = None
        prob_q = 0.0
        if self.QubitGrid.nodes[control_j].counter_beta != 0:
            prob_q = self.QubitGrid.nodes[control_j].quasimsmtsum / self.QubitGrid.nodes[control_j].counter_beta*1.0
            w_q = 0.5*forgetting_factor

        w_p = None
        prob_p =0.0
        if self.QubitGrid.nodes[control_j].counter_tau !=0:
            prob_p = self.QubitGrid.nodes[control_j].physcmsmtsum / self.QubitGrid.nodes[control_j].counter_tau*1.0
            w_p = 0.5 + 0.5*(1 - forgetting_factor)

        if w_p is None and w_q is None:
            print "NONE returned in sample_prob_from_msmts"
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
        '''docstring'''
        born_prob = np.cos(map_val / 2.0)**2
        return born_prob


    @staticmethod
    def inverse_born(born_prob):
        map_val = np.arccos(2.0*born_prob  - 1.0)
        return map_val

#       ------------------------------------------------------------------------
#       SMEARING ACTION VIA QUASI MEASUREMENTS
#       ------------------------------------------------------------------------
    
    def update_qubitgrid_via_quasimsmts(self, control_j, posterior_state):
        '''docstring
        this funciton should only be applied after lengthscales have been discovered
        (alpha, beta_alpha) particles carried over or are resampled with 
        sufficient alpha diversity; weights are set to uniform. 
        '''
        
        # print "In update qubitgrid via quasi msmts, posterior state = "
        # print posterior_state

        posterior_radius = posterior_state[self.QubitGrid.number_of_nodes*3 + control_j]
        SmearParticle = BetaParticle(control_j, posterior_state, posterior_radius)
        self.generate_beta_neighbourhood(SmearParticle)

        for idx in range(len(SmearParticle.neighbourhood_qj)):

            neighbour_q = SmearParticle.neighbourhood_qj[idx]
            quasi_phase_q = SmearParticle.smeared_phases_qj[idx]

            if quasi_phase_q >= 0.0 and quasi_phase_q <= np.pi:
                born_prob_q = ParticleFilter.born_rule(quasi_phase_q)
                quasi_msmt = np.random.binomial(1, born_prob_q)
                self.QubitGrid.nodes[neighbour_q].quasimsmtsum = quasi_msmt

            elif quasi_phase_q < 0.0 or quasi_phase_q > np.pi:
                print "quasi-phase posterior at q=", neighbour_q
                print "...was invalid, q_phase ", quasi_phase_q
                print "... no quasi_msmts were added."

#       ------------------------------------------------------------------------
#       SUPPORT FUNCTION 1: INITIALISE AND SAMPLE AT t = 0
#       ------------------------------------------------------------------------

    def InitializeParticles(self):
        '''Intialise properties of QubitGrid and AlphaSet for t = 0 time step'''

        for alphaparticle in self.AlphaSet.particles:
            self.set_init_alphaparticle(alphaparticle)

        self.QubitGrid.state_vector = self.AlphaSet.posterior_state


    def set_init_alphaparticle(self, alphaparticle):
        '''docstring'''

        sample_x = self.QubitGrid.get_all_nodes(["x_state"])
        sample_y = self.QubitGrid.get_all_nodes(["y_state"])
        #sample_s = np.random.normal(loc=self.INITIALDICT["MU_W"],
         #                           scale=self.INITIALDICT["SIG2_W"],
          #                          size=self.QubitGrid.number_of_nodes*2)

        sample_f = self.QubitGrid.get_all_nodes(["f_state"]) # np.random.uniform(low=0.,
                                     # high=np.pi, # INITIAL COND
                                     # size=self.QubitGrid.number_of_nodes)

        # sample_r = np.random.uniform(low=self.R_min,
        #                              high=self.R_max, # INITIAL COND
        #                              size=self.QubitGrid.number_of_nodes)
        # sample_r = self.QubitGrid.get_all_nodes(["r_state"]) # np.ones(self.QubitGrid.number_of_nodes)*ParticleFilter.sample_radii(self, self.R_min)
        sample_r = np.ones(self.QubitGrid.number_of_nodes)*ParticleFilter.sample_radii(self, self.R_min)
        alphaparticle.particle = np.concatenate((sample_x, sample_y, sample_f, sample_r),
                                                axis=0)

        # print "sample _r in Initialise particles", sample_r
        alphaparticle.pset_beta = self.INITIALDICT["P_BETA"]
        # alphaparticle.mean_radius_j = self.dgrid*1.5 # this should be seet as the posterior mean or mode of beta particles
        # TODO:1: alphaparticle.SIG2_MEASR = np.random.uniform(low=0.000000001, high=(0.1*np.pi)**2)

#       ------------------------------------------------------------------------
#       SUPPORT FUNCTION 2: RECEIVE MSMT
#       ------------------------------------------------------------------------

    def ReceiveMsmt(self, control_j, next_phys_msmt_j):
        '''docstring'''

        self.QubitGrid.nodes[control_j].physcmsmtsum = next_phys_msmt_j
        prob_j = self.sample_prob_from_msmts(control_j)

        for alpha_particle in self.AlphaSet.particles:
            alpha_particle.pset_beta = self.INITIALDICT["P_BETA"]
            alpha_particle.node_j = control_j
            map_index = self.QubitGrid.number_of_nodes*3 + alpha_particle.node_j 
            # TODO: TODO:1: alpha_particle.mean_radius_j = self.AlphaSet.posterior_state[map_index] 
            # TODO: TODO:1: neighbourbood mean radius update
        update_alpha_dictionary(next_phys_msmt_j, prob_j)

#       ------------------------------------------------------------------------
#       SUPPORT FUNCTION 3: PROPAGATE (ALPHA) STATES
#       ------------------------------------------------------------------------

    def PropagateState(self, control_j):
        '''docstring '''
        for alpha_particle in self.AlphaSet.particles:
            self.sample_from_transition_dist(alpha_particle, control_j)
 
    def sample_from_transition_dist(self, alpha_particle, control_j):
        '''docstring t > 1'''

        # sample_x = self.QubitGrid.get_all_nodes(["x_state"])
        # sample_y = self.QubitGrid.get_all_nodes(["y_state"])
        # #np.random.normal(loc=INITIALDICT["MU_W"],
        # #                            scale=INITIALDICT["SIG2_W"],
        # #                            size=self.QubitGrid.number_of_nodes*2)

        # sample_f = self.QubitGrid.get_all_nodes(["f_state"])
        # #sample_f += np.random.normal(loc=INITIALDICT["MU_F"],
        # #                            scale=INITIALDICT["SIG2_F"],
        # #                            size=self.QubitGrid.number_of_nodes)

        # sample_r = self.QubitGrid.get_all_nodes(["r_state"])
        # sample_r += np.random.normal(loc=INITIALDICT["MU_R"],
        #                             scale=INITIALDICT["SIG2_R"],
        #                             size=self.QubitGrid.number_of_nodes)

        # alpha_particle.particle = np.concatenate((sample_x, sample_y, sample_f, sample_r),
        #                                          axis=0)
 
#       ------------------------------------------------------------------------
#       SUPPORT FUNCTION 4: COMPUTE ALPHA WEIGHTS; GENERATE BETA WEIGHTS
#       ------------------------------------------------------------------------

    def ComputeWeights(self, control_j):
        '''docstring'''
 
        new_alpha_weights = self.AlphaSet.calc_weights_set() # Normlaised
        self.AlphaSet.weights_set = new_alpha_weights

        posterior_weights = []

        for alpha_particle in self.AlphaSet.particles:
            self.update_alpha_map_via_born_rule(control_j)
            beta_alpha_j_weights = self.generate_beta_layer(alpha_particle)

            # print "... ... Raw posterior beta alpha weights are", beta_alpha_j_weights
            posterior_weights.append(alpha_particle.weight*beta_alpha_j_weights)

        posterior_weights = np.asarray(posterior_weights).flatten()
        normalisation = np.sum(posterior_weights)
        normalised_posterior_weights = posterior_weights*(1.0/normalisation)

        return  normalised_posterior_weights


    def sample_radii(self, previous_length_scale):
        '''docstring'''

        if previous_length_scale < 0:
            print "Previous length scale is less than min:", previous_length_scale
            raise RuntimeError
        lower_bound = (previous_length_scale + self.R_min)*0.1 + self.R_min

        # if lower_bound < self.R_min:
        #     lower_bound = self.R_min
        # sample = np.random.uniform(low=self.R_max, high=self.R_max)
        # good_sample = False
        # while good_sample is False:
        #     sample = np.random.normal(loc=0,scale=self.dgrid*3.0) + previous_length_scale
        #     if sample > self.R_min:
        #         good_sample = True
        # sample = np.random.uniform(low=lower_bound, high=previous_length_scale+ delta)
        # sample = np.random.uniform(low=0.0, high=self.R_max)
        # sample = np.random.uniform(low=0.9, high=1.1)
        
        sample = np.random.uniform(low=lower_bound, high=self.R_max)
        # print "sample_radii gives a sample raidus = ", sample
        
        return sample #(self.R_min + self.R_max)*0.5


    def generate_beta_layer(self, alpha_particle):
        '''docstring'''

        len_idx = self.QubitGrid.number_of_nodes*3 + alpha_particle.node_j
        parent_alpha = alpha_particle.particle.copy()
        new_beta_state = parent_alpha * 1.0 # get property

        # posterior_r_at_j_prev_msmt = self.QubitGrid.get_all_nodes(["r_state"])[alpha_particle.node_j] # previous lengthscale doesnt matter right now

        list_of_parent_states = []
        list_of_length_samples = []
        for idx_beta in range(alpha_particle.pset_beta):
            # new_beta_state[len_idx] = self.sample_radii(parent_alpha[len_idx]) # shouldn't the argument be the posterior r at node j ?
            new_length_sample = self.sample_radii(0.0)
            new_beta_state[len_idx] = new_length_sample*1.0
            list_of_parent_states.append(new_beta_state.copy())
            list_of_length_samples.append(new_length_sample)

        # print "generate generate_beta_pset in generate_beta_layer from list of parents and radii "
        alpha_particle.generate_beta_pset(list_of_parent_states, list_of_length_samples)

        for beta_particle_object in alpha_particle.BetaAlphaSet_j.particles:
            # print "... ... ... for NEXT beta particle"
            # print "... ... ... ... Make a neighbourhood of qubits for each beta using new lengthscale sample..."
            self.generate_beta_neighbourhood(beta_particle_object) # compute smeared phases for each beta particle

        # print "... ... ... Calculate bset weights based on beta likelihood functiion..."
        beta_alpha_j_weights = alpha_particle.BetaAlphaSet_j.calc_weights_set()
        # print "... ... ... Check the shape of output ", beta_alpha_j_weights.shape
        # print "... ... ... Check normalisation beta_alpha_j_weights ", np.sum(beta_alpha_j_weights)
        # print "... ... ... RETURN beta_alpha_j_weights ", beta_alpha_j_weights.shape
        # print beta_alpha_j_weights
        # print
        # print "... ... ... EXIT generate_beta_layer"
        # print
        return beta_alpha_j_weights # these weights are normalised


    def generate_beta_neighbourhood(self, BetaParticle):
        '''docstring'''

        # BetaParticle.mean_radius = new_neighbourhood_L
        # print "... ... ... ... In generate_beta_neighbourhood..."
        # print "... ... ... ...  Pull information about the posterior state from t-1 QubitGrid..."
        NEIGHBOURDICT = {"prev_posterior_f_state" : self.QubitGrid.get_all_nodes(["f_state"]),
                         "prev_counter_tau_state" : self.QubitGrid.get_all_nodes(["counter_tau"]),
                         "lambda_" : self.INITIALDICT["LAMBDA"],
                         "kernel_function": self.INITIALDICT["kernel_function"]}

        # print 
        # print "... ... ... ...  The mean radius is... [before calling smear_fj_on_neighbours]", BetaParticle.mean_radius
        # print "... ... ... ...  List of neighbours are [before calling smear_fj_on_neighbours]:", BetaParticle.neighbourhood_qj
        # print "... ... ... ...  Smear phases in beta particle neighbood..."
        BetaParticle.smear_fj_on_neighbours(**NEIGHBOURDICT)
        # print
        # print "... ... ... ...  The mean radius is...", BetaParticle.mean_radius
        # print "... ... ... ...  List of neighbours are:", BetaParticle.neighbourhood_qj
        # print "... ... ... ...  List of smeared phases are:", BetaParticle.smeared_phases_qj
        # print "... ... ... ... EXIT generate_beta_neighbourhood"
        # print

#       ------------------------------------------------------------------------
#       SUPPORT FUNCTION 5: RESAMPLE AND BETA COLLAPSE
#       ------------------------------------------------------------------------

    def ResampleParticles(self, posterior_weights):
        '''docstring'''
        # print 
        # print "In ResampleParticles ..."
        resampled_idx = np.arange(self.pset_alpha*self.pset_beta)

        # print " ... Checking resamplign threshold using effective_particle_size ..."
        # if self.resample_thresh > self.effective_particle_size(posterior_weights):
        resampled_idx = ParticleFilter.resample_constant_pset_alpha(posterior_weights, self.pset_alpha, self.pset_beta)
        new_alpha_subtrees = self.get_subtrees(resampled_idx, self.pset_beta)
        new_alpha_list = self.collapse_beta(new_alpha_subtrees, resampled_idx)
        # print
        # self.pset_alpha = len(new_alpha_list)
        # print "...The new list of alphat particles should be same as previously: self.pset_alpha = %s, new list = %s" %(self.pset_alpha, len(new_alpha_list))

        # print "... Setting AlphaSet.particles  to new particles without a beta layer; and uniform weights"
        self.AlphaSet.particles = new_alpha_list # garanteed to be pset_alpha with no second layer
        self.AlphaSet.weights_set = (1.0/self.pset_alpha)*np.ones(self.pset_alpha)
        # print "post resampling and collapse, the alpha particles r_states are:", [alpha.particle[ self.QubitGrid.number_of_nodes*3 + alpha.node_j] for alpha in self.AlphaSet.particles]
        # print
        # print "EXIT ResampleParticles"
        # print

    def collapse_beta(self, subtree_list, resampled_indices):
        '''docstring'''
        # print 
        # print "... In collapse_beta, resampled indices are:", resampled_indices
        state_update = 0.
        new_alpha_particle_list = []
        for subtree in subtree_list:

            leaves_of_subtree = resampled_indices[subtree[0]:subtree[1]]
            leaf_count = float(len(leaves_of_subtree))
            # print "... ... The subtree is defined by the endpoint index boundaries", subtree
            # print "... ... The leaves of the subtree are ", leaves_of_subtree

            if leaf_count != 0:

                normaliser = (1./leaf_count)
                alpha_node = ParticleFilter.get_alpha_node_from_treeleaf(leaves_of_subtree[0], self.pset_beta)
                               # resampled_indices[subtree[0]], self.pset_beta)
                r_est_index = self.QubitGrid.number_of_nodes*3 + self.AlphaSet.particles[alpha_node].node_j
                beta_alpha_nodes = [ParticleFilter.get_beta_node_from_treeleaf(leafy, self.pset_beta) for leafy in leaves_of_subtree]
                               # resampled_indices[subtree[0]:subtree[1]]]
                # print "... ... The subtree has alpha node of: ", alpha_node
                # print "... ... The leaves of the subtree are labeled by beta indices", beta_alpha_nodes
                r_est_subtree_list = []

                for node in beta_alpha_nodes:
                
                    beta_state = self.AlphaSet.particles[alpha_node].BetaAlphaSet_j.particles[node].particle.copy()
                    beta_lengthscale = beta_state[r_est_index]*1.0
                    if np.isnan(beta_lengthscale): 
                        # print "A resampled beta_lengthscale has an invalid value"
                        raise RuntimeError
                    # r_est_subtree += normaliser*beta_lengthscale # collapses beta layer by taking the mean
                    r_est_subtree_list.append(beta_lengthscale)
                
                # print "in collapse beta, r_est_subtree_list = ", r_est_subtree_list
                r_est_subtree = ParticleFilter.calc_posterior_lengthscale(np.asarray(r_est_subtree_list)) # new posterior for alpha based on mode
                parent = self.AlphaSet.particles[alpha_node].particle.copy()*1.0
                parent[r_est_index] = r_est_subtree

                if np.any(np.isnan(parent)):
                    # print "A resampled parent particle has an invalid value"
                    raise RuntimeError

                # Beta Layer Collapsed
                # print "... ... Collapsing the beta layers for new alpha particles"
                # print "Mode length", r_est_subtree
                # print "Posterior parent = ", parent
                self.AlphaSet.particles[alpha_node].particle = parent*1.0
                self.AlphaSet.particles[alpha_node].BetaAlphaSet_j = None

                # New Alphas Stored
                new_alpha_particle_list.append(self.AlphaSet.particles[alpha_node])

        # print "... EXIT collapse_beta"
        # print
        return new_alpha_particle_list


    def effective_particle_size(self, posterior_weights):
        '''docstring'''

        p_size = 1.0/ np.sum(posterior_weights**2)
        # lower_b = (1.0 / float(self.pset_alpha + self.pset_beta))
        # if p_size < lower_b:
        #     p_size =  lower_b
        # if p_size > (self.pset_alpha + self.pset_beta):
        #     p_size = 1
        # self.L_factor = p_size*self.dgrid
        # print "... ... The p_size of these posterior weights is ", p_size
        # print "... ... The self.L_factor of these posterior weights is ", self.L_factor
        # print 
        # print "... EXIT effective_particle_size "
        # print
        return p_size


#       ------------------------------------------------------------------------
#       SUPPORT FUNCTION 6: STATIC METHODS
#       ------------------------------------------------------------------------
    @staticmethod
    def calc_posterior_lengthscale(r_lengthscales_array):
        '''docstring'''
        r_posterior_post_beta_collapse = ParticleFilter.calc_skew(r_lengthscales_array)
        return r_posterior_post_beta_collapse

    @staticmethod
    def calc_skew(r_lengthscales_array):
        '''docstring'''
        totalcounts = len(r_lengthscales_array)
        mean_ = np.mean(r_lengthscales_array)
        mode_, counts = mode(r_lengthscales_array)
        median_ = np.sort(r_lengthscales_array)[int(totalcounts/2) - 1] # not sure if this is median

        if mean_ < mode_ and counts > 1:
            # print "calc_skew gives skew left r, returning mode", mode_
            return mode_
        if mean_ > mode_ and counts > 1:
            # print "calc_skew gives skew left r, returning mode", mode_
            return mode_
        # print "calc_skew gives returns numerical mean", mean_
        return mean_


    @staticmethod
    def compute_dist(one_pair):
        '''doctring'''
        xval, yval = one_pair
        return np.sqrt((xval[0] - yval[0])**2 + (xval[1] - yval[1])**2)

    @staticmethod
    def get_distances(list_of_positions):
        '''docstring'''
        distance_pairs = [a_pair for a_pair in combinations(list_of_positions, 2)]
        distances = [ParticleFilter.compute_dist(one_pair)for one_pair in distance_pairs]
        return distances, distance_pairs

    @staticmethod
    def find_max_distance(list_of_positions):
        '''docstring'''
        distances, distance_pairs = ParticleFilter.get_distances(list_of_positions)
        return max(distances), distance_pairs[np.argmax(distances)]

    @staticmethod
    def find_min_distance(list_of_positions):
        '''docstring'''
        distances, distance_pairs = ParticleFilter.get_distances(list_of_positions)
        return min(distances), distance_pairs[np.argmin(distances)]

    @staticmethod
    def set_uniform_prior_for_correlation(list_of_positions, multiple=10):
        d_grid,_ = ParticleFilter.find_max_distance(list_of_positions)
        R_max = multiple*d_grid
        R_min,_  = ParticleFilter.find_min_distance(list_of_positions)
        return d_grid, R_min, R_max # d_grid, R_min, R_max

    @staticmethod
    def get_alpha_node_from_treeleaf(leaf_index, pset_beta):
        '''docstring'''
        alpha_node = int(leaf_index//float(pset_beta))
        return alpha_node

    @staticmethod
    def get_beta_node_from_treeleaf(leaf_index, pset_beta):
        '''docstring'''
        beta_node = int(leaf_index - int(leaf_index//float(pset_beta))*pset_beta)
        return beta_node

    @staticmethod
    def resample_from_weights(posterior_weights, number_of_samples):
        '''docstring'''
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
                    # # print('Break - max particle index reached during sampling')
                    break   # clip at max particle index, plus zero
            resampled_idx.append(j-1) # sgift down to match python indices

        return resampled_idx

    @staticmethod
    def resample_constant_pset_alpha(posterior_weights, pset_alpha, pset_beta):
        '''Returns indicies for particles picked after sampling from posterior'''
        # DO WEIGHTS NEED TO BE SORTED? (No, the CDF will be
        # monotnically increasing; and the shape of the CDF will determine
        # the frequency at which (x,y) are sampled if y is uniform )

        sufficient_sample = False
        num_of_samples = pset_alpha
        total_particles = len(posterior_weights)

        if total_particles != int(INITIALDICT["P_ALPHA"]*INITIALDICT["P_BETA"]):
            # print "In resample_constant_pset_alpha, total weights != P_alpha * P_beta"
            raise RuntimeError

        while sufficient_sample is False:

            num_of_samples += 5
            # print "In resample_constant_pset_alpha, number of sufficient samples are:", num_of_samples
            resampled_indices = ParticleFilter.resample_from_weights(posterior_weights, num_of_samples)
            resampled_alphas = [ParticleFilter.get_alpha_node_from_treeleaf(leafy, pset_beta) for leafy in resampled_indices]
            unique_alphas = set(list(resampled_alphas))
            if len(unique_alphas) == pset_alpha:
                sufficient_sample = True

        return resampled_indices

    @staticmethod
    def get_subtrees(resampled_indices, pset_beta):
        '''docstring'''

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
            # end_counter += 1
            new_sub_trees.append([strt_counter, end_counter])
        return new_sub_trees
    # def collapse_beta(self, new_sub_trees, resampled_indices):
    #     '''docstring'''

    #     state_update = 0.
    #     pset_beta = self.pset_beta

    #     new_alpha_particle_list = []
    #     for pairs in new_sub_trees:

    #         subtree = resampled_indices[pairs[0]:pairs[1]]
    #         leaf_count = float(len(subtree))

    #         # print
    #         print "In collapse_beta, and checking variables"
    #         print "subtree", subtree
    #         print "leaf_count", leaf_count


    #         if leaf_count != 0:

    #             normaliser = (1./leaf_count)
    #             alpha_node = ParticleFilter.get_alpha_node_from_treeleaf(pairs[0], pset_beta)
    #             beta_alpha_nodes = [ParticleFilter.get_beta_node_from_treeleaf(leafy, pset_beta) for leafy in subtree]
