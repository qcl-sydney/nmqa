import qslamr as qs
import numpy as np

#### THIS MODULE IS THE SAME AS SINGLERUN_V2 EXCEPT THAT SINGLERUN_V2
# SAVES 2 ADDITIONAL DATA OUTPUTS (POSTERIOR AND MARGINALISED STATES)
# AND REPETITIONS ARE EXCLUDED FROM SINGLERUN_V2, WHILE THEY ARE ENDOGENOUS HERE.

class SingleRunAnalysis(object):

    def __init__(self, 
                 SAMPLE_GLOBAL_MODEL,true_map_,repts,
                 measurements_controls_= None,
                 autocontrol_="ON",
                 var_thres_=1.0,
                 numofnodes =25,
                 beta_expansion_mode=False,
                 beta_skew_adjust=False):

        self.SAMPLE_GLOBAL_MODEL = SAMPLE_GLOBAL_MODEL
        self.true_map_ = true_map_
        self.repts = repts
        self.measurements_controls_=measurements_controls_
        self.autocontrol_=autocontrol_
        self.var_thres_=var_thres_
        self.numofnodes=len(SAMPLE_GLOBAL_MODEL["GRIDDICT"])
        self.beta_expansion_mode=beta_expansion_mode
        self.beta_skew_adjust=beta_skew_adjust
        
    def call_qslam(self):

        qslamobj=0
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

    def run_analysis(self, filename='./SingleRunAnalysis'):

        max_num_iterations=self.SAMPLE_GLOBAL_MODEL["MODELDESIGN"]["MAX_NUM_ITERATIONS"]

        # --------------------------------------------------------------------
        # Make true error matrix and rate of change matrix for a single run
        # --------------------------------------------------------------------

        absolute_errors_matrix = np.zeros((self.repts, max_num_iterations, self.numofnodes ))
        rate_of_change_of_map = np.zeros((self.repts, max_num_iterations, self.numofnodes ))
        rate_of_change_of_len = np.zeros((self.repts, max_num_iterations, self.numofnodes ))

        for idx_rept in range(self.repts):

            qslamobj = self.call_qslam()

            for idxt in range(0, max_num_iterations): 

                f_strt = self.numofnodes * 2  # fixed by QSLAM design of state vector
                f_end = self.numofnodes * 3
                r_strt = self.numofnodes * 3

                absolute_errors_matrix[idx_rept, idxt, :]  = abs(qslamobj.save_posterior_state[idxt][f_strt:f_end] - self.true_map_)

                if idxt > 0:

                    rate_of_change_of_map[idx_rept, idxt, :] = abs(qslamobj.save_posterior_state[idxt][f_strt:f_end] - qslamobj.save_posterior_state[idxt-1][f_strt:f_end])
                    rate_of_change_of_len[idx_rept, idxt, :] = abs(qslamobj.save_posterior_state[idxt][r_strt:] - qslamobj.save_posterior_state[idxt-1][r_strt:])

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

            posterior_weights[idxt, :] = np.asarray([qslamobj.save_alpha_posterior[idxt][idx].weight for idx in range(P_ALPHA)])

        np.savez(filename,
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
                 posterior_weights=posterior_weights)

        return

