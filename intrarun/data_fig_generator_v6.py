'''
This python script stores a dictionary object to call, analyze and plot final simulation results, as stored in './data_v3'.
The output of this file is to generate images and store time in folder structure './truth_flag_*/

Truth flag type (five options)
# 0 - Linear 1D field, 25 qubit array
# 1 - Square 2D field, 25 qubit array
# 2 - Gaussian 2D array, 25 qubit array
# 3 - Square 2D field, 9 qubit array
# 4 - Square 2D field, 16 qubit array

'''

import numpy as np
import matplotlib
import sys
import os
sys.path.append('../qslam/')
from visualiserisk import *

###########################################################
# USER TO SPECIFY
###########################################################
# User to specify truth flag
idx_prefix = int(sys.argv[1])
path = './data_v3/'
prefix = '_idx_prefix_'+str(idx_prefix)+'_'

###########################################################
# DATA SOURCES AND SAVE PATHS
###########################################################
likelihoodparams_scan = np.load('./random_variances.npz')
lambda_pairs = np.load('./lambda_pairs_2.npz')
# savetopath = './truth_flag_'+str(idx_prefix)+'/'
savetopath = './data_figures/'

###########################################################
# SIMULATION SETTINGS (DO NOT CHANGE)
###########################################################
datatype=['Uni_R', 'Trunc_R']
keys =  ['marginalised_weights', 'predictive_weights', 'alpha_labels', 'beta_labels', 'rate_of_change_of_len', 'joint_labels', 'posterior_weights', 'joint_weights', 'leaf_weights']
particleconfigs = [ [3,2], [9,6], [15,10] , [21,14], [30, 20]]
paramlist =['optimal', 'zerolambda']

SIMULATIONSDICT = {}

for idx_truth_flag in range(5):
    
    SIMULATIONSDICT[idx_truth_flag] = {}
    
    SIMULATIONSDICT[idx_truth_flag]["num_of_nodes"] = 25
    SIMULATIONSDICT[idx_truth_flag]["linear"] = False
    SIMULATIONSDICT[idx_truth_flag]["repts"] = 50
    SIMULATIONSDICT[idx_truth_flag]["max_iterations"] = 75
    
    for idx_expandtype in ["Uniform", "TruncGauss"]:
    
        SIMULATIONSDICT[idx_truth_flag][idx_expandtype] = {}
        SIMULATIONSDICT[idx_truth_flag][idx_expandtype]["optimal"] = {}
        SIMULATIONSDICT[idx_truth_flag][idx_expandtype]["optimal"]["idx_1"] = None
        SIMULATIONSDICT[idx_truth_flag][idx_expandtype]["optimal"]["idx_2"] = None   
    
        SIMULATIONSDICT[idx_truth_flag][idx_expandtype]["zerolambda"] = {}
        SIMULATIONSDICT[idx_truth_flag][idx_expandtype]["zerolambda"]["idx_1"] = None
        SIMULATIONSDICT[idx_truth_flag][idx_expandtype]["zerolambda"]["idx_2"] = 0 

    SIMULATIONSDICT[idx_truth_flag]["GRAPHMIN"] = -2.5# -1.9 # Figure plotting
    SIMULATIONSDICT[idx_truth_flag]["GRAPHMAX"] = 0.7 #1 # Figure plotting
    SIMULATIONSDICT[idx_truth_flag]["R_MAX"] = 150 # Figure plotting
    SIMULATIONSDICT[idx_truth_flag]["F_MAX"] = 8.5 # Figure plotting
    

# Specific settings
SIMULATIONSDICT[0]["linear"] = True
SIMULATIONSDICT[3]["num_of_nodes"] = 9
SIMULATIONSDICT[4]["num_of_nodes"] = 16

SIMULATIONSDICT[0]["R_MAX"] = 300 # Figure plotting
SIMULATIONSDICT[1]["R_MAX"] = 80 # Figure plotting
SIMULATIONSDICT[2]["R_MAX"] = 80 # Figure plotting
SIMULATIONSDICT[3]["R_MAX"] = 20 # Figure plotting
SIMULATIONSDICT[4]["R_MAX"] = 75 # Figure plotting

###########################################################
# NUMERICAL OPTIMISATION RESULTS
###########################################################

# trunc
SIMULATIONSDICT[0]["TruncGauss"]["optimal"]["idx_1"] = 5
SIMULATIONSDICT[0]["TruncGauss"]["optimal"]["idx_2"] = 18 
SIMULATIONSDICT[0]["TruncGauss"]["zerolambda"]["idx_1"] = 5
SIMULATIONSDICT[0]["TruncGauss"]["zerolambda"]["idx_2"] = 0 

SIMULATIONSDICT[0]["Uniform"]["optimal"]["idx_1"] = 29
SIMULATIONSDICT[0]["Uniform"]["optimal"]["idx_2"] = 18 
SIMULATIONSDICT[0]["Uniform"]["zerolambda"]["idx_1"] = 29
SIMULATIONSDICT[0]["Uniform"]["zerolambda"]["idx_2"] = 0 

# trunc
SIMULATIONSDICT[1]["TruncGauss"]["optimal"]["idx_1"] = 16
SIMULATIONSDICT[1]["TruncGauss"]["optimal"]["idx_2"] = 18 
SIMULATIONSDICT[1]["TruncGauss"]["zerolambda"]["idx_1"] = 16
SIMULATIONSDICT[1]["TruncGauss"]["zerolambda"]["idx_2"] = 0 

SIMULATIONSDICT[1]["Uniform"]["optimal"]["idx_1"] = 10
SIMULATIONSDICT[1]["Uniform"]["optimal"]["idx_2"] = 18 
SIMULATIONSDICT[1]["Uniform"]["zerolambda"]["idx_1"] = 10
SIMULATIONSDICT[1]["Uniform"]["zerolambda"]["idx_2"] = 0 

# uniform
SIMULATIONSDICT[2]["Uniform"]["optimal"]["idx_1"] = 29 # Min Trunc: 21
SIMULATIONSDICT[2]["Uniform"]["optimal"]["idx_2"] = 13 # Min Trunc: 13
SIMULATIONSDICT[2]["Uniform"]["zerolambda"]["idx_1"] = 29 # Min Trunc: 21
SIMULATIONSDICT[2]["Uniform"]["zerolambda"]["idx_2"] = 0 

SIMULATIONSDICT[2]["TruncGauss"]["optimal"]["idx_1"] = 21 
SIMULATIONSDICT[2]["TruncGauss"]["optimal"]["idx_2"] = 13 
SIMULATIONSDICT[2]["TruncGauss"]["zerolambda"]["idx_1"] = 21 
SIMULATIONSDICT[2]["TruncGauss"]["zerolambda"]["idx_2"] = 0 

# trunc
SIMULATIONSDICT[3]["TruncGauss"]["optimal"]["idx_1"] = 4
SIMULATIONSDICT[3]["TruncGauss"]["optimal"]["idx_2"] = 8
SIMULATIONSDICT[3]["TruncGauss"]["zerolambda"]["idx_1"] = 4
SIMULATIONSDICT[3]["TruncGauss"]["zerolambda"]["idx_2"] = 0 

SIMULATIONSDICT[3]["Uniform"]["optimal"]["idx_1"] = 10
SIMULATIONSDICT[3]["Uniform"]["optimal"]["idx_2"] = 12
SIMULATIONSDICT[3]["Uniform"]["zerolambda"]["idx_1"] = 10
SIMULATIONSDICT[3]["Uniform"]["zerolambda"]["idx_2"] = 0 

# uniform
SIMULATIONSDICT[4]["Uniform"]["optimal"]["idx_1"] = 18 # Min Trunc: 18
SIMULATIONSDICT[4]["Uniform"]["optimal"]["idx_2"] = 18 # Min Trunc: 12
SIMULATIONSDICT[4]["Uniform"]["zerolambda"]["idx_1"] = 18 # Min Trunc: 18
SIMULATIONSDICT[4]["Uniform"]["zerolambda"]["idx_2"] = 0 

SIMULATIONSDICT[4]["TruncGauss"]["optimal"]["idx_1"] = 18 
SIMULATIONSDICT[4]["TruncGauss"]["optimal"]["idx_2"] = 12 
SIMULATIONSDICT[4]["TruncGauss"]["zerolambda"]["idx_1"] = 18 
SIMULATIONSDICT[4]["TruncGauss"]["zerolambda"]["idx_2"] = 0 


###########################################################
# FIGURE GENERATION
###########################################################
max_iterations = SIMULATIONSDICT[idx_prefix]["max_iterations"]
num_of_nodes = SIMULATIONSDICT[idx_prefix]["num_of_nodes"]
GRAPHMIN = SIMULATIONSDICT[idx_prefix]["GRAPHMIN"]
GRAPHMAX = SIMULATIONSDICT[idx_prefix]["GRAPHMAX"]
R_MAX = SIMULATIONSDICT[idx_prefix]["R_MAX"]
F_MAX = SIMULATIONSDICT[idx_prefix]["F_MAX"]

for idx_figure in range(2):

    # Plots optimal case (idx_figure=0), and then zero lambda case (idx_figure=1)
    
    idx_paramtype = paramlist[idx_figure]
    idx_1_U = SIMULATIONSDICT[idx_prefix]["Uniform"][idx_paramtype]["idx_1"]
    idx_2_U = SIMULATIONSDICT[idx_prefix]["Uniform"][idx_paramtype]["idx_2"]
    idx_1 = SIMULATIONSDICT[idx_prefix]["TruncGauss"][idx_paramtype]["idx_1"]
    idx_2 = SIMULATIONSDICT[idx_prefix]["TruncGauss"][idx_paramtype]["idx_2"]
    
    BASELINE = False
    if idx_prefix < 3  and idx_figure==0 :
        BASELINE = True

    #################################
    # DATA ANALYSIS
    #################################
    
    error_scaling_matrix_0 = 100.*np.ones((len(particleconfigs), max_iterations))
    error_scaling_matrix_1 = 100.*np.ones((len(particleconfigs), max_iterations))


    for idx_pconfig in range(len(particleconfigs)):

        test_data_0 = np.load(path+datatype[0]+prefix+'rand_'+str(idx_1_U)+'_'+str(idx_2_U)+'_'+str(idx_pconfig)+'.npz')
        error_scaling_matrix_0[idx_pconfig, :] = np.mean(np.mean(test_data_0['absolute_errors_matrix']**2, axis=0), axis=1)

        test_data_1 = np.load(path+datatype[1]+prefix+'rand_'+str(idx_1)+'_'+str(idx_2)+'_'+str(idx_pconfig)+'.npz')
        error_scaling_matrix_1[idx_pconfig, :] = np.mean(np.mean(test_data_1['absolute_errors_matrix']**2, axis=0), axis=1)

    particle_number = np.asarray([idp[0] for idp in particleconfigs] ) 
    
    pick_t_regime = 74
    idx_pconfig_T = np.argmin(error_scaling_matrix_1[:, pick_t_regime])
    idx_pconfig_U = np.argmin(error_scaling_matrix_0[:, pick_t_regime])
    
    if BASELINE:
        
        # Only Uniform expansion used in baseline:
        # baseline_0['absolute_errors_matrix'].shape is (50, SIMULATIONSDICT[idx_truth_flag]["max_iterations"], 25) (repts, iter, grid)
        baseline_error_scaling_0 = 100.*np.ones((1, max_iterations))
        baseline_0 = np.load(path+datatype[0]+prefix+'baseline_.npz')
        baseline_error_scaling_0[0, :] = np.mean(np.mean(baseline_0['absolute_errors_matrix']**2, axis=0), axis=1)
        
    #################################
    # FIGURE : expansion behaviour
    #################################

    gslayout = gs(10,1, top = 0.99, bottom =0.05, left = 0.15, right = 0.98, wspace = 0.05, hspace = 0.4) # MAIN TEXT
    # gslayout = gs(1, 1, top = 0.99, bottom =0.1, left = 0.1, right = 0.98, wspace = 1.0, hspace = 0.01)  # APPENDIX

    fig = plt.figure(figsize=(cm2inch(6.),cm2inch(8))) # MAIN TEXT
    # fig = plt.figure(figsize=(cm2inch(7.8),cm2inch(4.5))) # APPENDIX 

    plt.suptitle('TRUNC (g1, g2) = (%s, %s) and (l1, l2) = (%s, %s) UNI (g1, g2) = (%s, %s) and (l1, l2) = (%s, %s)' %(likelihoodparams_scan['g1var'][idx_1], 
                                                                 likelihoodparams_scan['g2var'][idx_1],
                                                                 lambda_pairs['lambda_1'][idx_2],
                                                                 lambda_pairs['lambda_2'][idx_2],
                                                                 likelihoodparams_scan['g1var'][idx_1_U], 
                                                                 likelihoodparams_scan['g2var'][idx_1_U],
                                                                 lambda_pairs['lambda_1'][idx_2_U],
                                                                 lambda_pairs['lambda_2'][idx_2_U]),
                )

    ax_err = fig.add_subplot(gslayout[0:7, 0])  # MAIN TEXT
    # ax_err = fig.add_subplot(gslayout[0, 0])  # APPENDIX
    ax1 = fig.add_subplot(gslayout[7:, 0]) # MAIN TEXT
    
    #################################
    # top: error trajectory
    #################################
    
    if BASELINE:
        ax_err.plot(range(SIMULATIONSDICT[idx_truth_flag]["max_iterations"])[::2],np.log(baseline_error_scaling_0[0, ::2]), '-', 
                    lw=1., ms=4, c='k', 
                    label='Baseline')
        

    ax_err.plot(range(SIMULATIONSDICT[idx_truth_flag]["max_iterations"])[::2],np.log(error_scaling_matrix_0[idx_pconfig_U, ::2]), 'x--', 
                lw=1., ms=4, c='indianred', 
                label='Uniform')



    ax_err.plot(range(SIMULATIONSDICT[idx_truth_flag]["max_iterations"])[::2],np.log(error_scaling_matrix_1[idx_pconfig_T, ::2]), 'o--', 
                lw=1., ms=4, 
                c='steelblue',
                markerfacecolor='None', label='Trunc. Gauss', alpha=1.)

    ax_err.set_xlabel(r'$t$  (num)')
    ax_err.set_ylabel(r'$\log(\mathbb{E}[|\pi^n_t - \pi_\infty|])$  (a.u.)')
    ax_err.set_ylim([GRAPHMIN, GRAPHMAX])
    ax_err.axvline(x=num_of_nodes, ls='-', color='grey', linewidth=0.5, label='Num. of qubits')
    ax_err.legend(loc=0)
    ax_err.xaxis.set_ticklabels([]) # MAIN TEXT
    
    #################################
    # bottom: Rate of change of R
    #################################

    test_data_0 = np.load(path+datatype[0]+prefix+'rand_'+str(idx_1_U)+'_'+str(idx_2_U)+'_'+str(idx_pconfig_T)+'.npz')
    test_data_1 = np.load(path+datatype[1]+prefix+'rand_'+str(idx_1)+'_'+str(idx_2)+'_'+str(idx_pconfig_T)+'.npz')
    rate_of_change_of_len_0 = test_data_0['rate_of_change_of_len'][np.random.randint(low=0, high=50)]
    rate_of_change_of_len_1 = test_data_1['rate_of_change_of_len'][np.random.randint(low=0, high=50)]
    
    # MAIN TEXT
    ax1.plot(np.sum(rate_of_change_of_len_0, axis=1)[1:], ls='--', c='indianred', lw=1.1, label='Uniform')
    ax1.plot(np.sum(rate_of_change_of_len_1, axis=1)[1:], ls='-', c='steelblue', lw=1.2, label='Trunc. Gauss')
    ax1.set_ylim([-0.1, R_MAX])
    ax1.legend(loc=0)
    
    plt.savefig(savetopath+'fig_'+str(idx_paramtype)+'_expansion_v6_'+str(idx_prefix)+'.svg', format='svg', dpi=800)
    plt.show()
    #################################
    # FIGURE: Scaling Behaviour
    #################################
    
    fig = plt.figure(figsize=(cm2inch(7.8),cm2inch(4.5))) # MAIN TEXT # APPENDIX
    # fig = plt.figure(figsize=(cm2inch(3.5),cm2inch(3.))) # ZEROLAMBDA
    
    gslayout = gs(1, 1, top = 0.99, bottom =0.1, left = 0.1, right = 0.98, wspace = 1.0, hspace = 0.01)
    
    #################################
    # left: Delta vs. t behaviour
    #################################
    ax_scl = fig.add_subplot(gslayout[0, 0])
    
    delta_trace_0 = []
    delta_trace_1 = []
    
    for idx_t in range(SIMULATIONSDICT[idx_truth_flag]["max_iterations"]):

        slope_0, intercept_0 = np.polyfit(np.log(particle_number), np.log(error_scaling_matrix_0[:, idx_t]), 1)
        slope_1, intercept_1 = np.polyfit(np.log(particle_number), np.log(error_scaling_matrix_1[:, idx_t]), 1)
        delta_trace_0.append(slope_0)
        delta_trace_1.append(slope_1)
    
    delta_trace_0 = np.asarray(delta_trace_0)
    delta_trace_1 = np.asarray(delta_trace_1)
    ax_scl.plot(range(SIMULATIONSDICT[idx_truth_flag]["max_iterations"])[::2], delta_trace_0[::2], 'x--', 
                lw=1., ms=2, c='indianred', 
                label='Uniform')

    ax_scl.plot(range(SIMULATIONSDICT[idx_truth_flag]["max_iterations"])[::2], delta_trace_1[::2], 'o--', 
                lw=1., ms=2, 
                c='steelblue',
                markerfacecolor='None', label='Trunc. Gauss', alpha=1.)
    
    ax_scl.set_ylim([-0.3, 0.75]) # MAIN TEXT
    ax_scl.set_yticks(np.arange(-0.3, 0.76, 0.15)) # MAIN TEXT
    # ax_scl.set_ylim([-0.3, 0.8]) # APPENDIX 
    # ax_scl.set_yticks(np.arange(-0.3, 0.81, 0.15))   # APPENDIX
    # ax_scl.set_ylim([-0.3, 0.45]) # ZEROLAMBDA 
    # ax_scl.set_yticks(np.arange(-0.3, 0.451, 0.15))   # ZEROLAMBDA
    ax_scl.set_xlim([0, 75])
    ax_scl.axhline(y=0.0, ls='-', color='k', lw=0.5)
    
    
    plt.savefig(savetopath+'fig_'+str(idx_paramtype)+'_scalingbehav_v6_'+str(idx_prefix)+'.svg', format='svg', dpi=800)
    plt.show()
    
    #################################
    # inset: Error scaling vs particle number and best fit at t=75
    #################################
    
    fig = plt.figure(figsize=(cm2inch(2.75),cm2inch(1.9)))
    gslayout = gs(1, 1, top = 0.99, bottom =0.1, left = 0.1, right = 0.98, wspace = 1.0, hspace = 0.01)
    ax_inset = fig.add_subplot(gslayout[0, 0])
    
    
    # Uniform
    slope_0, intercept_0 = np.polyfit(np.log(particle_number), np.log(error_scaling_matrix_0[:, pick_t_regime]), 1)
    ax_inset.plot(np.log(particle_number), np.log( error_scaling_matrix_0[:, pick_t_regime]), 
                'x', ms=3.2, markerfacecolor='None', lw=1., c='indianred', label= 'Uniform')
    trendpoly_0 = np.poly1d([slope_0,intercept_0]) 
    ax_inset.plot(np.log(particle_number),trendpoly_0(np.log(particle_number)), '--',  lw=1., c='indianred', label= r'$\Delta = %s$'%(np.round(slope_0, 2)))
    ax_inset.axvline(x=np.log(particleconfigs[idx_pconfig_U][0]), ls='-', lw=1., c='indianred', alpha=0.5)
    
    # Trunc Gaussian
    slope_1, intercept_1 = np.polyfit(np.log(particle_number), np.log(error_scaling_matrix_1[:, pick_t_regime]), 1)
    ax_inset.plot(np.log(particle_number), np.log( error_scaling_matrix_1[:, pick_t_regime]), 
                'o', ms=3.2, markerfacecolor='None', lw=1., c='steelblue', label= 'Trunc. Gauss')
    trendpoly_1 = np.poly1d([slope_1,intercept_1]) 
    ax_inset.plot(np.log(particle_number),trendpoly_1(np.log(particle_number)), '--',  lw=1., c='steelblue', label= r'$\Delta = %s$'%(np.round(slope_1, 2)))
    ax_inset.axvline(x=np.log(particleconfigs[idx_pconfig_T][0]), ls='-', lw=1., c='steelblue', alpha=0.5)
    
    # Both
    ax_inset.set_ylim([-1.9, -1.2]) # Zoom into the t=75 case. # MAIN TEXT
    ax_inset.set_yticks(np.arange(-1.9, -1.21, 0.3))  # MAIN TEXT
    # ax_inset.set_ylim([-2.5, -1.]) # Zoom into the t=75 case. # APPENDIX #ZEROLAMBDA
    # ax_inset.set_yticks(np.arange(-2.5, -1., 0.5)) # APPENDIX #ZEROLAMBDA
    ax_inset.legend(loc=1, fontsize=8)
    
    # ax_inset.yaxis.set_ticklabels([])
    
    plt.savefig(savetopath+'fig_'+str(idx_paramtype)+'_scalinginset_v6_'+str(idx_prefix)+'.svg', format='svg', dpi=800)
    plt.show()
    


