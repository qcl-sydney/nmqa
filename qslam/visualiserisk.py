'''
Created on Thu Apr 20 19:20:43 2017
@author: riddhisw

.. module:: visualiserisk

    :synopsis: Container for plotting preamble, functions and custom data handling
        classes  to visualise Bayes Risk for QSLAM and Naive Brute Force measurement
        comparisons.

    Module Level Classes:
    ---------------------
        Metric: Collection of static methods to define a cost metric
            for Bayes Risk analysis.
        DataCube: Class object to handle updates to model parameters for QSLAM
            and EngineeredTruth during simulations.
        qPlotter: Class object for aiding visualisation of Bayes analysis.

    Module Level Dictionaries:
    -------------------------

    Data handling dictionaries

        DATAVIEWKEYS : Selects QSLAM data to plot based on output of
            Bayes Risk analysis:
                "truth" : Set of true maps from Bayes risk trials.
                "pred_f" : Set of posterior f_states from Bayes risk trials.
                "pred_r" : Set of posterior r_states from Bayes risk trials.
                "path" : Set of control paths from Bayes risk trials.
                "errs" : Set of residual errors from Bayes risk trials.

        FLAG : Distinguishes between QSLAM ('q') and naive ('n')
        NPZFLAG : Distinguishes between QSLAM ('q') and naive ('n')
            data file name suffix.
        PATHDICT : Specifies path to file and file name for a datafile
            "pdir": path to directory for outputs from Bayes Risk analysis
            "fle": file name for a specific .npz output from Bayes Risk analysis

    Plotting dictionaries

        PKWG : Dictionary object for control path plotting parameters
            "pc":       path base color
            "palpha"    path color alpha value
            "pss":      path start / stop size
            "ac":       arrow base color
            "aalpha":   arrow color alpha value
            "adxdy":    arrow dx and dy as ratio of path segment
            "alw":      arrow line width
            "ahw":      arrow head width

        HEATMAP : Dictionary object for controlling heatmap parameters.
            Keys: "vmin", "vmax", "cmap", "origin" for matplotlib
        CORRMAP :  Dictionary object for controlling correlation map parameters.
            Keys: "vmin", "vmax", "cmap", "origin" for matplotlib

    Module Level Functions:
    ----------------------

    path_to_file : Return path to file
    check_savefig: Return True if parameter regimes for plotting match controls'
    get_control_path : Returns qubit positions (x,y) visited based on a control path
         and input dictionary of qubits
    cm2inch :  Return values in inches for input values in cm

.. moduleauthor:: Riddhi Gupta <riddhi.sw@gmail.com>
'''
import numpy as np
import matplotlib
import sys
import os

################################################################################
# COST METRIC DEFINITIONS
################################################################################

class Metric(object):

    ''' Collection of static methods to define a cost metric
        for Bayes Risk analysis. Cost metrics are a function of
        residual error between a true map and an algorithmic
        estimate of f_state.

    Static Methods:
    ---------------
        infnorm :
            Return the inf-norm error metric over many trials via e_type.
        
        ssim :
            Return the expected value of SSIM for a set of true maps and map estimates;
            as well as the list of SSIM scores for a set of true maps and map estimates over
            which the expectation was calculated.
        score_ssim : Return SSIM score as a deviation from ideal. Helper function
            to ssim().

        rms :
            Return the expected value of root mean square error metric
            for a set of true maps and map estimates,  normalised to maximal error
            of pi radians per pixel and total number of grid points for each map.

        singlemap_rmse : Return the root mean square error metric for
            a true map and a map estimate normalised to maximal error of pi radians.
            Helper function for rms().

        singlemap_err : Identical to  singlemap_rmse [NOT USED]
        original_err_metric [NOT USED]
        err [NOT USED]
    '''

    def __init__(self):
        '''
        Creates a Metric instance.
        '''
        pass
    
    
    @staticmethod
    def infnorm(macro_residuals, e_type="maxinf"):
        '''
        Return the inf-norm error metric over many trials via e_type.
        
        '''
        residuals =  macro_residuals[0, :, :]
        dims = residuals.shape
        error = np.zeros(dims[0])
        
        # ignore np.nan values used as "fillers" in scipy interpolation functions
        # np.mean --> np.nanmean
        # np.max --> np.nanmax
        
        for idx_trial in range(dims[0]):
            error[idx_trial] = np.nanmax(abs(residuals[idx_trial, :]))
        
        if e_type =='maxinf':
            return np.max(error)
        
        if e_type =='expinf':
            return np.nanmean(error) 
            
        print("Invalid e_type")
        raise RuntimeError
        
        
    @staticmethod
    def rms(macro_residuals):
        ''' Return the expected value of root mean square error metric
            for a set of true maps and map estimates,  normalised to maximal error
            of pi radians per pixel and total number of grid points for each map.
        '''
        rmse_vector = Metric.singlemap_rmse(macro_residuals[0, :, :], axis=0)
        normaliser = 1.0 / np.sqrt(rmse_vector.shape[0])
        error = normaliser*np.linalg.norm(rmse_vector)
        return error

    @staticmethod
    def singlemap_rmse(residualsdata, axis):
        ''' Return the root mean square error metric for
            a true map and a map estimate normalised to maximal error of pi radians.
            Helper function for rms().

        Normalisation: Since a single qubit projective measurement can maximally only resolve
        a np.pi relative phase between 0 and 1 states; the maximal error on a single
        pixel is of np.pi radians. Hence, we normalise single pixel error by np.pi.
        '''
        result = np.sqrt(np.mean(residualsdata**2, axis=axis)) / np.pi
        return result

    @staticmethod
    def ssim(dataobj, Cone=0.01, Ctwo=0.01):
        '''
            Return the expected value of SSIM for a set of true maps and map estimates;
            as well as the list of SSIM scores for a set of true maps and map estimates over
            which the expectation was calculated.

        Parameter:
        ----------
            dataobj (`type`| .npz object):
                An npz file loaded via np.load() command. We access the following
                files in the compressed object:

                dataobj["macro_predictions"] : A set of estimate posterior f_states
                    produced by an algorithm.

                dataobj["macro_true_fstate"] A set of true maps generated by
                    EngineeredTruth or through experimental measurement with large
                    datasets.

                For both numpy arrays:

                    Axis = 0: Empty axis [NOT USED]
                    Axis = 1: Number of Bayes trials
                    Axis = 2: Vectorised map. Dimensions corresspond to
                        total number of qubits on the grid.

        Returns:
        -------

        result (`type` | float scalar):
            The expected value of SSIM score over all trials.
            This takes values between 0 - ideal (f_state and true map are identical)
            and 1 - maximal error with respect to the SSIM metric.
        scores  (`type` | list object):
            The list of SSIM scores for a set of true maps and map estimates over
            which the expectation was calculated.
        '''
        trials = dataobj["macro_predictions"].shape[1]
        scores = np.zeros(trials)

        for idx_run in range(trials):

            pred = dataobj["macro_predictions"][0, idx_run, :]
            truth = dataobj["macro_true_fstate"][0, idx_run, :]
            scores[idx_run] = Metric.score_ssim(pred, truth, Cone=Cone, Ctwo=Ctwo)

        result = np.mean(scores) # returing deviations abs(1 - ssim)

        return result, scores

    @staticmethod
    def score_ssim(sigx, sigy, Cone=0.01, Ctwo=0.01):
        '''Return SSIM score as a deviation from ideal. Helper function
            to ssim().
            Based on  Equation 13 in Reference [1].

            Note that score_ssim(sigx, sigy) != score_ssim(sigy, sigyx) as the
            ssim metric is not symmetric. Hence, the convention used for all
            analysis is that:

                sigx : vectorised image representing an algorithmic estimated of
                    f_state
                sigy : vectorised image representing the true map

        Parameters:
        -----------
            sigx (`dtype` | numpy array):
                Vectorised image array representing the estimate of f_state.
                sigx is shorthand for Sigma_x, following notation of Reference [1].
            sigy (`dtype` | numpy array):
                Vectorised image array representing the true f_state.
                sigy is shorthand for Sigma_y, following notation of Reference [1].
            Cone (`dtype` : float):
                Numerical stablisers for the SSIM score for data with means
                    and variances close to zero.
                Default value: 0.01.
            Ctwo (`dtype` : float):
                Numerical stablisers for the SSIM score for data with means
                    and variances close to zero.
                Default value: 0.01.

        Return:
        ------
            deviation : Deviation from the ideal ssim score of unity. This enables
                an interpretation where 0 represents the ideal (estimated f_state is exactly
                the true map); and 1 represents maximal error.

        Reference:
        ----------
            [1] Wang, Z., Bovik, A. C., Sheikh, H. R., & Simoncelli, E. P. (2004).
            Image quality assessment: from error visibility to structural similarity.
            IEEE transactions on image processing, 13(4), 600-612.

        '''
        mu_x = np.mean(sigx)
        mu_y = np.mean(sigy)
        stdx = np.std(sigx, ddof=1)
        stdy = np.std(sigx, ddof=1)
        
        Nminus1 = sigx.shape[0] - 1
        covarxy = (1.0 / Nminus1)*np.sum((sigx - mu_x)*(sigy - mu_y))

        ssim = (2.0*mu_x*mu_y + Cone) * (2.0*covarxy + Ctwo)
        ssim /= (mu_x**2 + mu_y**2 + Cone) * (stdx**2 + stdy**2 + Ctwo)

        deviation = abs(1.0 - ssim)  # returing deviation
        return deviation # ssim

    @staticmethod
    def original_err_metric(macro_residuals):
        ''' [NOT USED] '''
        result = (1.0/ np.sqrt(macro_residuals.shape[2]))*np.mean(np.linalg.norm(macro_residuals[0, :, :], axis=1), axis=0)
        return result / np.pi

    @staticmethod
    def err(macro_residuals):
        ''' [NOT USED] '''
        error = np.mean(Metric.singlemap_err(macro_residuals[0, :, :], axis=1), axis=0)
        return error

    @staticmethod
    def singlemap_err(residualsdata, axis):
        ''' [NOT USED] '''
        normaliser = 1.0 / np.sqrt(residualsdata.shape[axis])
        result = normaliser * np.linalg.norm(residualsdata, axis=axis) / np.pi
        return result


################################################################################
# DATA HANDLING FOR ALGORITHM BAYES RISK SIMULATIONS
################################################################################

class DataCube(object):
    '''Class object to handle updates to model parameters for QSLAM
    and EngineeredTruth during simulations. '''
    def __init__(self, loops_dictionary):
        ''' Creates a DataCube instance.

        Parameters:
        -----------
            loops_dictionary (`type` | dictionary object):

                Scans different parameter regimes for physical and algorithmic
                    degrees of freedom in QSLAM simulations.

                Dictionary object with the following keys | elements :

                True Field parameters:
                    meta_truth_floor_scan | List of "OneStepdfloorarea" values to scan.
                        See riskanalysis.EngineeredTruth class definition.
                    truth_step_scan | List of "OneStepdheight" values to scan.
                        See riskanalysis.EngineeredTruth class definition.
                    truth_step_scan | List of "OneStepdheight" values to scan.
                        See riskanalysis.EngineeredTruth class definition.

                Physical Experiment parameters:
                    meta_max_iter_scan | List of "MAX_NUM_ITERATIONS" values to scan.
                        See qslamdesignparams.GLOBALDICT dictionary definition.
                    msmt_per_qubit_scan | List of "MSMTS_PER_NODE" values to scan.
                        See qslamdesignparams.GLOBALDICT dictionary definition.

                Measurement Noise parameters:
                    meta_noisevar_scan | List of "noise_type" values to scan.
                        See engineerednoise.EngineeredNoise class definition.
                    meta_noisevar_scan | List of "prob_hit" values to scan.
                        See engineerednoise.EngineeredNoise class definition.

                QSLAM Free parameters:
                    lambda_scan |  List of [LAMBDA_1", "LAMBDA_2"] values to scan.
                        See qslamdesignparams.GLOBALDICT dictionary definition.

        Class Attributes:
        ----------------
            Identical to dictionary keys in loops_dictionary.

        '''
        for key in loops_dictionary.keys():
            setattr(self, key, loops_dictionary[key])

    def meta_loop_update(self, vardict, vardict_truth, idx_prevar, idx_msmt_var, idx_noise_var,
                         flag='floor'):

        if flag == 'floor':
            vardict_truth["OneStepdfloorarea"] = self.meta_truth_floor_scan[idx_prevar]

        if flag == 'height':
            vardict_truth["OneStepdheight"]["low"] = self.truth_step_scan[idx_prevar][0]
            vardict_truth["OneStepdheight"]["high"] = self.truth_step_scan[idx_prevar][1]

        vardict["MODELDESIGN"]["MAX_NUM_ITERATIONS"] = self.meta_max_iter_scan[idx_msmt_var]

        vardict["ADDNOISE"]["args"]["noise_type"] = self.meta_noisevar_scan[idx_noise_var][0]
        vardict["ADDNOISE"]["args"]["prob_hit"] = self.meta_noisevar_scan[idx_noise_var][1]

        return vardict, vardict_truth

    def inner_loop_update(self, vardict, idx_var, flag='weights'):

        if flag == 'weights':
            vardict["MODELDESIGN"]["LAMBDA_1"] = self.lambda_scan[idx_var][0]
            vardict["MODELDESIGN"]["LAMBDA_2"] = self.lambda_scan[idx_var][1]

        if flag == 'msmtinfo':
            vardict["MODELDESIGN"]["MSMTS_PER_NODE"] = self.msmt_per_qubit_scan[idx_var]

        return vardict


################################################################################
# MATPLOTLIB GLOBAL CONFIGURATIONS
################################################################################

fsize = 8
Fsize = 8
# Set global parameters
matplotlib.rcParams['font.size'] = fsize # global
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = 'Arial'
matplotlib.rcParams['axes.linewidth'] = 0.5

matplotlib.rcParams['mathtext.default'] = 'regular'
# makes mathtext mode Arial. note mathtext is used as ticklabel font in log plots

# Set global tick mark parameters
matplotlib.rcParams['xtick.major.width'] = 0.5
matplotlib.rcParams['ytick.major.width'] = 0.5
matplotlib.rcParams['xtick.labelsize'] = fsize
matplotlib.rcParams['ytick.labelsize'] = fsize
# matplotlib.rcParams['xtick.minor.visible'] = False
# # need to show mior grid as white outline in maps
# matplotlib.rcParams['ytick.minor.visible'] = False
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'

from matplotlib.gridspec import GridSpec as gs
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches
from matplotlib import cm
import matplotlib.ticker as plticker
from matplotlib.colors import ListedColormap
Path = mpath.Path


# DATA HANDLING

GRIDLW = 0.5

DATAVIEWKEYS = {"truth" : "macro_true_fstate",
                "pred_f" : "macro_predictions",
                "pred_r" : "macro_correlations",
                "path" : "macro_paths",
                "errs" : "macro_residuals"
               }

FLAG = {"q" : "qslamr",
        "n" : "naive"
       }

NPZFLAG = {FLAG["q"]: '_BR_Map.npz',
           FLAG["n"]: '_NE_Map.npz'
          }

PATHDICT = {"pdir": './data',
            "fle": 'empty'}

def path_to_file(PATHDICT, flag="q"):
    'Return path to file'
    pathtofile = os.path.join(PATHDICT["pdir"],
                              PATHDICT["fle"] + NPZFLAG[FLAG[flag]])
    return pathtofile

def check_savefig(current_indices, controls=None):
    'Return True if parameter regimes for plotting match controls'
    if controls is None:
        return False
    for item in controls:
        result = np.all(np.asarray(current_indices) == np.asarray(item))
        if result:
            return result

def cm2inch(cm):
    ''' Return values in inches for input values in cm'''
    return cm / 2.54

existingmap = cm.get_cmap('BuPu', 256)
newcolors = existingmap(np.linspace(0.2, 0.8, 256))
newcmap = ListedColormap(newcolors)

PKWG = {"pc": 'k',#'w', # path base color
        "palpha": 1., # path color alpha value
        "pss": 1.5, # path start / stop size
        "ac": 'k', # arrow base color
        "aalpha": 1., # arrow color alpha value
        "adxdy": 0.05, # arrow dx and dy as ratio of path segment
        "alw": 1.0, # arrow line width
        "ahw": 0.05, # arrow head width
       }

HEATMAP = {"vmin" : 0.0,
           "vmax" : np.pi,
           "cmap" : newcmap,# 'pink', 'bone',#'viridis',
           "origin": 'lower'
          }

CORRMAP = {"vmin" : 0.0,
           "vmax" : 30.,
           "cmap" : 'Purples',
           "origin": 'lower'
          }

class qPlotter(object):
    '''Class object for aiding visualisation of Bayes analysis.'''

    def __init__(self,
                 userPKWG=None,
                 userHEATMAP=None,
                 userCORRMAP=None):

        self.PKWG = userPKWG
        self.HEATMAP = userHEATMAP
        self.CORRMAP = userCORRMAP

        if self.PKWG is None:
            self.PKWG = PKWG

        if self.HEATMAP is None:
            self.HEATMAP = HEATMAP

        if self.CORRMAP is None:
            self.CORRMAP = CORRMAP



    def get_single_run(self, dataobj, viewtype, pickone):
        ''' Return a randomly chosen Bayes trial.
        Helper function for plotting results.
        '''

        totalruns = dataobj[DATAVIEWKEYS[viewtype]].shape[1]

        if pickone is None:
            pickone = np.random.randint(low=0, high=totalruns)

        data = dataobj[DATAVIEWKEYS[viewtype]][0, pickone, :]
        return data


    def show_map(self, ax, dataobj, viewtype, pickone=None, linear=False):
        ''' Return the axes handle of a figure where a true_map is plotted
        as a heat map (2D or 1D).
        Helper function for plotting results.
        '''
        if viewtype != 'expt':
            statedata = self.get_single_run(dataobj, viewtype, pickone)

        if viewtype == 'expt':
            statedata = dataobj * 1.0

        if linear is False:
            if statedata.shape[0] < 4:
                cax = ax.imshow(statedata[np.newaxis, :], **self.HEATMAP)
            if statedata.shape[0] >= 4:
                mapdims = int(np.sqrt(statedata.shape))
                mapdata = statedata.reshape(mapdims, mapdims)
                cax = ax.imshow(mapdata, **self.HEATMAP)
        if linear is True:
            mapdims = statedata.shape[0]
            mapdata = np.zeros((1, mapdims))
            for idx in range(1): # make imshow fatter # undone - do with aspect ratio.
                mapdata[idx, : ] = statedata
            cax = ax.imshow(mapdata, **self.HEATMAP)

        # Show all ticks...
        ax.set_xticks(np.arange(mapdata.shape[1]), minor=True)
        ax.set_yticks(np.arange(mapdata.shape[0]), minor=True)

        # Get rid of bounding box
        # for edge, spine in ax.spines.items():
        #    spine.set_visible(False)

        # Make a grid and put labels on the center
        axis_list = [ax.yaxis, ax.xaxis]
        for idx_axis in range(2):

            labels = range(1, mapdata.shape[idx_axis] + 1, 1)
            locs = np.arange(len(labels))

            axis = axis_list[idx_axis]
            axis.set_ticks(locs + 0.5, minor=True)

            if linear == True:
                locs = locs[4::5]
                labels = labels[4::5]

            axis.set(ticks=locs, ticklabels=labels)

        # Make the grid white
        ax.grid(which="minor", color="w", linestyle='-', linewidth=GRIDLW)
        ax.tick_params(axis='both', which='both', color='w')
        ax.tick_params(axis='both', which='major', color='k', direction='out')

        # Turn off ticks in linear plots
        # if linear is True:
            # ax.tick_params(axis='both', which='both', colors='None')

        return ax, cax


    def show_control_path(self, ax,
                          dataobj,
                          GRIDDICT,
                          viewtype="path",
                          linear=False,
                          pickone=None):
        ''' Return the axes handle of a figure where the control path of a single
        QSLAM run is plotted from start to finish using pointed, connected arrows.
        Helper function for plotting results.
        '''
        if viewtype != 'expt':
            controlpath = self.get_single_run(dataobj, viewtype, pickone)

        elif viewtype == 'expt':
            controlpath = dataobj

        points = self.get_control_path(controlpath, GRIDDICT)

        if linear is True:

            # plot 1D path as expanded verticle steps in y axis
            points = [(points[idx][0], idx*-1.0) for idx in np.arange(len(points))]

        codes = [Path.LINETO] * len(points)
        codes[0] = Path.MOVETO

        path = mpath.Path(points, codes)
        patch = mpatches.PathPatch(path,
                                   facecolor='None',
                                   edgecolor=self.PKWG["pc"],
                                   alpha=self.PKWG["palpha"])
        ax.add_patch(patch)

        ax.plot(points[0][0], points[0][1], 'o',
                color=self.PKWG["pc"],
                # edgecolor=self.PKWG["pc"],
                ms=self.PKWG["pss"])
        ax.plot(points[-1][0], points[-1][1], 's',
                color=self.PKWG["pc"],
                # edgecolor=self.PKWG["pc"],
                ms=self.PKWG["pss"])

        for idp in range(len(points)-1):

            ax.arrow(points[idp][0], points[idp][1],
                     self.PKWG["adxdy"]*(points[idp+1][0] - points[idp][0]),
                     self.PKWG["adxdy"]*(points[idp+1][1] - points[idp][1]),
                     shape='full',
                     lw=self.PKWG["alw"],
                     facecolor=self.PKWG["ac"],
                     edgecolor=self.PKWG["ac"],
                     alpha=self.PKWG["aalpha"],
                     length_includes_head=True,
                     width=self.PKWG["ahw"])

        return ax


    def get_control_path(self, path, griddict):
        '''Returns qubit positions (x,y) visited based on a control path and
            input dictionary of qubits

        '''

        result = []

        for qubit in path:

            if qubit + 1 < 10:
                result.append(griddict["QUBIT_0" + str(qubit + 1)])

            if qubit + 1 >= 10:
                result.append(griddict["QUBIT_" + str(qubit + 1)])

        if len(path) == len(result):
            return result

        raise RuntimeError

