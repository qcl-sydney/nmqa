'''
Created on Thu Sep 28 17:40:43 2017
@author: riddhisw

.. module:: EngineeredNoise

    :synopsis: Adds engineered noise to ideal single qubit measurements.

    Module Level Classes:
    ----------------------
        EngineeredNoise : Creates an EngineeredNoise class object to add
        non Gaussian bit-flip noise to ideal qubit (0 or 1) measurements.
        Supported noise types: salt and pepper noise, dark noise.

.. moduleauthor:: Riddhi Gupta <riddhi.sw@gmail.com>
'''
import numpy as np

class EngineeredNoise(object):
    ''' Creates an EngineeredNoise class object to add
        non Gaussian bit-flip noise to ideal qubit (0 or 1) measurements.
        Supported noise types: salt and pepper noise, dark noise.

    Static Methods:
    --------------
    alwaysdark : Return dark noise (sets any input to zero).
    spnoise : Return salt and pepper noise, with user defined dark probability (lightdarksplit)
    noiseless : Return ideal (noiseless) measurement

    Class Methods:
    --------------
    add_noise : Return noisy measurements from input ideal measurements and engineered noise type.
    '''

    def __init__(self):
        '''Creates a noise class object'''
        self.NOISE = {"alwaysdark": EngineeredNoise.alwaysdark,
                      "spnoise": EngineeredNoise.spnoise,
                      "noiseless":EngineeredNoise.noiseless}

    @staticmethod
    def alwaysdark(msmt):
        '''Return dark noise (sets any input to zero).'''
        return np.zeros_like(msmt)

    @staticmethod
    def spnoise(msmt, lightdarksplit=0.5):
        '''Return salt and pepper noise, with user defined dark probability (lightdarksplit).

        Parameters:
        ----------
        lightdarksplit: Dark probability (a value between 0 and 1). Default value 0.5
        '''

        rand = np.random.uniform(0, 1, size=1)

        if rand <= lightdarksplit:
            return EngineeredNoise.alwaysdark(msmt)

        return np.ones_like(msmt)

    @staticmethod
    def noiseless(msmt):
        '''Return a noiseless ideal measurement.'''
        return msmt


    def add_noise(self, msmts,
                  prob_hit=0.1,
                  noise_type='noiseless'):
        ''' Return noisy measurements from input ideal measurements and engineered noise.

        Parameters:
        ----------
        msmts : Input array of ideal single qubit measurement (0 or 1).
        noise_type : Type of engineered noise ('alwaysdark', 'spnoise' or 'noiseless')
                     Default value: 'noiseless'.
        prob_hit : Engineered noise strength (0 - noiseless; 1 - maximal noise)
                   Default value: '0.1'.
        Return:
        ------
        msmts : noisy measurements affected by EngineeredNoise
        '''

        msmts = np.asarray(msmts, dtype=float)
        original_shape = msmts.shape
        msmts = msmts.flatten()

        totalmsmts = msmts.shape[0]
        rand = np.random.uniform(0, 1, size=totalmsmts)

        for idx_msmt in range(msmts.shape[0]):

            if rand[idx_msmt] <= prob_hit:

                msmts[idx_msmt] = self.NOISE[noise_type](msmts[idx_msmt])

        return msmts.reshape(original_shape)
