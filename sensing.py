''' Module: sensing.py

The purpose of this module is to ...

The module contains the following classes:
Robot(object):
Scanner(Map): Inherits from Robot.
'''

import numpy as np
import sys
sys.path.append('../../GIT/')
from qif.common import projected_msmt

class Robot(object):
    '''docstring'''
    def __init__(self, x=0, y=0, corr_r=0, sigma=0, R=0, phi=None, 
                 localgridcoords=[1,1]):

        self.r_pose = np.asarray([x, y, corr_r])
        self.r_msmtnoise = R
        self.r_motionnoise = sigma

        self.r_dims = len(self.r_pose)
        self.r_phi = np.eye(self.r_dims) if phi is None else phi

        self.r_questbk = np.zeros(localgridcoords)
        self.r_guestbk_counter = np.zeros(localgridcoords)

    def r_xy(self):
        '''docstring'''
        return (self.r_pose[0], self.r_pose[1])

    def r_move(self, u_cntrl, noisy=False):
        '''docstring'''
        dynamics = np.dot(self.r_pose, self.r_phi) + u_cntrl

        noise = np.zeros(self.r_dims)
        if noisy:
            stdev = np.sqrt(self.r_motionnoise)
            noise = np.eye(self.r_dims)*np.random.normal(loc=0,
                                                        scale=stdev)
            noise[self.r_dims-1] = 0.0 # corr_r unaaffected

        self.r_pose = dynamics + noise
        return self.r_pose, dynamics, noise

    def r_measure(self, mapval):
        '''docstring'''
        msmt = projected_msmt([mapval]) # imported from qif. mapval can be a list
        return msmt

    def r_addguest(self, pos_x, pos_y, msmt):
        '''Updates the guestbook of physical measurements at each nodes'''
        n_x, n_y = self.return_position(pos_x), self.return_position(pos_y)

        old_prob = self.r_questbk[n_x, n_y]*1.0
        old_counter = self.r_guestbk_counter[n_x, n_y]*1.0
        self.r_questbk[n_x, n_y] = (old_prob*old_counter + msmt)/(old_counter+1)
        self.r_guestbk_counter[n_x, n_y] += 1
    
    @staticmethod
    def return_position(pos):
        print pos
        return int(pos)
    
class Scanner(Robot):
    '''docstring
    '''

    def __init__(self, localgridcoords, x=0, y=0, corr_r=1.42, sigma=0, R=0,
                 phi=None):
        '''
        Defines a 'measurement scan' that a robot makes. This consists of a
        physical measurement taken at the robot position (x,y), and quasi
        msmts taken on k nearest neighbours to (x,y). The quasi-msmt procedure
        and the size of the neighbourhood depends on an effective correlation
        length, corr_r.

        localgrid: assumed known and unchanging (array of qubits on hardware)
        corr_r = 1.42 : For k=1, nearest diagnoal neighbours are included
        '''
        Robot.__init__(self,
                       localgridcoords=localgridcoords,
                       x=x,
                       y=y,
                       corr_r=corr_r,
                       sigma=sigma,
                       R=R,
                       phi=phi)

    def r_corr_length(self, v_sep, type='Gaussian'):
        '''Returns the correlation strenght between physical and quasi
        measurements based on a separation distance, v_sep
        '''
        return np.exp(-(v_sep)**2 / 2.0*self.r_pose[2]**2)

    def r_get_quasi_msmts(self, knn_list, born_est):
        '''Returns quasi measurements on the neighbours of the robot,
        using the past measurement record at the robot pose
        '''
        quasi_msmts = []

        for neighbour in knn_list:
            vsep = np.linalg.norm(np.subtract(neighbour, self.r_xy()))
            quasi_born = self.r_corr_length(vsep)*born_est
            quasi_msmts.append(self.r_measure(quasi_born))
        return quasi_msmts

    def r_scan_local(self, mapval, knn_list):
        '''This function takes Born probability estimate based on physical
        msmts in the guestbook, and then blurs each msmt on its nearest
        neighbours to approximate a scan
        '''

        pose_x, pose_y = self.r_xy()
        print("PRINT POSITION", pose_x, pose_y)

        msmt = self.r_measure(mapval)
        self.r_addguest(pose_x, pose_y, msmt)
        born_est = self.r_questbk[pose_x, pose_y]

        scan_msmts = [msmt]  +  self.r_get_quasi_msmts(knn_list, born_est)
        scan_posxy = [self.r_xy()] + knn_list

        return zip(scan_posxy, scan_msmts)

    