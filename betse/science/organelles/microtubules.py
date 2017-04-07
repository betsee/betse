#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

"""
**Microtubules** (i.e., vectors aligned in the direction of endogenous current
flux through cells) functionality.
"""

import numpy as np


class Mtubes(object):
    '''
    Object encapsulating all cellular microtubules for the cell cluster
    simulated at the current time step.

    Attributes
    ----------
    mtubes_x : ndarray
        One-dimensional Numpy array indexing each cell membrane such that each
        element is the normalized X component of the microtubule unit vector
        spatially situated at the midpoint of that membrane for this time step.
    mtubes_y : ndarray
        One-dimensional Numpy array indexing each cell membrane such that each
        element is the normalized Y component of the microtubule unit vector
        spatially situated at the midpoint of that membrane for this time step.
    mtubes_xo : ndarray
        One-dimensional Numpy array indexing each cell membrane such that each
        element is the non-normalized X component of the microtubule vector
        spatially situated at the midpoint of that membrane for this time step.
    mtubes_yo : ndarray
        One-dimensional Numpy array indexing each cell membrane such that each
        element is the non-normalized Y component of the microtubule vector
        spatially situated at the midpoint of that membrane for this time step.
    '''

    def __init__(self, cells, p, alpha_noise = 1.0):

        self.alpha_noise = alpha_noise

        # basic parameters for microtubule units:
        tubulin_mass = 1.66e-27 * 1e3 * 100    # kg
        tubulin_length = 4.5e-9                # m

        # microtubule moment of inertia :
        self.Imit = (tubulin_mass * (cells.R.mean() / tubulin_length)) * (cells.R.mean()) ** 2

        # microtubule dipole moment [C m] (NOTE: dipole moment of microtubule may be as high as 34000 D):
        self.pmit = p.mt_dipole_moment * 3.33e-30   # 1740

        # normalize the cell radial vectors:
        norm_rads = np.sqrt(cells.rads[:, 0] ** 2 + cells.rads[:, 1] ** 2)

        self.radx = cells.rads[:, 0] / norm_rads
        self.rady = cells.rads[:, 1] / norm_rads

        # initialize microtubule dipole vectors:
        self.mtubes_xo = self.radx * self.pmit
        self.mtubes_yo = self.rady * self.pmit

        mtcn = np.sqrt(self.mtubes_xo ** 2 + self.mtubes_yo ** 2)

        # normalized microtubule vectors from the cell centre point:
        self.mtubes_x = self.mtubes_xo / mtcn
        self.mtubes_y = self.mtubes_yo / mtcn


    def update_mtubes(self, cells, sim, p):

        # print(self.pmit/3.33e-30)

        # calculate average field at axial regions of the cell from currents:
        # jjx = (cells.nn_tx*sim.Jn + sim.J_cell_x[cells.mem_to_cells])/2
        # jjy = (cells.nn_ty*sim.Jn + sim.J_cell_y[cells.mem_to_cells])/2

        MM = np.column_stack((self.mtubes_xo, self.mtubes_yo))
        JJ = np.column_stack((sim.J_cell_x[cells.mem_to_cells], sim.J_cell_y[cells.mem_to_cells]))

        torque = np.cross(MM, JJ)

        dtheta = torque / self.Imit

        # noise term:
        noisex = np.random.random(sim.mdl) - 0.5
        noisey = np.random.random(sim.mdl) - 0.5
        noiseTheta = np.arctan2(noisey, noisex)

        theta = np.arctan2(self.mtubes_yo, self.mtubes_xo) + 0.1*dtheta*p.dt - 0.1*noiseTheta*self.alpha_noise*p.dt

        rmit = np.sqrt(self.mtubes_xo ** 2 + self.mtubes_yo ** 2)

        self.mtubes_xo = rmit * np.cos(theta)
        self.mtubes_yo = rmit * np.sin(theta)

        mtcn = np.sqrt(self.mtubes_xo ** 2 + self.mtubes_yo ** 2)

        # normalized microtubule vectors from the cell centre point:
        self.mtubes_x = self.mtubes_xo / mtcn
        self.mtubes_y = self.mtubes_yo / mtcn

    def mtubes_to_cell(self, cells, p):

        # determine the microtubules base electroosmotic velocity:
        uxmto = self.mtubes_x
        uymto = self.mtubes_y

        uxmt = np.dot(cells.M_sum_mems, uxmto*cells.mem_sa)/cells.cell_sa
        uymt = np.dot(cells.M_sum_mems, uymto*cells.mem_sa)/cells.cell_sa

        # uumt = np.sqrt(uxmt**2 + uymt**2)

        return uxmt, uymt

    def remove_mtubes(self, target_inds_mem, cells, sim, p):

        # remove microtubules from the lists and reassign objects:
        mtubesxo2 = np.delete(self.mtubes_xo, target_inds_mem)
        self.mtubes_xo = mtubesxo2*1

        mtubesyo2 = np.delete(self.mtubes_yo, target_inds_mem)
        self.mtubes_yo = mtubesyo2*1

        mtubesx2 = np.delete(self.mtubes_x, target_inds_mem)
        self.mtubes_x = mtubesx2*1

        mtubesy2 = np.delete(self.mtubes_y, target_inds_mem)
        self.mtubes_y = mtubesy2*1


