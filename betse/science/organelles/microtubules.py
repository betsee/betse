#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

"""
**Microtubules** (i.e., vectors aligned in the direction of endogenous current
flux through cells) functionality.
"""

import numpy as np
from betse.exceptions import BetseSimInstabilityException
from betse.science import sim_toolbox as stb


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

        self.alpha_noise = p.mtube_noise

        # basic parameters for microtubule units:
        tubulin_mass = 1.66e-27 * 1e3 * 100    # kg
        tubulin_length = 4.5e-9                # m

        # microtubule moment of inertia :
        self.Imit = (tubulin_mass * (cells.R.mean() / tubulin_length)) * (cells.R.mean()) ** 2

        # microtubule dipole moment [C m] (NOTE: dipole moment of microtubule may be as high as 34000 D):
        self.pmit = p.mt_dipole_moment * 3.33e-30   # 1740

        # microtubule end charge:
        self.z = (self.pmit/(cells.R_rads*p.q))

        # total microtubule concentration (determines sensitivity to ordering)
        if p.mtube_noise <= 0.0:

            p.mtube_noise = 1.0e-4

        self.sensitivity = (1/p.mtube_noise)

        # microtubule "plus end" concentration initialized:
        self.cp = np.ones(len(cells.mem_i))

        # microtubule "negative end" concentration initialized:
        self.cn = np.ones(len(cells.mem_i))

        # microtubule diffusion constant:
        self.D = p.D_mtube

        # normalize the cell radial vectors:
        norm_rads = np.sqrt(cells.rads[:, 0] ** 2 + cells.rads[:, 1] ** 2)

        self.radx = cells.rads[:, 0] / norm_rads
        self.rady = cells.rads[:, 1] / norm_rads

        self.mt_n = np.zeros(len(cells.mem_i))

        # initialize microtubule dipole vectors:
        self.mtubes_xo = self.radx * self.pmit
        self.mtubes_yo = self.rady * self.pmit

        # normalized microtubule vectors from the cell centre point:
        self.mtubes_x = cells.mem_vects_flat[:,2]
        self.mtubes_y = cells.mem_vects_flat[:,3]


    def reinit(self, cells, p):

        # microtubule diffusion constant:
        self.D = p.D_mtube

        # microtubule dipole moment [C m] (NOTE: dipole moment of microtubule may be as high as 34000 D):
        self.pmit = p.mt_dipole_moment * 3.33e-30  # 1740

        # microtubule end charge:
        self.z = (self.pmit / (cells.R_rads * p.q))


        if p.mtube_noise <= 0.0:

            p.mtube_noise = 1.0e-4

        # total microtubule concentration (determines sensitivity to ordering)
        self.sensitivity = (1/p.mtube_noise)

    def update_mtubes_o(self, cells, sim, p):

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

    def update_mtubes(self, cells, sim, p):

        cav = 1.0  # concentration at cell centre
        cpi = self.cp  # concentration at membrane
        z = self.z  # charge of ion
        Do = self.D  # diffusion constant of ion

        cap = (cav + cpi) / 2  # concentration at midpoint between cell centre and membrane
        cgp = (cpi - cav) / cells.R_rads  # concentration gradients

        cfluxpo = -Do*cgp + ((Do * p.q * cav * z)/(p.kb * sim.T))*sim.Ec

        # as no net mass must leave this intracellular movement, make the flux divergence-free:
        cfluxp = stb.single_cell_div_free(cfluxpo, cells)

        # calculate the actual concentration at membranes by unpacking to concentration vectors:
        self.cp = cpi + cfluxp*(cells.mem_sa/cells.mem_vol)*p.dt


        #-----calculate a "negative end" concentration that has equal and opposite value of z:
        can = (1.0 + self.cn) / 2  # concentration at midpoint between cell centre and membrane
        cgn = (self.cn - 1.0) / cells.R_rads  # concentration gradients

        cfluxno = -Do*cgn - ((Do * p.q * cav * z)/(p.kb * sim.T))*sim.Ec

        # as no net mass must leave this intracellular movement, make the flux divergence-free:
        cfluxn = stb.single_cell_div_free(cfluxno, cells)

        # calculate the actual concentration at membranes by unpacking to concentration vectors:
        self.cn = self.cn + cfluxn*(cells.mem_sa/cells.mem_vol)*p.dt

        # deal with the fact that our coarse diffusion model may leave some sub-zero concentrations:
        indsZ = (self.cp < 0.0).nonzero()

        if len(indsZ[0]):
            raise BetseSimInstabilityException(
                "A microtubule calculation has lead to simulation instability.")

        # define microtubule direction vectors in terms of density difference between plus and central minus end:
        # component normal to the membrane:

        mtno = (self.cp - self.cn)*self.sensitivity

        mtx = np.dot(cells.M_sum_mems, mtno*cells.mem_vects_flat[:, 2]*cells.mem_sa) / cells.cell_sa
        mty = np.dot(cells.M_sum_mems, mtno*cells.mem_vects_flat[:, 3]*cells.mem_sa) / cells.cell_sa

        self.mtubes_xo = cells.mem_vects_flat[:, 2] + mtx[cells.mem_to_cells]
        self.mtubes_yo = cells.mem_vects_flat[:, 3] + mty[cells.mem_to_cells]

        mtmag = np.sqrt(sim.mtubes.mtubes_xo ** 2 + sim.mtubes.mtubes_yo ** 2)

        mtmag[mtmag == 0.0] = 1.0

        # normalized microtubule vectors from the cell centre point:
        self.mtubes_x = self.mtubes_xo / mtmag
        self.mtubes_y = self.mtubes_yo / mtmag

    def mtubes_to_cell(self, cells, p):

        # determine the microtubules base electroosmotic velocity:
        uxmto = self.mtubes_x
        uymto = self.mtubes_y

        uxmt = np.dot(cells.M_sum_mems, uxmto*cells.mem_sa)/cells.cell_sa
        uymt = np.dot(cells.M_sum_mems, uymto*cells.mem_sa)/cells.cell_sa

        # uumt = np.sqrt(uxmt**2 + uymt**2)

        return uxmt, uymt

    def remove_mtubes(self, target_inds_mem, target_inds_cell, cells, sim, p):

        # remove microtubules from the lists and reassign objects:
        mtubesxo2 = np.delete(self.mtubes_xo, target_inds_mem)
        self.mtubes_xo = mtubesxo2*1

        mtubesyo2 = np.delete(self.mtubes_yo, target_inds_mem)
        self.mtubes_yo = mtubesyo2*1

        mtubesx2 = np.delete(self.mtubes_x, target_inds_mem)
        self.mtubes_x = mtubesx2*1

        mtubesy2 = np.delete(self.mtubes_y, target_inds_mem)
        self.mtubes_y = mtubesy2*1

        cp2 = np.delete(self.cp, target_inds_mem)
        self.cp = cp2*1


        cn2 = np.delete(self.cn, target_inds_mem)
        self.cn = cn2*1


        z2 = np.delete(self.z, target_inds_mem)
        self.z = z2*1


