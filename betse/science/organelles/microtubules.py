#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

"""
**Microtubules** (i.e., vectors aligned in the direction of endogenous current
flux through cells) functionality.
"""

import numpy as np
from betse.exceptions import BetseSimInstabilityException
from betse.util.io.log import logs
from betse.science.math import modulate as mods
# from betse.science.math import finitediff as fd
from scipy.ndimage.filters import gaussian_filter


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

    def __init__(self, sim, cells, p):

        # basic parameters for microtubule units:
        self.rmt = p.mt_radius # radius of a microtubule

        self.L = p.mt_length*(cells.R_rads/cells.R_rads.mean())

        self.tubulin_N = (222/1.0e-6)*self.L # total number of tubulin molecules per microtubule

        self.charge_mtube = (p.length_charge*p.q*self.L)/1.0e-6  # charge on the microtubule in Coulombs

        self.visc = p.cytoplasm_viscocity   # viscocity of cytoplasm

        # microtubule dipole moment [C m] (assumed to be sum of tubulin dimer dipole moments, which are 1740 D each):
        self.p_mtube = self.tubulin_N*p.tubulin_dipole*3.33e-30

        self.p_ind = self.tubulin_N*p.tubulin_polar*1.0e-40 # inducible polarizability of tubulin C m2/V

        # calculate drag co-efficients (from drag force on cylinders Broersma et al, Hunt et al 1994)
        self.drag_C(cells)

        # radial drag force on the microtubule from Stokes-Einstein relation "Rotational Diffusion" theory:
        self.drag_r = self.C_rad*p.cytoplasm_viscocity*(self.L**3)

        # microtubule rotational diffusion constant:
        self.Dr = ((p.kb*p.T)/self.drag_r)

        # initial angle of microtubules:
        self.mt_theta = np.arctan2(cells.mem_vects_flat[:,3], cells.mem_vects_flat[:,2])

        # normalized microtubule vectors from the cell centre point:
        self.mtubes_x = cells.mem_vects_flat[:,2]
        self.mtubes_y = cells.mem_vects_flat[:,3]

        # microtubule density function initialized:
        mtdx = np.dot(cells.M_sum_mems, self.mtubes_x*cells.mem_sa) / cells.cell_sa
        mtdy = np.dot(cells.M_sum_mems, self.mtubes_y*cells.mem_sa) / cells.cell_sa

        self.mtdf = ((mtdx[cells.mem_to_cells]*cells.mem_vects_flat[:,2] +
                                         mtdy[cells.mem_to_cells]*cells.mem_vects_flat[:,3]))

        self.modulator = np.ones(sim.mdl)  # initialize a modulator structure for the microtubule dynamics

        # Initialize microtubules with mini-simulation to define initial state

        gXo = np.zeros(sim.cdl)
        gYo = np.zeros(sim.cdl)

        if p.init_mtx != 'None' and p.init_mtx is not None:
            Ixo, _ = getattr(mods, p.init_mtx)(cells.cell_i, cells, p)

            gXo = (Ixo / Ixo.max())

            gX = gXo[cells.mem_to_cells]

        else:
            p.init_mtx = None
            gX = gXo[cells.mem_to_cells]

        if p.init_mty != 'None' and p.init_mty is not None:
            Iyo, _ = getattr(mods, p.init_mty)(cells.cell_i, cells, p)

            gYo = (Iyo / Iyo.max())

            gY = gYo[cells.mem_to_cells]

        else:
            p.init_mty = None
            gY = gYo[cells.mem_to_cells]

        if p.init_mtx is not None or p.init_mty is not None:

            logs.log_info("Running initialization sequence for microtubules...")

            for tt in range(0, 300):
                # cross-product with director
                Fmt = (self.mtubes_x*gY - self.mtubes_y*gX)/(2*np.pi)

                # calculate rotational flux of microtubule:
                gc_mtdf = np.dot(cells.gradTheta, self.mtdf)

                flux_theta = Fmt - gc_mtdf*1.0e-8
                # flux_theta = Fmt

                self.mt_theta +=  flux_theta*0.15

                # recalculate the new cell microtubule density function:
                mtdx = np.dot(cells.M_sum_mems, self.mtubes_x * cells.mem_sa) / cells.cell_sa
                mtdy = np.dot(cells.M_sum_mems, self.mtubes_y * cells.mem_sa) / cells.cell_sa

                self.mtdf = (mtdx[cells.mem_to_cells] * cells.mem_vects_flat[:, 2] +
                        mtdy[cells.mem_to_cells] * cells.mem_vects_flat[:, 3])

                self.mtubes_x = np.cos(self.mt_theta)
                self.mtubes_y = np.sin(self.mt_theta)

    def drag_C(self, cells):

        # calculate drag co-efficients for the microtubules:

        v = 1 / (np.log(self.L/self.rmt))

        g_para = -0.114 - 0.15 * v - 13.5 * (v ** 2) + 37 * (v ** 3) - 22 * (v ** 4)
        g_perp = 0.866 - 0.15 * v - 8.1 * (v ** 2) + 18 * (v ** 3) - 9 * (v ** 4)
        g_rad = -0.446 - 0.2 * v - 16 * (v ** 2) + 63 * (v ** 3) - 62 * (v ** 4)

        self.C_para = (2 * np.pi) / (np.log(self.L / (2 * self.rmt)) + g_para)
        self.C_perp = (4 * np.pi) / (np.log(self.L / (2 * self.rmt)) + g_perp)
        self.C_rad = ((1 / 3) * np.pi) / (np.log(self.L / (2 * self.rmt)) + g_rad)

    def reinit(self, cells, p):

        """
        Reinitialize key microtubule parameters that may have changed in config file since init.

        """

        # basic parameters for microtubule units:
        self.rmt = p.mt_radius # radius of a microtubule

        self.L = p.mt_length*(cells.R_rads/cells.R_rads.mean()) # length of microtubule

        self.tubulin_N = (222/1.0e-6)*self.L # total number of tubulin molecules per microtubule

        self.charge_mtube = (p.length_charge*p.q*self.L)/1.0e-6  # charge on the microtubule in Coulombs

        self.visc = p.cytoplasm_viscocity   # viscocity of cytoplasm

        # microtubule dipole moment [C m] (assumed to be sum of tubulin dimer dipole moments, which are 1740 D each):
        self.p_mtube = self.tubulin_N*p.tubulin_dipole*3.33e-30

        self.p_ind = self.tubulin_N*p.tubulin_polar*1.0e-40 # inducible polarizability of tubulin C m2/V

        # calculate drag co-efficients (from drag force on cylinders Broersma et al, Hunt et al 1994)
        self.drag_C(cells)

        # radial drag force on the microtubule from Stokes-Einstein relation "Rotational Diffusion" theory:
        self.drag_r = self.C_rad*p.cytoplasm_viscocity*(self.L**3)

        # microtubule rotational diffusion constant:
        self.Dr = ((p.kb*p.T)/self.drag_r)

    def update_mtubes(self, cells, sim, p):

        # microtubule radial vectors:
        ui = self.mtubes_x*self.L
        vi = self.mtubes_y*self.L

        ui_hat = self.mtubes_x
        vi_hat = self.mtubes_y

        Ex = sim.E_cell_x[cells.mem_to_cells]
        Ey = sim.E_cell_y[cells.mem_to_cells]

        gEx = (Ex[cells.nn_i] - Ex[cells.mem_i]) / (cells.nn_len)

        gExx = gEx*cells.nn_tx
        gExy = gEx*cells.nn_ty

        gEy = (Ey[cells.nn_i] - Ey[cells.mem_i]) / (cells.nn_len)

        gEyx = gEy*cells.nn_tx
        gEyy = gEy*cells.nn_ty

        q_tube = self.charge_mtube

        # if p.fluid_flow:
        #
        #     # force from any fluid flow mapped to membranes:
        #     Fdux = sim.u_cells_x[cells.mem_to_cells]*self.C_perp*self.L*p.cytoplasm_viscocity
        #     Fduy = sim.u_cells_y[cells.mem_to_cells]*self.C_perp*self.L*p.cytoplasm_viscocity
        #
        # else:

        Fdux = np.zeros(sim.mdl)
        Fduy = np.zeros(sim.mdl)

        if p.tethered_tubule is False:

            torque_tether = np.zeros(sim.mdl)

            # gradient of the field will torque the monopole by applying different forces at ends:
            torque_gradient = (q_tube * (ui) * (gEyx.ravel() * ui + gEyy.ravel() * vi) -
                               q_tube * (vi) * (gExx.ravel() * ui + gExy.ravel() * vi))

            # fiber will also align such that ends are at the same voltage:
            torque_monopole = (q_tube*(ui)*Ex.ravel() + q_tube*(vi)*Ey.ravel()) + (ui*Fduy.ravel() + vi*Fdux.ravel())

            # fiber will also align via its dipole in the electric field:
            torque_dipole = -(self.p_ind * ui_hat * Ey.ravel() - self.p_ind * vi_hat * Ex.ravel())


        # if fiber is tethered, any perpendicular force will represent a torque:
        else:

            torque_tether = (q_tube * ui * Ey.ravel() - q_tube * vi * Ex.ravel()) + (ui*Fduy.ravel() - vi*Fdux.ravel())

            # gradient of the field will torque the monopole by applying different forces at ends:
            torque_gradient = np.zeros(sim.mdl)

            # fiber will also align such that ends are at the same voltage:
            torque_monopole = np.zeros(sim.mdl)

            # fiber will also align via its dipole in the electric field:
            torque_dipole = (self.p_ind * ui_hat * Ey.ravel() - self.p_ind * vi_hat * Ex.ravel())

        flux_theta = (
            + (torque_tether / self.drag_r)
            + (torque_dipole / self.drag_r)
            + (torque_monopole / self.drag_r)
            + (torque_gradient / self.drag_r)
            + ((p.kb * p.T) / self.drag_r) * (0.5 - np.random.rand(len(self.mt_theta)))
        )

        # update the microtubule angle:
        self.mt_theta = self.mt_theta + flux_theta*p.dt*p.dilate_mtube_dt*self.modulator

        # update the microtubule coordinates with the new angle:
        self.mtubes_x = np.cos(self.mt_theta)
        self.mtubes_y = np.sin(self.mt_theta)

    def mtubes_to_cell(self, cells, p):

        # determine the microtubules base electroosmotic velocity:
        uxmto = self.mtubes_x
        uymto = self.mtubes_y

        uxmt = np.dot(cells.M_sum_mems, uxmto*cells.mem_sa)/cells.cell_sa
        uymt = np.dot(cells.M_sum_mems, uymto*cells.mem_sa)/cells.cell_sa

        return uxmt, uymt

    def remove_mtubes(self, target_inds_mem, target_inds_cell, cells, sim, p):


        mtubesx2 = np.delete(self.mtubes_x, target_inds_mem)
        self.mtubes_x = mtubesx2*1

        mtubesy2 = np.delete(self.mtubes_y, target_inds_mem)
        self.mtubes_y = mtubesy2*1

        mtdf2 = np.delete(self.mtdf, target_inds_mem)
        self.mtdf = mtdf2*1

        th2 = np.delete(self.mt_theta, target_inds_mem)
        self.mt_theta = th2*1

        mod2 = np.delete(self.modulator, target_inds_mem)
        self.modulator = mod2*1

        L2 = np.delete(self.L, target_inds_mem)
        self.L = L2*1

        cmt2 = np.delete(self.charge_mtube, target_inds_mem)
        self.charge_mtube = cmt2*1

        tn2 = np.delete(self.tubulin_N, target_inds_mem)
        self.tubulin_N = tn2*1

        pi2 = np.delete(self.p_ind, target_inds_mem)
        self.p_ind = pi2*1

        dra2 = np.delete(self.drag_r, target_inds_mem)
        self.drag_r = dra2*1

        cpa2 = np.delete(self.C_perp, target_inds_mem)
        self.C_perp = cpa2*1

        d2 = np.delete( self.Dr, target_inds_mem)
        self.Dr = d2*1



