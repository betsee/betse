#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

"""
**Microtubules** (i.e., vectors aligned in the direction of endogenous current
flux through cells) functionality.
"""

# ....................{ IMPORTS                            }....................
import numpy as np
# from betse.exceptions import BetseSimUnstableException
from betse.util.io.log import logs
from betse.science.math import modulate as mods
# from betse.science.math import finitediff as fd
# from scipy.ndimage.filters import gaussian_filter

# ....................{ CLASSES                            }....................
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

        # inhomogeneous density:
        if p.mt_space_density is not None:

            self.mt_density, _ = mods.gradient_bitmap(cells.mem_i, cells, p, bitmap_filename=p.mt_space_density)

        else:

            self.mt_density = 1.0

        # initial angle of microtubules:
        self.mt_theta = np.arctan2(cells.mem_vects_flat[:,3], cells.mem_vects_flat[:,2])

        # normalized microtubule vectors from the cell centre point:
        self.mtubes_x = cells.mem_vects_flat[:,2]*self.mt_density
        self.mtubes_y = cells.mem_vects_flat[:,3]*self.mt_density

        # microtubule density function initialized:
        mtdx = np.dot(cells.M_sum_mems, self.mtubes_x*cells.mem_sa) / cells.cell_sa
        mtdy = np.dot(cells.M_sum_mems, self.mtubes_y*cells.mem_sa) / cells.cell_sa

        self.mtdf = ((mtdx[cells.mem_to_cells]*cells.mem_vects_flat[:,2] +
                                         mtdy[cells.mem_to_cells]*cells.mem_vects_flat[:,3]))

        self.modulator = np.ones(sim.mdl) # initialize modulator as chemical factor altering microtubule dynamics

        #-----------------------------------------------------------------------------------------------------------

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

        # smoothing weights for membrane and central values:
        nfrac = p.smooth_cells
        self.smooth_weight_mem = ((nfrac*cells.num_mems[cells.mem_to_cells] -1)/(nfrac*cells.num_mems[cells.mem_to_cells]))
        self.smooth_weight_o = 1/(nfrac*cells.num_mems[cells.mem_to_cells])

        # amt = self.L/(2*self.rmt)
        # self.drag_r = (np.pi*p.cytoplasm_viscocity*(self.L**3))/(np.log(amt - 0.662 + (0.917/amt) - (0.05/amt**2)))

        # microtubule rotational diffusion constant:
        self.Dr = ((p.kb*p.T)/self.drag_r)

        # inhomogeneous density:
        if p.mt_space_density is not None:

            self.mt_density, _ = mods.gradient_bitmap(cells.mem_i, cells, p, bitmap_filename=p.mt_space_density)
            # self.mt_denstiy = self.mt_density[cells.mem_to_cells]


        else:

            self.mt_density = 1.0

        if p.mtube_init_x is not None:

            # Initialize microtubules with mini-simulation to define initial state, which may include
            self.Phi_orient_x, _ = mods.gradient_bitmap(cells.mem_i, cells, p, bitmap_filename=p.mtube_init_x)

        else:
            # self.Phi_orient_x = np.zeros(len(cells.mem_i))
            self.Phi_orient_x = cells.mem_vects_flat[:, 2]

        if p.mtube_init_y is not None:

            self.Phi_orient_y, _ = mods.gradient_bitmap(cells.mem_i, cells, p, bitmap_filename=p.mtube_init_y)

        else:
            self.Phi_orient_y = cells.mem_vects_flat[:, 3]
            # self.Phi_orient_y = np.zeros(len(cells.mem_i))
        # if self.Phi_orient_x != cells.mem_vects_flat[:, 2] or self.Phi_orient_y != cells.mem_vects_flat[:,3]:

        # if there's a user-requested gradient for initial condition alignment, then perform the pre-simulation:
        self.presim(self.Phi_orient_x, self.Phi_orient_y, cells, p)


        self.uxmt, self.uymt = self.mtubes_to_cell(cells, p)

    def update_mtubes(self, cells, sim, p):

        # microtubule radial vectors:
        ui = np.cos(self.mt_theta)*self.L
        vi = np.sin(self.mt_theta)*self.L

        # nx = cells.mem_vects_flat[:,2]
        # ny = cells.mem_vects_flat[:,3]

        Ex = sim.E_cell_x[cells.mem_to_cells]
        Ey = sim.E_cell_y[cells.mem_to_cells]

        gEx = (Ex[cells.nn_i] - Ex[cells.mem_i]) / (cells.nn_len)

        gExx = gEx*cells.nn_tx
        gExy = gEx*cells.nn_ty

        gEy = (Ey[cells.nn_i] - Ey[cells.mem_i]) / (cells.nn_len)

        gEyx = gEy*cells.nn_tx
        gEyy = gEy*cells.nn_ty

        q_tube = self.charge_mtube


        if p.tethered_tubule is False:

            torque_tether = np.zeros(sim.mdl)

            # gradient of the field will torque the monopole by applying different forces at ends:
            torque_gradient = (q_tube * (ui) * (gEyx.ravel() * ui + gEyy.ravel() * vi) -
                               q_tube * (vi) * (gExx.ravel() * ui + gExy.ravel() * vi))

            # fiber will also align such that ends are at the same voltage:
            torque_monopole = (q_tube*(ui)*Ex.ravel() + q_tube*(vi)*Ey.ravel())


        # if fiber is tethered, any perpendicular force will represent a torque:
        else:

            torque_tether = (q_tube * ui * Ey.ravel() - q_tube * vi * Ex.ravel())

            # gradient of the field will torque the monopole by applying different forces at ends:
            torque_gradient = np.zeros(sim.mdl)

            # fiber will also align such that ends are at the same voltage:
            torque_monopole = np.zeros(sim.mdl)

            # fiber will also align via its dipole in the electric field:
            # torque_dipole = (self.p_ind * ui_hat * Ey.ravel() - self.p_ind * vi_hat * Ex.ravel())

        flux_theta = (
            + (torque_tether / self.drag_r)
            + (torque_monopole / self.drag_r)
            + (torque_gradient / self.drag_r)
            # + ((p.kb * p.T) / self.drag_r) * (0.5 - np.random.rand(len(self.mt_theta)))
        )

        # update the microtubule coordinates with the new angle:
        if p.dilate_mtube_dt > 0.0:
            # normalized correlation length of the microtubules
            lenmt = np.sqrt(self.uxmt[cells.mem_to_cells] ** 2 + self.uymt[cells.mem_to_cells] ** 2) + p.kb * sim.T

            stdev = np.sqrt(2 * p.dt * p.dilate_mtube_dt * ((p.kb * p.T) / self.drag_r) * lenmt * self.L ** 2)

            noise = np.random.normal(loc=0.0, scale=stdev, size=sim.mdl)

            # update the microtubule angle:
            self.mt_theta = self.mt_theta + flux_theta * p.dt * p.dilate_mtube_dt * self.modulator + noise

            mtubes_xo = np.cos(self.mt_theta)*self.mt_density
            mtubes_yo = np.sin(self.mt_theta)*self.mt_density

            self.mtubes_x, self.mtubes_y = cells.single_cell_div_free(mtubes_xo, mtubes_yo)

            self.uxmt, self.uymt = self.mtubes_to_cell(cells, p)

    def mtubes_to_cell(self, cells, p):

        # determine the microtubules base electroosmotic velocity:
        uxmto = self.mtubes_x
        uymto = self.mtubes_y

        # averages of mtube field at the cell centres:
        # uxmt = (np.dot(cells.M_sum_mems, uxmto*cells.mem_sa)/cells.cell_sa)
        # uymt = (np.dot(cells.M_sum_mems, uymto*cells.mem_sa)/cells.cell_sa)

        uxmt = (np.dot(cells.M_sum_mems, uxmto)/cells.num_mems)
        uymt = (np.dot(cells.M_sum_mems, uymto)/cells.num_mems)

        # average the mtube field to the centre of pie-shaped midpoints of each individual cell:
        # uxmti = (uxmt[cells.mem_to_cells] + uxmto)/2
        # uymti = (uymt[cells.mem_to_cells] + uymto)/2

        uxmti = uxmt[cells.mem_to_cells]
        uymti = uymt[cells.mem_to_cells]

        # uxmti, uymti = cells.single_cell_div_free(uxmti, uymti)

        # Store the normal component of microtubule alignment field mapped to membranes:
        self.umtn = (uxmti * cells.mem_vects_flat[:, 2] + uymti * cells.mem_vects_flat[:, 3])


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

        mtd2 = np.delete( self.mt_density, target_inds_mem)
        self.mt_density = mtd2*1

        # mtux2 = np.delete(self.uxmt, target_inds_cell)
        # self.uxmt = mtux2*1
        #
        # mtuy2 = np.delete(self.uymt, target_inds_cell)
        # self.uymt = mtuy2*1

        # mtuu2 = np.delete(self.umtn, target_inds_cell)
        # self.umtn = mtuu2*1

    def presim(self, gFxo, gFyo, cells, p):

        """
        Initialize microtubule orientation using gradient of a vector field
        :param sim:
        :param cells:
        :return:
        """


        if len(gFxo) == len(cells.cell_i):
            gFxo = gFxo[cells.mem_to_cells]
            gFyo = gFyo[cells.mem_to_cells]

        if p.mtube_init_x is not None and p.mtube_init_y is not None:

            mssg = ("Preinitializing microtubule x- and y- coorinates with {} and {}").format(p.mtube_init_x,
                                                                                              p.mtube_init_y)
            logs.log_info("-------------------------------")
            logs.log_info(mssg)
            logs.log_info("-------------------------------")

            # rotate the axis of the model:
            rotangle = (p.mtube_init_rotangle*np.pi)/180.0

            gFx = gFxo*np.cos(rotangle) - gFyo*np.sin(rotangle)
            gFy = gFxo*np.sin(rotangle) + gFyo*np.cos(rotangle)

            # magnitude of the orienting field:
            magF = (np.sqrt(gFx ** 2 + gFy ** 2)).max() + 1.0e-15

            # set the microtubule vectors with the field values:
            mtubes_xo = -(gFx/magF) * self.mt_density
            mtubes_yo = -(gFy/magF) * self.mt_density

            self.mtubes_x, self.mtubes_y = cells.single_cell_div_free(mtubes_xo, mtubes_yo)

            # initial angle of microtubules:
            self.mt_theta = np.arctan2(self.mtubes_y, self.mtubes_x)

            self.uxmt, self.uymt = self.mtubes_to_cell(cells, p)



