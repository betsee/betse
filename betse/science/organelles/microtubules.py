#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
**Microtubules** (i.e., vectors aligned in the direction of endogenous current
flux through cells) functionality.
'''

# ....................{ IMPORTS                           }....................
import numpy as np
from betse.lib.numpy import nparray
from betse.science.math import modulate
from betse.util.io.log import logs

# ....................{ CLASSES                           }....................
class Mtubes(object):
    '''
    High-level **microtubules modeller** (i.e., object encapsulating all
    cellular microtubules for the cell cluster) simulated at the current
    simulation time step.

    Attributes (Scalars)
    ----------
    Dr : float
        Microtubule rotational diffusion constant.
    charge_mtube : float
        Charge in Coulombs (C) on each microtubule.
    drag_r : float
        Radial drag force on the microtubule from the Stokes-Einstein relation
        and "Rotational Diffusion" theory.
    p_ind : float
        Inducible polarizability of tubulin in units of [C m^2/V].
    p_mtube : float
        Microtubule dipole moment in units of [C m]. This moment is assumed to
        be the sum of tubulin dimer dipole moments, which are 1740 D each.
    rmt : float
        Radius in meters (m) of each microtubule.
    tubulin_N : float
        Total number of tubulin molecules per microtubule.

    Attributes (Arrays)
    ----------
    mt_density : {ndarray, float}
        Inhomogeneous density. If :attr:`mt_space_density` is ``None``, this is
        a Numpy array; else, this is a float.

    Attributes (Arrays: Vector)
    ----------
    mtubes_x : ndarray
        One-dimensional Numpy array indexing each cell membrane such that each
        item is the normalized X component of the microtubule unit vector
        spatially situated at the midpoint of that membrane for this time step.
    mtubes_y : ndarray
        One-dimensional Numpy array indexing each cell membrane such that each
        item is the normalized Y component of the microtubule unit vector
        spatially situated at the midpoint of that membrane for this time step.
    mtubes_xo : ndarray
        One-dimensional Numpy array indexing each cell membrane such that each
        item is the non-normalized X component of the microtubule vector
        spatially situated at the midpoint of that membrane for this time step.
    mtubes_yo : ndarray
        One-dimensional Numpy array indexing each cell membrane such that each
        item is the non-normalized Y component of the microtubule vector
        spatially situated at the midpoint of that membrane for this time step.
    '''

    # ..................{ INITIALIZORS                      }..................
    def __init__(self, sim, cells, p) -> None:

        # Initialize all core scalar parameters pertaining to microtubules.
        self._init_scalars(cells, p)

        # initial angle of microtubules:
        self.mt_theta = np.arctan2(
            cells.mem_vects_flat[:,3], cells.mem_vects_flat[:,2])

        # normalized microtubule vectors from the cell centre point:
        self.mtubes_x = cells.mem_vects_flat[:,2]*self.mt_density
        self.mtubes_y = cells.mem_vects_flat[:,3]*self.mt_density

        # microtubule density function initialized:
        mtdx = np.dot(cells.M_sum_mems, self.mtubes_x*cells.mem_sa) / cells.cell_sa
        mtdy = np.dot(cells.M_sum_mems, self.mtubes_y*cells.mem_sa) / cells.cell_sa

        self.mtdf = ((mtdx[cells.mem_to_cells]*cells.mem_vects_flat[:,2] +
                                         mtdy[cells.mem_to_cells]*cells.mem_vects_flat[:,3]))

        # Initialize the modulator as a chemical factor altering microtubule
        # dynamics.
        self.modulator = np.ones(sim.mdl)


    def _init_scalars(self, cells, p) -> None:
        '''
        Initialize all core scalar parameters pertaining to microtubules.
        '''

        # Radius of a microtubule.
        self.rmt = p.mt_radius

        self.L = self.rmt * (cells.R_rads/cells.R_rads.mean())

        # Total number of tubulin molecules per microtubule.
        self.tubulin_N = (222/1.0e-6)*self.L

        # Charge on the microtubule in Coulombs.
        self.charge_mtube = (p.length_charge*p.q*self.L)/1.0e-6

        # Microtubule dipole moment [C m] (assumed to be sum of tubulin dimer
        # dipole moments, which are 1740 D each).
        self.p_mtube = self.tubulin_N*p.tubulin_dipole*3.33e-30

        # Inducible polarizability of tubulin C m2/V.
        self.p_ind = self.tubulin_N*p.tubulin_polar*1.0e-40

        # Calculate microtubule drag co-efficients.
        self._init_scalars_drag(cells)

        # Radial drag force on the microtubule from Stokes-Einstein relation
        # "Rotational Diffusion" theory.
        self.drag_r = self.C_rad*p.cytoplasm_viscocity*(self.L**3)

        # Microtubule rotational diffusion constant.
        self.Dr = ((p.kb*p.T)/self.drag_r)

        # Inhomogeneous density.
        if p.mt_space_density is not None:
            self.mt_density, _ = modulate.gradient_bitmap(
                cells.mem_i, cells, p, bitmap_filename=p.mt_space_density)
        else:
            self.mt_density = 1.0


    def _init_scalars_drag(self, cells) -> None:
        '''
        Calculate microtubule drag co-efficients from the drag force on
        cylinders given by Broersma, et al. and Hunt, et al. (1994).
        '''

        v = 1 / (np.log(self.L/self.rmt))

        g_para = -0.114 - 0.15 * v - 13.5 * (v ** 2) + 37 * (v ** 3) - 22 * (v ** 4)
        g_perp = 0.866 - 0.15 * v - 8.1 * (v ** 2) + 18 * (v ** 3) - 9 * (v ** 4)
        g_rad = -0.446 - 0.2 * v - 16 * (v ** 2) + 63 * (v ** 3) - 62 * (v ** 4)

        self.C_para = (2 * np.pi) / (np.log(self.L / (2 * self.rmt)) + g_para)
        self.C_perp = (4 * np.pi) / (np.log(self.L / (2 * self.rmt)) + g_perp)
        self.C_rad = ((1 / 3) * np.pi) / (np.log(self.L / (2 * self.rmt)) + g_rad)


    def reinit(self, cells, p) -> None:
        '''
        Reinitialize key microtubule parameters that may have changed in the
        simulation configuration file since the prior initialization (if any).
        '''

        # Initialize all core scalar parameters pertaining to microtubules.
        self._init_scalars(cells, p)

        # smoothing weights for membrane and central values:
        nfrac = p.smooth_cells
        self.smooth_weight_mem = ((nfrac*cells.num_mems[cells.mem_to_cells] -1)/(nfrac*cells.num_mems[cells.mem_to_cells]))
        self.smooth_weight_o = 1/(nfrac*cells.num_mems[cells.mem_to_cells])

        # amt = self.L/(2*self.rmt)
        # self.drag_r = (np.pi*p.cytoplasm_viscocity*(self.L**3))/(np.log(amt - 0.662 + (0.917/amt) - (0.05/amt**2)))

        if p.mtube_init_x is not None:
            # Initialize microtubules with mini-simulation to define initial
            # state.
            self.Phi_orient_x, _ = modulate.gradient_bitmap(
                cells.mem_i, cells, p, bitmap_filename=p.mtube_init_x)
        else:
            # self.Phi_orient_x = np.zeros(len(cells.mem_i))
            self.Phi_orient_x = cells.mem_vects_flat[:, 2]

        if p.mtube_init_y is not None:
            self.Phi_orient_y, _ = modulate.gradient_bitmap(
                cells.mem_i, cells, p, bitmap_filename=p.mtube_init_y)
        else:
            self.Phi_orient_y = cells.mem_vects_flat[:, 3]
            # self.Phi_orient_y = np.zeros(len(cells.mem_i))
        # if self.Phi_orient_x != cells.mem_vects_flat[:, 2] or self.Phi_orient_y != cells.mem_vects_flat[:,3]:

        # if there's a user-requested gradient for initial condition alignment, then perform the pre-simulation:
        self.presim(self.Phi_orient_x, self.Phi_orient_y, cells, p)

        self.uxmt, self.uymt = self.mtubes_to_cell(cells, p)

    # ..................{ UPDATERS                          }..................
    def update_mtubes(self, cells, sim, p) -> None:

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


    def mtubes_to_cell(self, cells, p) -> None:

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


    def remove_mtubes(
        self, target_inds_mem, target_inds_cell, cells, sim, p) -> None:
        '''
        Remove all microtubules from the cells and cell membranes with the
        passed indices.
        '''

        self.mtubes_x     = np.delete(self.mtubes_x, target_inds_mem)
        self.mtubes_y     = np.delete(self.mtubes_y, target_inds_mem)
        self.mtdf         = np.delete(self.mtdf, target_inds_mem)
        self.mt_theta     = np.delete(self.mt_theta, target_inds_mem)
        self.modulator    = np.delete(self.modulator, target_inds_mem)
        self.L            = np.delete(self.L, target_inds_mem)
        self.charge_mtube = np.delete(self.charge_mtube, target_inds_mem)
        self.tubulin_N    = np.delete(self.tubulin_N, target_inds_mem)
        self.p_ind        = np.delete(self.p_ind, target_inds_mem)
        self.drag_r       = np.delete(self.drag_r, target_inds_mem)
        self.C_perp       = np.delete(self.C_perp, target_inds_mem)
        self.Dr           = np.delete(self.Dr, target_inds_mem)

        # If the "mt_density" instance variable has been conditionally defined
        # by the _init_scalars() method to be a Numpy array rather than a
        # float, remove the required items from this array.
        if nparray.is_array(self.mt_density):
            self.mt_density = np.delete(self.mt_density, target_inds_mem)

        # mtux2 = np.delete(self.uxmt, target_inds_cell)
        # self.uxmt = mtux2*1
        #
        # mtuy2 = np.delete(self.uymt, target_inds_cell)
        # self.uymt = mtuy2*1

        # mtuu2 = np.delete(self.umtn, target_inds_cell)
        # self.umtn = mtuu2*1


    def presim(self, gFxo, gFyo, cells, p) -> None:
        '''
        Initialize microtubule orientation using gradient of a vector field
        '''

        if len(gFxo) == len(cells.cell_i):
            gFxo = gFxo[cells.mem_to_cells]
            gFyo = gFyo[cells.mem_to_cells]

        if p.mtube_init_x is not None and p.mtube_init_y is not None:
            logs.log_info('-------------------------------')
            logs.log_info(
                'Preinitializing microtubule x- and y- coordinates with '
                '%f and %f', p.mtube_init_x,  p.mtube_init_y)
            logs.log_info('-------------------------------')

            # rotate the axis of the model:
            rotangle = (p.mtube_init_rotangle*np.pi)/180.0

            gFx = gFxo*np.cos(rotangle) - gFyo*np.sin(rotangle)
            gFy = gFxo*np.sin(rotangle) + gFyo*np.cos(rotangle)

            # magnitude of the orienting field:
            magF = (np.sqrt(gFx ** 2 + gFy ** 2)) + 1.0e-15

            # set the microtubule vectors with the field values:
            mtubes_xo = (gFx/magF) * self.mt_density
            mtubes_yo = (gFy/magF) * self.mt_density

            self.mtubes_x, self.mtubes_y = cells.single_cell_div_free(mtubes_xo, mtubes_yo)

            # initial angle of microtubules:
            self.mt_theta = np.arctan2(self.mtubes_y, self.mtubes_x)

            self.uxmt, self.uymt = self.mtubes_to_cell(cells, p)
