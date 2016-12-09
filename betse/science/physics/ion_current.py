#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

import numpy as np
from betse.science import finitediff as fd
from scipy.ndimage.filters import gaussian_filter
from betse.science import sim_toolbox as stb

def get_current(sim, cells, p):


    # membrane normal vectors short form:
    nx = cells.mem_vects_flat[:, 2]
    ny = cells.mem_vects_flat[:, 3]

    # calculate membrane current density (- as fluxes were defined into cell)
    Jn = -np.dot(sim.zs * p.F, sim.fluxes_mem)

    # calculate current density across cell membranes via gap junctions:
    Jgj = np.dot(sim.zs * p.F, sim.fluxes_gj)

    # add the free current sources together into a single transmembrane current:
    sim.Jn = Jn + Jgj

    # multiply final result by membrane surface area to obtain current (direction into cell is +)
    sim.I_mem = -sim.Jn*cells.mem_sa

    # components of current:
    Jnx = sim.Jn * nx
    Jny = sim.Jn * ny

    # average to cell centres:
    sim.J_cell_x = np.dot(cells.M_sum_mems, Jnx) / cells.num_mems
    sim.J_cell_y = np.dot(cells.M_sum_mems, Jny) / cells.num_mems


    # Current in the environment --------------------------------------------------------------------------------------
    if p.sim_ECM is True:

        # current densities in the environment:
        J_env_x_o = np.dot(p.F*sim.zs, sim.fluxes_env_x)
        J_env_y_o = np.dot(p.F*sim.zs, sim.fluxes_env_y)

        # # base current in environment is equivalent to transmembrane current of cells:
        # J_env_x_o[cells.map_cell2ecm] = J_env_x_o[cells.map_cell2ecm] + sim.J_cell_x*(
        #     cells.memSa_per_envSquare[cells.map_cell2ecm]/cells.delta**2)
        #
        # J_env_y_o[cells.map_cell2ecm] = J_env_y_o[cells.map_cell2ecm] + sim.J_cell_y*(
        #     cells.memSa_per_envSquare[cells.map_cell2ecm]/cells.delta**2)

        # reshape the matrix:
        J_env_x_o = J_env_x_o.reshape(cells.X.shape)
        J_env_y_o = J_env_y_o.reshape(cells.X.shape)



        #---METHOD 1: Calculate potential to make currents divergence-free ---------------------

        # smooth the currents:
        # if p.smooth_level > 0.0:
        #     J_env_x_o = gaussian_filter(J_env_x_o, p.smooth_level)
        #     J_env_y_o = gaussian_filter(J_env_y_o, p.smooth_level)

        # conductivity in the media is modified by the environmental diffusion weight matrix:
        sigma = (1/p.media_rho)*sim.D_env_weight   # general conductivity

        div_Jo = fd.divergence((J_env_x_o / sigma), (J_env_y_o / sigma), cells.delta, cells.delta)

        # add-in any boundary conditions pertaining to an applied (i.e. external) voltage:
        div_Jo[:,0] = sim.bound_V['L']*(1/cells.delta**2)
        div_Jo[:,-1] = sim.bound_V['R']*(1/cells.delta**2)
        div_Jo[0,:] = sim.bound_V['B']*(1/cells.delta**2)
        div_Jo[-1,:] = sim.bound_V['T']*(1/cells.delta**2)

        # calculate the voltage balancing the divergence of the currents:
        Phi = np.dot(cells.lapENVinv, div_Jo.ravel())

        # the global environmental voltage is equal to Phi:
        sim.v_env = Phi

        # calculate the gradient of Phi:
        gPhix, gPhiy = fd.gradient(Phi.reshape(cells.X.shape), cells.delta)

        sim.E_env_x = -gPhix
        sim.E_env_y = -gPhiy

        # # correct the currents using Phi:
        # sim.J_env_x = J_env_x_o + sim.E_env_x*sigma
        # sim.J_env_y = J_env_y_o + sim.E_env_y*sigma

        #Helmholtz-Hodge decomposition to obtain divergence-free projection of currents (zero at boundary):
        _, sim.J_env_x, sim.J_env_y, _, _, _ = stb.HH_Decomp(J_env_x_o,
                                                             J_env_y_o, cells)

        # sim.J_env_x = sim.E_env_x*(1/p.media_rho)
        # sim.J_env_y = sim.E_env_y*(1/p.media_rho)


        #
        # # smooth the currents:
        # if p.smooth_level > 0.0:
        #     sim.J_env_x = gaussian_filter(sim.J_env_x, p.smooth_level)
        #     sim.J_env_y = gaussian_filter(sim.J_env_y, p.smooth_level)



        #--METHOD 2: reconstruct field from currents using Helmholtz-Hodge decomposition-----
        # bounds = {'L': -sim.bound_V['L'], 'R': -sim.bound_V['R'], 'T': -sim.bound_V['T'], 'B': -sim.bound_V['B'] }
        # AA, Ax, Ay, BB, Bx, By = stb.HH_Decomp(J_env_x_o, J_env_y_o, cells, bounds = bounds)
        #
        # # the environmental voltage is the negative of the curl-free field:
        # sim.v_env = -(1/p.media_rho)*BB.ravel()
        # sim.E_env_x = Bx
        # sim.E_env_y = By      #
        # # sim.Bz = AA           # magnetic (H) field
        #
        # # fluxes are already individually divergence-corrected; save final currents:
        # sim.J_env_x = J_env_x_o
        # sim.J_env_y = J_env_y_o












