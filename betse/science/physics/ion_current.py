#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

import numpy as np
from betse.science import finitediff as fd
from scipy.ndimage.filters import gaussian_filter
from betse.science import sim_toolbox as stb

def get_current(sim, cells, p):


    # calculate membrane current density (- as fluxes were defined into cell)
    sim.Jmem = -np.dot(sim.zs * p.F, sim.fluxes_mem)

    # calculate current density across cell membranes via gap junctions:
    sim.Jgj = np.dot(sim.zs * p.F, sim.fluxes_gj)

    # add the free current sources together into a single transmembrane current:
    sim.Jn = sim.Jmem + sim.Jgj

    # multiply final result by membrane surface area to obtain current (direction into cell is +)
    sim.I_mem = -sim.Jn*cells.mem_sa

    # components of GJ current:
    Jnx = sim.Jn * cells.nn_tx
    Jny = sim.Jn * cells.nn_ty


    # average intracellular current to cell centres:
    sim.J_cell_x = np.dot(cells.M_sum_mems, Jnx) / cells.num_mems
    sim.J_cell_y = np.dot(cells.M_sum_mems, Jny) / cells.num_mems



    # Current in the environment --------------------------------------------------------------------------------------
    if p.sim_ECM is True:

        # current densities in the environment:
        J_env_x_o = np.dot(p.F*sim.zs, sim.fluxes_env_x)
        J_env_y_o = np.dot(p.F*sim.zs, sim.fluxes_env_y)

        # reshape the matrix:
        J_env_x_o = J_env_x_o.reshape(cells.X.shape)
        J_env_y_o = J_env_y_o.reshape(cells.X.shape)

        if p.cell_polarizability == 0.0:

        # ---METHOD 1: Calculate potential to make currents divergence-free ---------------------

        # smooth the currents:
        # if p.smooth_level > 0.0:
        #     J_env_x_o = gaussian_filter(J_env_x_o, p.smooth_level)
        #     J_env_y_o = gaussian_filter(J_env_y_o, p.smooth_level)

            # conductivity in the media is modified by the environmental diffusion weight matrix:
            sigma = (1/p.media_rho)*sim.D_env_weight  # general conductivity

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

        else:

        # METHOD 2: Currents are already pre-corrected, but calculate sim.v_env-----------------------------
            # sim.v_env = np.sum(sim.Phi_vect, axis=0)

            Vbase = np.zeros(cells.X.shape)

            Vbase[:, 0] = sim.bound_V['L'] * (1 / cells.delta ** 2)
            Vbase[:, -1] = sim.bound_V['R'] * (1 / cells.delta ** 2)
            Vbase[0, :] = sim.bound_V['B'] * (1 / cells.delta ** 2)
            Vbase[-1, :] = sim.bound_V['T'] * (1 / cells.delta ** 2)

            Phi_base = np.dot(cells.lapENVinv, Vbase.ravel())


            sim.v_env = -np.dot(sim.zs, sim.Phi_vect) + Phi_base

            # sim.v_env = -np.dot(sim.zs, sim.Phi_vect)

            # corr = (cells.true_ecm_vol.min() / cells.ecm_vol)
            # sim.v_env = np.dot(cells.lapENVinv, -(corr * sim.rho_env.ravel()) / (p.eo * p.er))

            # calculate the gradient of Phi:
            gPhix, gPhiy = fd.gradient(sim.v_env.reshape(cells.X.shape), cells.delta)


            sim.E_env_x = -gPhix
            sim.E_env_y = -gPhiy

        #-------------------------------------------------------------------------------------------------------

        #Helmholtz-Hodge decomposition to obtain divergence-free projection of currents (zero n_hat at boundary):
        _, sim.J_env_x, sim.J_env_y, _, _, _ = stb.HH_Decomp(J_env_x_o,
                                                             J_env_y_o, cells)

        # sim.Jxe = J_env_x_o
        # sim.Jye = J_env_y_o

        # sim.J_env_x = J_env_x_o
        # sim.J_env_y = J_env_y_o















