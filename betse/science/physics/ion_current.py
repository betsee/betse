#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
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
    Jnx = sim.Jgj * cells.nn_tx
    Jny = sim.Jgj * cells.nn_ty

    # average intracellular current to cell centres:
    sim.J_cell_x = np.dot(cells.M_sum_mems, Jnx) / cells.num_mems
    sim.J_cell_y = np.dot(cells.M_sum_mems, Jny) / cells.num_mems

    # # calculate field in the cells resulting from intracellular current:
    # sigma = np.dot((((sim.zs**2)*p.q*p.F*sim.D_free)/(p.kb*p.T)), sim.cc_cells)*p.tissue_rho
    #
    # divJc = np.dot(cells.M_sum_mems, (sim.Jgj/sigma[cells.mem_to_cells])*cells.mem_sa)/cells.cell_vol
    #
    # sim.v_cell = np.dot(cells.lapGJ_P_inv, -divJc)


    # Current in the environment --------------------------------------------------------------------------------------
    if p.sim_ECM is True:

        # diffusive component of current densities in the environment:
        J_env_x_o = np.dot(p.F*sim.zs, sim.fluxes_env_x)
        J_env_y_o = np.dot(p.F*sim.zs, sim.fluxes_env_y)

        # reshape the matrix:
        J_env_x_o = J_env_x_o.reshape(cells.X.shape)
        J_env_y_o = J_env_y_o.reshape(cells.X.shape)

        # conductivity in the media is modified by the environmental diffusion weight matrix:
        sigma = np.dot((((sim.zs ** 2) * p.q * p.F) / (p.kb * p.T)), sim.cc_env*sim.D_env).reshape(cells.X.shape)

        #---Calculate divergences for concentration & transmembrane fluxes ---------------------------------------------
        # div_Jo = -fd.divergence(J_env_x_o/sigma, J_env_y_o/sigma, cells.delta, cells.delta)

        div_Jo = -fd.divergence(sim.conc_J_x.reshape(cells.X.shape)/sigma,
                                sim.conc_J_y.reshape(cells.X.shape)/sigma, cells.delta, cells.delta)

        # div_Jo = np.zeros(cells.X.shape)
        #
        # # determine finite divergence from cellular trans-membrane fluxes to the environmental space:
        # div_from_cells = (np.dot(cells.M_sum_mems,
        #                  (sim.Jmem/(sigma.ravel()[cells.map_mem2ecm])*cells.mem_sa))
        #                   /cells.cell_vol)*cells.cell2env_corrF
        #
        # div_from_cells_map = np.zeros(sim.edl)
        # div_from_cells_map[cells.map_cell2ecm] = div_from_cells
        # # div_from_cells_map = gaussian_filter(div_from_cells_map.reshape(cells.X.shape), p.smooth_level, mode = 'constant')
        # div_from_cells_map = div_from_cells_map.reshape(cells.X.shape)
        # div_Jo = div_Jo + div_from_cells_map

        div_Jo[:,0] = sim.bound_V['L']*(1/cells.delta**2)
        div_Jo[:,-1] = sim.bound_V['R']*(1/cells.delta**2)
        div_Jo[0,:] = sim.bound_V['B']*(1/cells.delta**2)
        div_Jo[-1,:] = sim.bound_V['T']*(1/cells.delta**2)

        # Calculate a voltage that resists the divergence:
        Phi = np.dot(cells.lapENVinv, (div_Jo).ravel())

        #Helmholtz-Hodge decomposition to obtain divergence-free projection of currents (zero n_hat at boundary):
        _, sim.J_env_x, sim.J_env_y, _, _, _ = stb.HH_Decomp(J_env_x_o,
                                                             J_env_y_o, cells)

        # sim.J_env_x = J_env_x_o
        # sim.J_env_y = J_env_y_o

        sim.v_env = Phi.reshape(cells.X.shape)

        if p.smooth_level > 0.0:

            sim.v_env = gaussian_filter(sim.v_env, p.smooth_level, mode = 'constant')

        sim.v_env = sim.v_env.ravel()

        #--------------------------------------------------------------------------------------------------------------

        # calculate the gradient of v_env:
        gPhix, gPhiy = fd.gradient(sim.v_env.reshape(cells.X.shape), cells.delta)

        sim.E_env_x = -gPhix
        sim.E_env_y = -gPhiy

        # sim.J_env_x = sim.E_env_x*sigma
        # sim.J_env_y = sim.E_env_y*sigma

        # sim.J_env_x = J_env_x_o - sim.conc_J_x.reshape(cells.X.shape)
        # sim.J_env_y = J_env_y_o - sim.conc_J_y.reshape(cells.X.shape)


















