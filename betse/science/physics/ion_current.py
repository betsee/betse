#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

import numpy as np
from scipy.ndimage.filters import gaussian_filter

from betse.science import sim_toolbox as stb
from betse.science.math import finitediff as fd


def get_current(sim, cells, p):


    # nmx = cells.mem_vects_flat[:,2]
    # nmy = cells.mem_vects_flat[:,3]

    # calculate membrane current density (- as fluxes were defined into cell)
    sim.Jmem = -np.dot(sim.zs * p.F, sim.fluxes_mem)

    # calculate current density across cell membranes via gap junctions:
    sim.Jgj = np.dot(sim.zs * p.F, sim.fluxes_gj)

    # add the free current sources together into a single transmembrane current:
    sim.Jn = sim.Jmem + sim.Jgj + sim.extra_J_mem

    # multiply final result by membrane surface area to obtain current (direction into cell is +)
    sim.I_mem = -sim.Jn*cells.mem_sa

    # components of membrane current:
    Jnx = sim.Jn * cells.nn_tx
    Jny = sim.Jn * cells.nn_ty

    # average intracellular current to cell centres:
    # sim.J_cell_x = np.dot(cells.M_sum_mems, Jnx) / cells.num_mems
    # sim.J_cell_y = np.dot(cells.M_sum_mems, Jny) / cells.num_mems

    sim.J_cell_x = np.dot(cells.M_sum_mems, Jnx*cells.mem_sa) / cells.cell_sa
    sim.J_cell_y = np.dot(cells.M_sum_mems, Jny*cells.mem_sa) / cells.cell_sa


    # Current in the environment --------------------------------------------------------------------------------------
    if p.sim_ECM is True:

        # diffusive component of current densities in the environment:
        J_env_x_o = np.dot(p.F*sim.zs, sim.fluxes_env_x)
        J_env_y_o = np.dot(p.F*sim.zs, sim.fluxes_env_y)

        # calculate div J mems mapped to ecm:
        J_mems_env = np.zeros(sim.edl)
        J_mems_env[cells.map_mem2ecm] = -sim.Jmem

        # div_cells = (J_mems_env * (cells.memSa_per_envSquare / cells.ecm_vol)).reshape(cells.X.shape)

        # reshape the matrix:
        J_env_x_o = J_env_x_o.reshape(cells.X.shape)
        J_env_y_o = J_env_y_o.reshape(cells.X.shape)

        # conductivity in the media is modified by the environmental diffusion weight matrix:
        sigma = np.dot((((sim.zs ** 2) * p.q * p.F) / (p.kb * p.T)), sim.cc_env*sim.D_env).reshape(cells.X.shape)

        #---Calculate divergences for concentration & transmembrane fluxes ---------------------------------------------
        div_Jo = fd.divergence(J_env_x_o/sigma, J_env_y_o/sigma, cells.delta, cells.delta)

        div_Jo[:,0] = sim.bound_V['L']*(1/cells.delta**2)
        div_Jo[:,-1] = sim.bound_V['R']*(1/cells.delta**2)
        div_Jo[0,:] = sim.bound_V['B']*(1/cells.delta**2)
        div_Jo[-1,:] = sim.bound_V['T']*(1/cells.delta**2)

        # Calculate a voltage that resists the divergence:
        Phi = np.dot(cells.lapENVinv, (-div_Jo).ravel())

        if p.smooth_level > 0.0:
            # smoothing of Phi:
            Phi = gaussian_filter(Phi.reshape(cells.X.shape), p.smooth_level, mode='constant')

        sim.v_env = Phi.ravel()

        #--------------------------------------------------------------------------------------------------------------
        # calculate the gradient of v_env:
        gPhix, gPhiy = fd.gradient(Phi.reshape(cells.X.shape), cells.delta)

        sim.E_env_x = -gPhix
        sim.E_env_y = -gPhiy

        #Helmholtz-Hodge decomposition to obtain divergence-free projection of actual currents (zero n_hat at boundary):
        _, sim.J_env_x, sim.J_env_y, _, _, _ = stb.HH_Decomp(J_env_x_o, J_env_y_o, cells)



















