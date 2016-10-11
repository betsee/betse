#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

import numpy as np
from betse.science import finitediff as fd
from scipy.ndimage.filters import gaussian_filter

def get_current(sim, cells, p):


    # # Basic method with no correction to currents -----------------------------------------------------------------
    #
    # # calculate current across cell membranes:
    # # (reverse polarity to account for direction of cell membrane normals)
    # Jn = -np.dot(sim.zs*p.F, sim.fluxes_mem)
    #
    # # calculate current across cell membranes via gap junctions:
    # Jgj = np.dot(sim.zs*p.F, sim.fluxes_gj)
    #
    # # for now, add the two current sources together into a single transmembrane current:
    # sim.Jn = Jn + Jgj
    #
    # # get components of corrected current at membrane mids:
    # Jnx = sim.Jn * cells.mem_vects_flat[:, 2]
    # Jny = sim.Jn * cells.mem_vects_flat[:, 3]
    #
    # # map current to the cell centers:
    # sim.J_cell_x = np.dot(cells.M_sum_mems, Jnx) / cells.num_mems
    # sim.J_cell_y = np.dot(cells.M_sum_mems, Jny) / cells.num_mems
    #
    # # multiply final result by membrane surface area to obtain current (negative assigns into cell + for plotting)
    # sim.I_mem = -sim.Jn*cells.mem_sa

    # Whole System Method: --------------------------------------------------------------------------------------

    # calculate current across cell membranes:
    # (reverse polarity to account for direction of cell membrane normals)
    Jn = -np.dot(sim.zs*p.F, sim.fluxes_mem)

    # calculate current across cell membranes via gap junctions:
    Jgj = np.dot(sim.zs*p.F, sim.fluxes_gj)

    # for now, add the two current sources together into a single transmembrane current:
    sim.Jn = Jn + Jgj

    # calculate divergence as the sum of this vector x each surface area, divided by cell volume:
    div_J = (np.dot(cells.M_sum_mems, sim.Jn * cells.mem_sa) / cells.cell_vol)

    Phi = np.dot(cells.lapGJ_P_inv, div_J)

    gPhi = (Phi[cells.cell_nn_i[:, 1]] - Phi[cells.cell_nn_i[:, 0]]) / (cells.nn_len)

    Jn_corr = sim.Jn - gPhi

    sim.gPhi = gPhi

    Jnx = Jn_corr*cells.mem_vects_flat[:, 2]
    Jny = Jn_corr*cells.mem_vects_flat[:, 3]

    # map current to the cell centers:
    sim.J_cell_x = np.dot(cells.M_sum_mems, Jnx) / cells.num_mems
    sim.J_cell_y = np.dot(cells.M_sum_mems, Jny) / cells.num_mems

    # multiply final result by membrane surface area to obtain current (negative assigns into cell + for plotting)
    sim.I_mem = -sim.Jn*cells.mem_sa


    # Current in the environment --------------------------------------------------------------------------------------
    if p.sim_ECM is True:

        # non divergence free current densities in the environment:

        # J_env_x_o = np.dot(sim.zs, sim.fluxes_env_x)
        # J_env_y_o = np.dot(sim.zs, sim.fluxes_env_y)

        J_env_x_o = np.zeros(sim.edl)
        J_env_y_o = np.zeros(sim.edl)

        J_env_x_o[cells.map_cell2ecm] = sim.J_cell_x
        J_env_y_o[cells.map_cell2ecm] = sim.J_cell_y

        if p.smooth_level > 0.0:
            J_env_x_o = gaussian_filter(J_env_x_o.reshape(cells.X.shape), p.smooth_level).ravel()
            J_env_y_o = gaussian_filter(J_env_y_o.reshape(cells.X.shape), p.smooth_level).ravel()


        # determine corrected current densities assuming
        # bulk electrolyte neutrality of entire system:

        # # Next, calculate the divergence of the environmental current density:
        div_J_env_o = fd.divergence(J_env_x_o.reshape(cells.X.shape), J_env_y_o.reshape(cells.X.shape),
            cells.delta, cells.delta)

        # Find the value of the environmental electric potential:
        Phi = np.dot(cells.lapENV_P_inv, div_J_env_o.ravel())
        Phi = Phi.reshape(cells.X.shape)

        sim.Phi_env = Phi

        # Take the grid gradient of the scaled internal potential:
        gPhix, gPhiy = fd.gradient(Phi, cells.delta)

        # subtract the potential term from the solution to yield the actual current density in the environment:
        sim.J_env_x = J_env_x_o.reshape(cells.X.shape) - gPhix
        sim.J_env_y = J_env_y_o.reshape(cells.X.shape) - gPhiy

        # if p.smooth_level > 0.0:
        #     sim.J_env_x = gaussian_filter(sim.J_env_x, p.smooth_level)
        #     sim.J_env_y = gaussian_filter(sim.J_env_y, p.smooth_level)

        # sim.E_env_x = sim.J_env_x*(sim.D_env_weight)
        # sim.E_env_y = sim.J_env_y*(sim.D_env_weight)





