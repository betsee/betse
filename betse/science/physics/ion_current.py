#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

import numpy as np
from betse.science import finitediff as fd
from scipy.ndimage.filters import gaussian_filter

def get_current(sim, cells, p):

    # calculate current across cell membranes:----------------------------
    sim.Jn = np.dot(sim.zs*p.F, sim.fluxes_mem)

    sim.Jn = -sim.Jn  # reverse polarity to account for direction of cell membrane normals

    # calculate a rate of change of current at the boundary prior to modification of sim.Jn:

    d_rho = np.zeros(sim.cdl)

    if p.sim_ECM is True:

        b_charge = sim.d_rho + sim.J_TJ

    else:

        b_charge = sim.d_rho

    cell_I = np.dot(cells.M_sum_mems, b_charge*cells.mem_sa)/cells.cell_vol

    d_rho[cells.bflags_cells] = cell_I[cells.bflags_cells]

    # The current density needs to be corrected to account for overall charge neutrality of system:
    sim.Jn, sim.J_cell_x, sim.J_cell_y = cells.zero_div_cell(sim.Jn, rho = d_rho, bc = b_charge[cells.bflags_mems],
                                                             open_bounds=False)

    # sim.Jn, sim.J_cell_x, sim.J_cell_y = cells.zero_div_cell(sim.Jn, rho = d_rho, bc = sim.d_rho[cells.bflags_mems],
    #                                                          open_bounds=False)

    # multiply final result by membrane surface area to obtain current (negative assigns into cell + for plotting)
    sim.I_mem = -sim.Jn*cells.mem_sa

    # Current in the environment --------------------------------------------------------------------------------------
    if p.sim_ECM is True:

        # non divergence free current densities in the environment:
        # J_env_x_o = np.dot(sim.zs*p.F, sim.fluxes_env_x)
        # J_env_y_o = np.dot(sim.zs*p.F, sim.fluxes_env_y)

        J_env_x_o = np.zeros(sim.edl)
        J_env_y_o = np.zeros(sim.edl)

        # add in current density contribution from cell membranes:
        # J_env_x_o[cells.map_mem2ecm] = J_env_x_o[cells.map_mem2ecm] + sim.Jn*cells.mem_vects_flat[:,2]
        # J_env_y_o[cells.map_mem2ecm] = J_env_y_o[cells.map_mem2ecm] + sim.Jn*cells.mem_vects_flat[:,3]

        # J_env_x_o[cells.map_mem2ecm] = sim.Jn*cells.mem_vects_flat[:,2]
        # J_env_y_o[cells.map_mem2ecm] = sim.Jn*cells.mem_vects_flat[:,3]

        J_env_x_o[cells.map_cell2ecm] = sim.J_cell_x
        J_env_y_o[cells.map_cell2ecm] = sim.J_cell_y

        # optional smoothing of current using gaussian:
        # if p.smooth_level > 0.0:
        #     J_env_x_o = gaussian_filter(J_env_x_o.reshape(cells.X.shape), p.smooth_level)
        #     J_env_y_o = gaussian_filter(J_env_y_o.reshape(cells.X.shape), p.smooth_level)

        # determine corrected current densities assuming
        # bulk electrolyte neutrality:

        # # Next, calculate the divergence of the environmental current density:
        div_J_env_o = fd.divergence(J_env_x_o.reshape(cells.X.shape), J_env_y_o.reshape(cells.X.shape),
            cells.delta, cells.delta)

        # Find the value of the environmental electric potential:
        Phi = np.dot(cells.lapENV_P_inv, div_J_env_o.ravel())
        Phi = Phi.reshape(cells.X.shape)

        # Take the grid gradient of the scaled internal potential:
        gPhix, gPhiy = fd.gradient(Phi, cells.delta)

        # subtract the potential term from the solution to yield the actual current density in the environment:
        sim.J_env_x = J_env_x_o.reshape(cells.X.shape) - gPhix
        sim.J_env_y = J_env_y_o.reshape(cells.X.shape) - gPhiy


        sim.E_env_x = sim.J_env_x*p.media_sigma
        sim.E_env_y = sim.J_env_y*p.media_sigma




