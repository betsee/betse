#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

import copy
import os
import os.path
import time
from random import shuffle

import numpy as np
from betse.science import finitediff as fd
from betse.science import sim_toolbox as stb

def get_current(sim, cells, p):

    # calculate current across cell membranes:----------------------------
    sim.Jn = np.dot(sim.zs*p.F, sim.fluxes_mem)

    sim.Jn = -sim.Jn  # reverse polarity to account for direction of cell membrane normals

    # calculate a rate of change of current at the boundary prior to modification of sim.Jn:

    d_rho = np.zeros(sim.cdl)

    cell_I = np.dot(cells.M_sum_mems, sim.d_rho*cells.mem_sa)/cells.cell_vol

    d_rho[cells.bflags_cells] = cell_I[cells.bflags_cells]

    # The current density needs to be corrected to account for overall charge neutrality of system:
    sim.Jn, sim.J_cell_x, sim.J_cell_y = cells.zero_div_cell(sim.Jn, rho = d_rho, bc = sim.d_rho[cells.bflags_mems],
                                                             open_bounds=False)

    # sim.Jn = sim.Jn + sim.d_rho

    # multiply final result by membrane surface area to obtain current (negative assigns into cell + for plotting)
    sim.I_mem = -sim.Jn*cells.mem_sa

    # Current in the environment --------------------------------------------------------------------------------------
    if p.sim_ECM is True:

        # non divergence free current densities in the environment:
        J_env_x_o = np.zeros(len(cells.xypts))
        J_env_y_o = np.zeros(len(cells.xypts))

        for flux_array, zi in zip(sim.fluxes_env_x, sim.zs):

            J_i = flux_array * zi * p.F

            J_env_x_o = J_env_x_o + J_i

        for flux_array, zi in zip(sim.fluxes_env_y, sim.zs):
            J_i = flux_array * zi * p.F

            J_env_y_o = J_env_y_o + J_i

        # add in current density contribution from cell membranes:
        J_env_x_o[cells.map_mem2ecm] = J_env_x_o[cells.map_mem2ecm] + J_mem_xo
        J_env_y_o[cells.map_mem2ecm] = J_env_y_o[cells.map_mem2ecm] + J_mem_yo

        # determine correction factor for current densities assuming
        # bulk electrolyte neutrality:

        # # Next, calculate the divergence of the environmental current density:
        div_J_env_o = fd.divergence(J_env_x_o.reshape(cells.X.shape), J_env_y_o.reshape(cells.X.shape),
            cells.delta, cells.delta)

        # add the rate of charge change to the divergence:
        # div_J_env_o = div_J_env_o + d_rho_env.reshape(cells.X.shape)

        # Find the value of the environmental electric potential:
        Phi = np.dot(cells.lapENVinv, -div_J_env_o.ravel())
        Phi = Phi.reshape(cells.X.shape)

        # Take the grid gradient of the scaled internal potential:
        gPhix, gPhiy = fd.gradient(Phi, cells.delta)

        # subtract the potential term from the solution to yield the actual current density in the environment:
        sim.J_env_x = J_env_x_o.reshape(cells.X.shape) - gPhix
        sim.J_env_y = J_env_y_o.reshape(cells.X.shape) - gPhiy

        # sim.J_env_x = J_env_x_o.reshape(cells.X.shape)
        # sim.J_env_y = J_env_y_o.reshape(cells.X.shape)

        # sim.E_env_x = -gvx + sim.J_env_x*(p.media_sigma)
        # sim.E_env_y = -gvy + sim.J_env_y*(p.media_sigma)


        # Currents are not zero at the boundary, but for plotting purposes set them to zero
        sim.J_env_x[:, 0] = 0
        sim.J_env_x[:, 1] = 0
        # right
        sim.J_env_x[:, -1] = 0
        sim.J_env_x[:, -2] = 0
        # top
        sim.J_env_x[-1, :] = 0
        sim.J_env_x[-2, :] = 0
        # bottom
        sim.J_env_x[0, :] = 0
        sim.J_env_x[1, :] = 0

        # left
        sim.J_env_y[:, 0] = 0
        sim.J_env_y[:, 1] = 0
        # right
        sim.J_env_y[:, -1] = 0
        sim.J_env_y[:, -2] = 0
        # top
        sim.J_env_y[-1, :] = 0
        sim.J_env_y[-2, :] = 0
        # bottom
        sim.J_env_y[0, :] = 0
        sim.J_env_y[1, :] = 0


