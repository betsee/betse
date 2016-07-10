#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

import copy
import os
import os.path
import time
from random import shuffle

import numpy as np
from scipy import interpolate as interp
from betse.science import finitediff as fd
from scipy.sparse.linalg import lsmr


def get_current(sim, cells, p):

    # calculate current across cell membranes:----------------------------
    sim.I_mem = np.zeros(len(cells.mem_i))
    for flux_array, zi in zip(sim.fluxes_mem, sim.zs):
        I_i = flux_array * zi * p.F * cells.mem_sa

        sim.I_mem = sim.I_mem + I_i

    # divide final result by membrane surface area to obtain a component of current density
    # component is negative as transmembrane fluxes point into the cell, but mem normals point out:
    J_trans_mem = -sim.I_mem/cells.mem_sa

    # calculate current density across gap junctions:
    J_gj_o = np.zeros(len(cells.mem_i))

    for flux_array, zi in zip(sim.fluxes_gj, sim.zs):
        J_i = flux_array * zi * p.F

        J_gj_o = J_gj_o + J_i

    # raw current density through the gap junctions:
    sim.J_gj_x = J_gj_o * cells.mem_vects_flat[:,2]
    sim.J_gj_y = J_gj_o * cells.mem_vects_flat[:,3]

    # total current density across the membranes (uncorrected by Continuity Equation):
    J_mem_xo = sim.J_gj_x + J_trans_mem * cells.mem_vects_flat[:,2]
    J_mem_yo = sim.J_gj_y + J_trans_mem * cells.mem_vects_flat[:,3]

    # J_mem_xo = sim.J_gj_x
    # J_mem_yo = sim.J_gj_y

    # First calculate rate of change of charge in cell:
    if len(sim.charge_cells_time) > 1:

        d_rho_cells = (sim.charge_cells_time[-1] - sim.charge_cells_time[-2]) / p.dt

    else:
        d_rho_cells = np.zeros(len(cells.cell_i))

    # Next, calculate the divergence of cell current density:

    # get the component normal to the membrane:
    J_mem_n = J_mem_xo * cells.mem_vects_flat[:, 2] + J_mem_yo * cells.mem_vects_flat[:, 3]

    # calculate divergence as the sum of this vector x each surface area, divided by cell volume:
    div_Jmem = (np.dot(cells.M_sum_mems, J_mem_n * cells.mem_sa) / cells.cell_vol)

    # calculate the reaction potential required to counter-balance the flow field:

    # if p.deformation is True:
    #
    #     b = -div_Jmem - d_rho_cells
    #
    #     V_react = lsmr(cells.lapGJ, b)[0]
    #
    # else:
    V_react = np.dot(cells.lapGJinv, -div_Jmem - d_rho_cells)

    # calculate its gradient:
    gradV_react = (V_react[cells.cell_nn_i[:, 1]] - V_react[cells.cell_nn_i[:, 0]]) / (cells.nn_len)

    gV_x = gradV_react * cells.mem_vects_flat[:, 2]
    gV_y = gradV_react * cells.mem_vects_flat[:, 3]

    #correct the current density by the reaction potential:
    # sim.J_mem_x = J_mem_xo - gV_x
    # sim.J_mem_y = J_mem_yo - gV_y

    sim.J_mem_x = J_mem_xo
    sim.J_mem_y = J_mem_yo

    # average the components at cell centres:
    sim.J_cell_x = np.dot(cells.M_sum_mems, sim.J_mem_x) / cells.num_mems
    sim.J_cell_y = np.dot(cells.M_sum_mems, sim.J_mem_y) / cells.num_mems

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

        # determine correction factor for current densities assuming
        # bulk electrolyte neutrality

        # First calculate rate of change of charge in environment:
        if len(sim.charge_env_time) > 1:

            d_rho_env = (sim.charge_env_time[-1] - sim.charge_env_time[-2])/p.dt

        else:
            d_rho_env = np.zeros(len(cells.xypts))

        # Next, calculate the divergence of the environmental current density:
        div_J_env_o = fd.divergence(J_env_x_o.reshape(cells.X.shape), J_env_y_o.reshape(cells.X.shape),
            cells.delta, cells.delta)

        # add the rate of charge change to the divergence:
        div_J_env_o = div_J_env_o + d_rho_env.reshape(cells.X.shape)

        # Find the value of the correcting potential field Phi:
        Phi = np.dot(cells.lapENVinv, div_J_env_o.ravel())
        Phi = Phi.reshape(cells.X.shape)

        # enforce zero normal gradient boundary conditions on P:
        Phi[:, 0] = Phi[:, 1]
        Phi[:, -1] = Phi[:, -2]
        Phi[0, :] = Phi[1, :]
        Phi[-1, :] = Phi[-2, :]

        # Take the grid gradient of the scaled internal potential:
        gPhix, gPhiy = fd.gradient(Phi, cells.delta)

        # subtract the potential term from the solution to yield a divergence-free flow field:
        sim.J_env_x = J_env_x_o.reshape(cells.X.shape) - gPhix
        sim.J_env_y = J_env_y_o.reshape(cells.X.shape) - gPhiy

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
