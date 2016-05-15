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
from scipy.ndimage.filters import gaussian_filter
from betse.science import finitediff as fd


def get_current(sim, cells, p):

    # calculate current density across gap junctions in x direction:
    J_gj_x_o = np.zeros(len(cells.mem_i))

    for flux_array, zi in zip(sim.fluxes_gj_x, sim.zs):
        J_i_x = flux_array * zi * p.F

        J_gj_x_o = J_gj_x_o + J_i_x

    # calculate current density across gap junctions in x direction:
    J_gj_y_o = np.zeros(len(cells.mem_i))

    for flux_array, zi in zip(sim.fluxes_gj_y, sim.zs):
        J_i_y = flux_array * zi * p.F

        J_gj_y_o = J_gj_y_o + J_i_y

    # need to calculate a divergence-free flow field through gap junctions, assuming
    # bulk electroneutrality of the electrolyte:

    # First calculate rate of change of charge in environment:
    if len(sim.charge_cells_time) > 1:

        d_rho_cells = (sim.charge_cells_time[-1] - sim.charge_cells_time[-2]) / p.dt

    else:
        d_rho_cells = 0


    # get the normal component to each cell membrane:
    J_gj = J_gj_x_o * cells.mem_vects_flat[:, 2] + J_gj_y_o * cells.mem_vects_flat[:, 3]


    # calculate divergence as the sum of this vector x each surface area, divided by cell volume:
    div_J_gj_o = (np.dot(cells.M_sum_mems, J_gj * cells.mem_sa) / cells.cell_vol)

    # add the rate of charge change to the divergence:
    div_J_gj_o = div_J_gj_o + d_rho_cells

    # calculate the reaction voltage required to counter-balance the flow field:
    Phi_gj = np.dot(cells.lapGJ_P_inv, div_J_gj_o)

    # calculate its gradient:
    gradPhi_gj = (Phi_gj[cells.cell_nn_i[:, 1]] - Phi_gj[cells.cell_nn_i[:, 0]]) / (cells.nn_len)

    gPhi_x = gradPhi_gj * cells.cell_nn_tx
    gPhi_y = gradPhi_gj * cells.cell_nn_ty

    J_gj_x = J_gj_x_o - gPhi_x
    J_gj_y = J_gj_y_o - gPhi_y

    sim.J_gj_x = interp.griddata((cells.mem_mids_flat[:,0],cells.mem_mids_flat[:,1]),J_gj_x,(cells.X,cells.Y),
                                  method=p.interp_type,fill_value=0)

    sim.J_gj_x = np.multiply(sim.J_gj_x,cells.maskECM)

    sim.J_gj_y = interp.griddata((cells.mem_mids_flat[:,0],cells.mem_mids_flat[:,1]),J_gj_y,(cells.X,cells.Y),
                                  method=p.interp_type,fill_value=0)


    sim.J_gj_y = np.multiply(sim.J_gj_y,cells.maskECM)


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
        Phi = np.dot(cells.lapENV_P_inv, div_J_env_o.ravel())
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

        # boundary conditions reinforced:
        sim.J_env_x[:, 0] = 0
        # right
        sim.J_env_x[:, -1] = 0
        # top
        sim.J_env_x[-1, :] = 0
        # bottom
        sim.J_env_x[0, :] = 0

        # left
        sim.J_env_y[:, 0] = 0
        # right
        sim.J_env_y[:, -1] = 0
        # top
        sim.J_env_y[-1, :] = 0
        # bottom
        sim.J_env_y[0, :] = 0
