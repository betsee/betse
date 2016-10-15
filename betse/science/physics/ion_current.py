#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

import numpy as np
from betse.science import finitediff as fd
from scipy.ndimage.filters import gaussian_filter

def get_current(sim, cells, p):


    # Whole System Method: --------------------------------------------------------------------------------------

    # calculate current density across cell membranes:
    # (reverse polarity to account for direction of cell membrane normals)
    Jn = -np.dot(sim.zs*p.F, sim.fluxes_mem)

    # calculate current density across cell membranes via gap junctions:
    Jgj = np.dot(sim.zs*p.F, sim.fluxes_gj)

    # add the current sources together into a single transmembrane current:
    Jn_o = Jn + Jgj

    drho = np.zeros(sim.cdl)

    # bsigma = np.zeros(sim.cdl)
    # bsigma[cells.bflags_cells] = p.D_tj

    # correct the current density using the continuity equation:

    # calculate divergence as the sum of this vector x each surface area, divided by cell volume:
    div_J = (np.dot(cells.M_sum_mems, Jn_o*cells.mem_sa) / cells.cell_vol)

    # this boundary condition allows us to accumulate a surface charge:
    # drho[cells.bflags_cells] = div_J[cells.bflags_cells]

    Phi = np.dot(cells.lapGJ_P_inv, div_J - drho)

    sim.Phi = Phi

    gPhi = (Phi[cells.cell_nn_i[:, 1]] - Phi[cells.cell_nn_i[:, 0]]) / (cells.nn_len)

    sim.gPhi = gPhi

    Jn_corr = Jn_o - gPhi

    Jnx = Jn_corr*cells.mem_vects_flat[:, 2]
    Jny = Jn_corr*cells.mem_vects_flat[:, 3]

    # map current to the cell centers:
    sim.J_cell_x = np.dot(cells.M_sum_mems, Jnx) / cells.num_mems
    sim.J_cell_y = np.dot(cells.M_sum_mems, Jny) / cells.num_mems

    # # reassign the normal current to the un_corrected component
    # sim.Jn = Jn_o*1

    # assign the normal current to the corrected component
    sim.Jn = Jn_corr*1

    # multiply final result by membrane surface area to obtain current (direction into cell is +)
    sim.I_mem = -sim.Jn*cells.mem_sa

    #------------------------------------------------------------------

    # # finally, obtain a weighting function describing individual membrane moments of current/field:
    # J_ave = np.dot(cells.M_sum_mems, sim.Jn) / cells.num_mems
    #
    # inds_zero = (J_ave == 0.0).nonzero()
    #
    # J_ave[inds_zero] = 1.0  # add some small number for the case that net current is zero
    #
    # sim.J_weight = sim.Jn / J_ave[cells.mem_to_cells]


    # Current in the environment --------------------------------------------------------------------------------------
    if p.sim_ECM is True:

        # non divergence free current densities in the environment:

        # J_env_x_o = np.dot(sim.zs, sim.fluxes_env_x)
        # J_env_y_o = np.dot(sim.zs, sim.fluxes_env_y)

        J_env_x_o = np.zeros(sim.edl)
        J_env_y_o = np.zeros(sim.edl)

        # due to electrolyte screening and fluxes, current in extracellular spaces opposite that of in cells:
        J_env_x_o[cells.map_cell2ecm] = sim.J_cell_x
        J_env_y_o[cells.map_cell2ecm] = sim.J_cell_y

        if p.smooth_level > 0.0:
            J_env_x_o = gaussian_filter(J_env_x_o.reshape(cells.X.shape), p.smooth_level).ravel()
            J_env_y_o = gaussian_filter(J_env_y_o.reshape(cells.X.shape), p.smooth_level).ravel()


        # determine corrected current densities assuming
        # bulk electrolyte neutrality of entire system:

        # # Next, calculate the divergence of the environmental current density:
        # div_J_env_o = fd.divergence(J_env_x_o.reshape(cells.X.shape), J_env_y_o.reshape(cells.X.shape),
        #     cells.delta, cells.delta)

        div_J_env_o = fd.divergence(J_env_x_o.reshape(cells.X.shape)/sim.D_env_weight,
                                    J_env_y_o.reshape(cells.X.shape)/sim.D_env_weight,
                                    cells.delta, cells.delta)

        # Find the value of the environmental electric potential:
        Phi = np.dot(cells.lapENV_P_inv, div_J_env_o.ravel())
        Phi = Phi.reshape(cells.X.shape)

        sim.Phi_env = Phi

        sim.v_env = np.copy(sim.Phi_env.ravel())
        sim.v_env[cells.inds_env] = 0.0
        sim.v_env[cells.ecm_bound_k] = 0.0

        # Take the grid gradient of the scaled internal potential:
        gPhix, gPhiy = fd.gradient(Phi, cells.delta)

        # subtract the potential term from the solution to yield the actual current density in the environment:
        sim.J_env_x = (J_env_x_o.reshape(cells.X.shape) - gPhix)
        sim.J_env_y = (J_env_y_o.reshape(cells.X.shape) - gPhiy)






