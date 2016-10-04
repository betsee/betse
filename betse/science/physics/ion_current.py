#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

import numpy as np
from betse.science import finitediff as fd
from scipy.ndimage.filters import gaussian_filter

def get_current(sim, cells, p):

    # Whole System Method: --------------------------------------------------------------------------------------

    # calculate current across cell membranes:
    # (reverse polarity to account for direction of cell membrane normals)
    Jn = -np.dot(sim.zs*p.F, sim.fluxes_mem)

    # calculate current across cell membranes via gap junctions:
    Jgj = np.dot(sim.zs*p.F, sim.fluxes_gj)

    # for now, add the two current sources together into a single transmembrane current:
    Jn = Jn + Jgj


    # calculate a rate of change of current at the boundary prior to modification of sim.Jn:
    # initialize an empty cell-centred storge vector
    d_rho = np.zeros(sim.cdl)

    # As it's very difficult to do this in terms of raw environmental fluxes, alter the all-important boundary
    # current in terms of the tight junction relative weight. When TJ_load is ~ 0.1, the system behaves like an
    # open boundary with low TEP. At TJ_load ~ 0.99, the system progresses towards a closed TJ boundary with high TEP
    if p.sim_ECM is True:

        # For sim_ECM, do this in a way that will be modified with cutting events and TJ interventions:
        TJ_load = sim.D_env_weight.ravel()[cells.map_cell2ecm][cells.bflags_cells]

        # TJ_load = 1 / (1 + TJ_load_o)

    else:
        # For no ECM, just use the default relative TJ diffusion constant:
        TJ_load = p.D_tj


    # # TJ_load = 1 / (1 + (p.D_tj/0.1))
    # TJ_load = 1.0

    # calculate divergence as the sum of this vector x each surface area, divided by cell volume:
    div_J = (np.dot(cells.M_sum_mems, Jn * cells.mem_sa) / cells.cell_vol)

    d_rho[cells.bflags_cells] = div_J[cells.bflags_cells]*TJ_load

    Phi = np.dot(cells.lapGJ_P_inv, div_J - d_rho)

    gPhi = (Phi[cells.cell_nn_i[:, 1]] - Phi[cells.cell_nn_i[:, 0]]) / (cells.nn_len)

    # make the field divergence-free:
    sim.Jn = Jn - gPhi*0.99

    # get components of corrected current at membrane mids:
    Jnx = sim.Jn * cells.mem_vects_flat[:, 2]
    Jny = sim.Jn * cells.mem_vects_flat[:, 3]

    # map current to the cell centers:
    sim.J_cell_x = np.dot(cells.M_sum_mems, Jnx) / cells.num_mems
    sim.J_cell_y = np.dot(cells.M_sum_mems, Jny) / cells.num_mems

    # multiply final result by membrane surface area to obtain current (negative assigns into cell + for plotting)
    sim.I_mem = -sim.Jn*cells.mem_sa


# alternative method --------------------------------------------------------------------------------------------------
#     # calculate current across cell membranes:
#     # (reverse polarity to account for direction of cell membrane normals)
#     Jn = -np.dot(sim.zs*p.F, sim.fluxes_mem)
#
#     # calculate current across cell membranes via gap junctions:
#     Jgj = np.dot(sim.zs*p.F, sim.fluxes_gj)
#
#     # # for now, add the two current sources together into a single transmembrane current:
#     # Jn = Jn + Jgj
#
#     # correct the gap junction current profile using a Phi potential defined between cells: --------------------------
#     div_Jgj = (np.dot(cells.M_sum_mems, Jgj * cells.mem_sa) / cells.cell_vol)
#
#     Phi_gj = np.dot(cells.lapGJ_P_inv, div_Jgj)
#
#     gPhi = (Phi_gj[cells.cell_nn_i[:, 1]] - Phi_gj[cells.cell_nn_i[:, 0]]) / (cells.nn_len)
#
#     # make the field divergence-free:
#     Jgj = Jgj - gPhi
#
#     # enforce the boundary condition:
#     Jgj[cells.bflags_mems] = 0.0
#
#     #--------cell membrane currents ----------------------------------------------------------------------------------
#
#     # correct the membrane ion currents using a Phi potential defined between inside and outside of each cell:
#     # calculate the divergence of the currents at each cell patch:
#     divJn = np.dot(cells.M_sum_mems, Jn * cells.mem_sa) / cells.cell_vol
#
#     # obtain the value of correcting potential gradient across each membrane segment of each cell:
#     sim.gradPhi_mem = np.dot(cells.divCell_inv, divJn)
#
#     # correct the membrane current to make it divergence free:
#     Jn = Jn + sim.gradPhi_mem
#
#     # Total currents -------------------------------------------------------------------------------------------------
#     # Final total current across membranes:
#     sim.Jn = Jn + Jgj
#     # sim.Jn = Jn
#
#     # get components of corrected current at membrane mids:
#     Jnx = sim.Jn * cells.mem_vects_flat[:, 2]
#     Jny = sim.Jn * cells.mem_vects_flat[:, 3]
#
#     # map current to the cell centers:
#     sim.J_cell_x = np.dot(cells.M_sum_mems, Jnx) / cells.num_mems
#     sim.J_cell_y = np.dot(cells.M_sum_mems, Jny) / cells.num_mems
#
#     # multiply final result by membrane surface area to obtain current (negative assigns into cell + for plotting)
#     sim.I_mem = -sim.Jn*cells.mem_sa

    # Current in the environment --------------------------------------------------------------------------------------
    if p.sim_ECM is True:

        # non divergence free current densities in the environment:
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

        # Take the grid gradient of the scaled internal potential:
        gPhix, gPhiy = fd.gradient(Phi, cells.delta)

        # subtract the potential term from the solution to yield the actual current density in the environment:
        sim.J_env_x = J_env_x_o.reshape(cells.X.shape) - gPhix
        sim.J_env_y = J_env_y_o.reshape(cells.X.shape) - gPhiy

        # if p.smooth_level > 0.0:
        #     sim.J_env_x = gaussian_filter(sim.J_env_x, p.smooth_level)
        #     sim.J_env_y = gaussian_filter(sim.J_env_y, p.smooth_level)

        sim.E_env_x = sim.J_env_x*(sim.D_env_weight)
        sim.E_env_y = sim.J_env_y*(sim.D_env_weight)



# Single equation method:
    # # calculate a rate of change of current at the boundary prior to modification of sim.Jn:
    # # initialize an empty cell-centred storge vector
    # d_rho = np.zeros(sim.cdl)
    #
    # # As it's very difficult to do this in terms of raw environmental fluxes, alter the all-important boundary
    # # current in terms of the tight junction relative weight. When TJ_load is ~ 0.1, the system behaves like an
    # # open boundary with low TEP. At TJ_load ~ 0.99, the system progresses towards a closed TJ boundary with high TEP
    # if p.sim_ECM is True:
    #
    #     # For sim_ECM, do this in a way that will be modified with cutting events and TJ interventions:
    #     TJ_load_o = sim.D_env_weight.ravel()[cells.map_mem2ecm]
    #
    #     TJ_load = 1 / (1 + TJ_load_o)
    #
    # else:
    #     # For no ECM, just use the default relative TJ diffusion constant:
    #     TJ_load = 1/(1+ p.D_tj)
    #

    # TJ_load = 1 / (1 + (p.D_tj/0.1))
    # d_rho_m = sim.Jn*TJ_load
    # d_rho_m = Jn

    #
    # cell_I = np.dot(cells.M_sum_mems, d_rho_m*cells.mem_sa)/cells.cell_vol
    #
    # d_rho[cells.bflags_cells] = cell_I[cells.bflags_cells]


    # if p.cluster_open is True:
    #
    #     # The current density needs to be corrected to account for overall charge neutrality of system:
    #     sim.Jn, sim.J_cell_x, sim.J_cell_y = cells.zero_div_cell(sim.Jn, rho = d_rho, bc = d_rho_m[cells.bflags_mems],
    #                                                              open_bounds=True)
    # else:

    #     # The current density needs to be corrected to account for overall charge neutrality of system:
    #     sim.Jn, sim.J_cell_x, sim.J_cell_y = cells.zero_div_cell(sim.Jn, rho = d_rho, bc = d_rho_m[cells.bflags_mems],
    #                                                              open_bounds=True)




