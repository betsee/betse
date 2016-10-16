#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

import numpy as np
from betse.science import finitediff as fd
from scipy.ndimage.filters import gaussian_filter

def get_current(sim, cells, p):

    # FIXME change current over to HH decomposition!

    # Helmholtz-Decomposition method ----------------------------------------------------------------------------

    # # membrane normal vectors short form:
    # nx = cells.mem_vects_flat[:, 2]
    # ny = cells.mem_vects_flat[:, 3]
    #
    # # calculate membrane current density (- as fluxes were defined into cell)
    # Jn = -np.dot(sim.zs * p.F, sim.fluxes_mem)
    #
    # Jnx = Jn * nx
    # Jny = Jn * ny
    #
    # # calculate current density across cell membranes via gap junctions:
    # Jgj = np.dot(sim.zs * p.F, sim.fluxes_gj)
    #
    # Jgjx = Jgj * nx
    # Jgjy = Jgj * ny
    #
    # # add the current sources together into a single transmembrane current:
    # # Jn_o = Jn + Jgj
    #
    # # obtain the divergence of the current density across membranes:
    # div_Jo = np.dot(cells.M_sum_mems, Jn * cells.mem_sa) / cells.cell_vol
    #
    # # calculate the scalar potential:
    # # note: if using only one component of the field, don't need the extra factor of 2
    # Phi = np.dot(cells.lapGJinv, (div_Jo) / (2*cells.geom_weight))
    # # Phi = np.dot(cells.lapGJinv, (div_Jo))
    #
    # # calculate the gradient of the scalar potential between cell centres:
    # gPhi = (Phi[cells.mem_to_cells][cells.nn_i] - Phi[cells.mem_to_cells][cells.mem_i]) / cells.nn_len
    #
    # gPhi[cells.bflags_mems] = Jn[cells.bflags_mems]
    #
    # # calculate the gradient wrt cell membranes and integrate wrt membrane thickness:
    # gPhi_vm = (Phi[cells.mem_to_cells][cells.nn_i] - Phi[cells.mem_to_cells][cells.mem_i])
    # # gPhi_vm[cells.bflags_mems] = 2*Jn_o[cells.bflags_mems]*cells.nn_len.mean()
    #
    # # components of curl-free current:
    # gPhix = gPhi * nx
    # gPhiy = gPhi * ny
    #
    # # pi/2 rotation matrix factor:-----------------------------------------------
    # a = np.sqrt(2) / 2
    #
    # # transform the current by rotating each vector by pi/2:
    # Jgjxr = a * (Jgjx - Jgjy)
    # Jgjyr = a * (Jgjx + Jgjy)
    #
    # # recalculate rotated current vectors membrane normal component
    # Jr = Jgjxr * nx + Jgjyr * ny
    #
    # # calculate the divergence:
    # div_Jr = np.dot(cells.M_sum_mems, Jr * cells.mem_sa) / cells.cell_vol
    #
    # # calculate the vector potential z-component:
    # A_cell = np.dot(cells.lapGJinv, div_Jr / (2*cells.geom_weight))
    # # A_cell = np.dot(cells.lapGJinv, div_Jr)
    #
    # # calcualte the gradient of the vector potential wrt cell centers (zero at outer boundary):
    # gA = (A_cell[cells.mem_to_cells][cells.nn_i] - A_cell[cells.mem_to_cells][cells.mem_i]) / cells.nn_len
    #
    # # and with respect to membranes:
    # gA_vm = (A_cell[cells.mem_to_cells][cells.nn_i] - A_cell[cells.mem_to_cells][cells.mem_i])
    #
    # # calculate comonents of the field:
    # gAxo = gA * nx
    # gAyo = gA * ny
    #
    # # calculate rotation of the field by Pi/2:
    # gAx = a * (gAxo - gAyo)
    # gAy = a * (gAxo + gAyo)
    #
    # J_cell_x = np.dot(cells.M_sum_mems, Jgjx) / cells.num_mems
    # J_cell_y = np.dot(cells.M_sum_mems, Jgjy) / cells.num_mems
    #
    # J_mems_x = np.dot(cells.M_sum_mems, gPhix) / cells.num_mems
    # J_mems_y = np.dot(cells.M_sum_mems, gPhiy) / cells.num_mems
    #
    # # alternative way to calculate the divergence-free current components:
    # # J_rot_x, J_rot_y, _ = cells.curl(0, 0, A_cell)
    #
    # J_gj_x = np.dot(cells.M_sum_mems, gAx) / cells.num_mems
    # J_gj_y = np.dot(cells.M_sum_mems, gAy) / cells.num_mems


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

    # Jn_corr = Jn_o - gPhi
    Jn_corr = Jn_o

    Jnx = Jn_corr*cells.mem_vects_flat[:, 2]
    Jny = Jn_corr*cells.mem_vects_flat[:, 3]

    # map current to the cell centers:
    sim.J_cell_x = np.dot(cells.M_sum_mems, Jnx) / cells.num_mems
    sim.J_cell_y = np.dot(cells.M_sum_mems, Jny) / cells.num_mems

    # # reassign the normal current to the un_corrected component
    sim.Jn = Jn_o*1

    # assign the normal current to the corrected component
    # sim.Jn = Jn_corr*1

    # multiply final result by membrane surface area to obtain current (direction into cell is +)
    sim.I_mem = -sim.Jn*cells.mem_sa

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






