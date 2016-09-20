#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

import numpy as np
from scipy import interpolate as interp
from scipy.ndimage.filters import gaussian_filter
from betse.science import finitediff as fd
from scipy.sparse.linalg import lsmr


def getFlow(sim, cells, p):   # FIXME env flow should use MACs grid formalism
    """
    Calculate the electroosmotic fluid flow in the cell and extracellular
    networks using Hagen–Poiseuille "pipe flow" equation.

    """
    # # treat osmotic flow, which is mass flow across membranes:
    # if p.deform_osmo is True:
    #
    #     # u_osmo is negative as it's defined + into the cell in pressures.py
    #     u_osmo_cells = -sim.u_net
    #
    #     if p.sim_ECM is True:
    #
    #         unetxo = sim.u_net*cells.mem_vects_flat[:,2]
    #         unetyo = sim.u_net*cells.mem_vects_flat[:,3]
    #
    #         unetx = unetxo[cells.map_mem2ecm]
    #         unety = unetyo[cells.map_mem2ecm]
    #
    #         u_osmo_envx = np.zeros(len(cells.xypts))
    #         u_osmo_envx[cells.envInds_inClust] = unetx
    #
    #         u_osmo_envy = np.zeros(len(cells.xypts))
    #         u_osmo_envy[cells.envInds_inClust] = unety
    #
    #         u_osmo_envx = u_osmo_envx.reshape(cells.X.shape)
    #         u_osmo_envy = u_osmo_envy.reshape(cells.X.shape)
    #
    # else:
    #     u_osmo_cells = np.zeros(sim.mdl)
    #
    #     if p.sim_ECM is True:
    #
    #         u_osmo_env = np.zeros(sim.edl)

    # First do extracellular space electroosmotic flow--------------------------------------------------------------

    if p.sim_ECM is True:

        btag = 'closed'

        # estimate the inverse viscosity for extracellular flow based on the diffusion constant weighting
        # for the world:
        alpha = ((p.cell_space ** 2) / (p.mu_water * 32)) * sim.D_env_weight

        if p.deform_electro is True:

            Qenv = sim.rho_env.reshape(cells.X.shape)

            Fe_x = Qenv * sim.E_env_x * sim.field_mod
            Fe_y = Qenv * sim.E_env_y * sim.field_mod

        else:

            Fe_x = np.zeros(cells.X.shape)
            Fe_y = np.zeros(cells.Y.shape)

        # sum the forces:
        Fx = Fe_x
        Fy = Fe_y

        # calculate the base fluid flow using the Hagen-Poiseuille equation:
        ux_ecm_o = gaussian_filter(Fx * alpha, p.smooth_level)
        uy_ecm_o = gaussian_filter(Fy * alpha, p.smooth_level)

        # calculate the divergence of the flow field as the sum of the two spatial derivatives:
        div_uo = fd.divergence(ux_ecm_o, uy_ecm_o, cells.delta, cells.delta)

        # calculate the alpha-scaled internal pressure from the divergence of the force:
        P = np.dot(cells.lapENVinv, div_uo.ravel())
        P = P.reshape(cells.X.shape)

        # enforce zero normal gradient boundary conditions on P:
        P[:, 0] = P[:, 1]
        P[:, -1] = P[:, -2]
        P[0, :] = P[1, :]
        P[-1, :] = P[-2, :]

        # Take the grid gradient of the scaled internal pressure:
        gPx, gPy = fd.gradient(P, cells.delta)

        # subtract the pressure term from the solution to yield a divergence-free flow field
        u_env_x = ux_ecm_o.reshape(cells.X.shape) - gPx
        u_env_y = uy_ecm_o.reshape(cells.X.shape) - gPy

        # velocities at grid-cell centres:
        sim.u_env_x = u_env_x[:]
        sim.u_env_y = u_env_y[:]

        # boundary conditions reinforced:
        sim.u_env_x[:, 0] = 0
        # right
        sim.u_env_x[:, -1] = 0
        # top
        sim.u_env_x[-1, :] = 0
        # bottom
        sim.u_env_x[0, :] = 0

        # left
        sim.u_env_y[:, 0] = 0
        # right
        sim.u_env_y[:, -1] = 0
        # top
        sim.u_env_y[-1, :] = 0
        # bottom
        sim.u_env_y[0, :] = 0

    # -------Next do flow through gap junction connected cells-------------------------------------------------------

    # calculate the inverse viscosity for the cell collection, which is scaled by gj state:
    # alpha_gj = (1 / (32 * p.mu_water)) * ((sim.gjopen * 8e-10) ** 2)

    # approximate radius of gap junctions:
    gj_rad = np.sqrt((cells.mem_sa*p.gj_surface)/3.14)

    alpha_gj = (1 / (32 * p.mu_water)) * ((sim.gjopen * gj_rad) ** 2)


    if p.deform_electro is True:

        Fe_cell_x = sim.F_gj_x
        Fe_cell_y = sim.F_gj_y

    else:

        Fe_cell_x = np.zeros(len(cells.mem_i))
        Fe_cell_y = np.zeros(len(cells.mem_i))

    if p.deform_osmo is True:

        F_osmo_x = sim.F_hydro_x_gj
        F_osmo_y = sim.F_hydro_y_gj

    else:

        F_osmo_x = np.zeros(len(cells.mem_i))
        F_osmo_y = np.zeros(len(cells.mem_i))

    # net force is the sum of electrostatic and hydrostatic pressure induced body forces:
    F_net_x = Fe_cell_x + F_osmo_x
    F_net_y = Fe_cell_y + F_osmo_y

    # Calculate flow under body forces:
    u_gj_xo = F_net_x * alpha_gj
    u_gj_yo = F_net_y * alpha_gj

    # calculate component of velocity normal to cell membranes:
    # u_gj = np.sqrt(u_gj_xo ** 2 + u_gj_yo ** 2)
    u_gj = u_gj_xo * cells.mem_vects_flat[:, 2] + u_gj_yo * cells.mem_vects_flat[:, 3]

    # calculate divergence as the sum of this vector x each surface area, divided by cell volume:
    div_u = (np.dot(cells.M_sum_mems, u_gj * cells.mem_sa) / cells.cell_vol)

    # calculate the reaction pressure required to counter-balance the flow field:

    # if p.deformation is False:
        # if we're not doing deformation, solve with dot product as it's much faster
    P_react = np.dot(cells.lapGJ_P_inv, div_u/2)

    # else:
    #     # if running deformation, solve with scipy's lsmr iterative solver:
    #     P_react = lsmr(cells.lapGJ_P, div_u)[0]

    # calculate its gradient:
    gradP_react = (P_react[cells.cell_nn_i[:, 1]] - P_react[cells.cell_nn_i[:, 0]]) / (cells.nn_len)

    gP_x = gradP_react * cells.mem_vects_flat[:,2]
    gP_y = gradP_react * cells.mem_vects_flat[:,3]

    # sim.u_gj_x = u_gj_xo - gP_x
    # sim.u_gj_y = u_gj_yo - gP_y

    sim.u_gj_x = u_gj_xo
    sim.u_gj_y = u_gj_yo

    # average the components at cell centres:
    sim.u_cells_x = np.dot(cells.M_sum_mems, sim.u_gj_x) / cells.num_mems
    sim.u_cells_y = np.dot(cells.M_sum_mems, sim.u_gj_y) / cells.num_mems

    # # enforce the boundary conditions:
    # sim.u_cells_x[cells.bflags_cells] = 0
    # sim.u_cells_y[cells.bflags_cells] = 0


#--------WASTELANDS-----------------------------------------------------------------------------------------------

# def getFlow_o(sim ,cells ,p):
#
#     """
#     Calculate the electroosmotic fluid flow in the cell and extracellular
#     networks using Stokes Flow equation.
#
#     This function is not presently used in BETSE.
#
#     """
#
#     if p.sim_ECM is True:
#
#         # force of gravity:
#         if p.closed_bound is True:
#
#             btag = 'closed'
#
#         else:
#
#             btag = 'open'
#
#         # estimate the inverse viscosity for extracellular flow based on the diffusion constant weighting
#         # for the world:
#         alpha = ( 1 /p.mu_water ) * sim. D_env_weight * 1.0e-6
#
#         E_ave_x = sim.E_env_x
#         E_ave_y = sim.E_env_y
#
#         if p.deform_electro is True:
#
#             # determine the geometric factor to map charge density from volume to surface:
#             # Qfactor = p.cell_space
#             #
#             # Fe_x = Qfactor*self.rho_env.reshape(cells.X.shape)*E_ave_x
#             # Fe_y = Qfactor*self.rho_env.reshape(cells.X.shape)*E_ave_y
#
#             # map charge density to rectangular grid volume:
#             Qenv = sim. rho_env * cells.ave2ecmV
#
#             # Qenv = Qenv.reshape(cells.X.shape)
#
#             Qenv = gaussian_filter(Qenv.reshape(cells.X.shape) ,p.smooth_level)
#
#             Fe_x = Qenv* E_ave_x
#             Fe_y = Qenv * E_ave_y
#
#
#         else:
#
#             Fe_x = np.zeros(cells.X.shape)
#             Fe_y = np.zeros(cells.Y.shape)
#
#         # sum the forces:
#         Fx = Fe_x
#         Fy = Fe_y
#
#         source_x = -Fx * alpha
#         # source_x = cells.grid_obj.grid_int(source_x,bounds=btag) # perform finite volume integration of source
#
#         source_y = -Fy * alpha
#         # source_y = cells.grid_obj.grid_int(source_y,bounds=btag) # perform finite volume integration of source
#
#         # # calculated the fluid flow using the time-independent Stokes Flow equation:
#         ux_ecm_o = np.dot(cells.lapENVinv, source_x.ravel())
#         uy_ecm_o = np.dot(cells.lapENVinv, source_y.ravel())
#
#         # calculate the divergence of the flow field as the sum of the two spatial derivatives:
#         div_uo = fd.divergence(ux_ecm_o.reshape(cells.X.shape), uy_ecm_o.reshape(cells.X.shape),
#             cells.delta, cells.delta)
#
#         # perform finite volume integration on the divergence:
#         # div_uo = fd.integrator(div_uo)
#         # div_uo = cells.grid_obj.grid_int(div_uo,bounds=btag)
#
#         # calculate the alpha-scaled internal pressure from the divergence of the force:
#         P = np.dot(cells.lapENV_P_inv, div_uo.ravel())
#         P = P.reshape(cells.X.shape)
#
#         # enforce zero normal gradient boundary conditions on P:
#         P[:, 0] = P[:, 1]
#         P[:, -1] = P[:, -2]
#         P[0, :] = P[1, :]
#         P[-1, :] = P[-2, :]
#
#         # Take the grid gradient of the scaled internal pressure:
#         gPx, gPy = fd.gradient(P, cells.delta)
#
#         # subtract the pressure term from the solution to yield a divergence-free flow field
#         u_env_x = ux_ecm_o.reshape(cells.X.shape) - gPx
#         u_env_y = uy_ecm_o.reshape(cells.X.shape) - gPy
#
#         # velocities at cell centres:
#         sim.u_env_x = u_env_x[:]
#         sim.u_env_y = u_env_y[:]
#
#         sim.u_env_x[:, 0] = 0
#         # right
#         sim.u_env_x[:, -1] = 0
#         # top
#         sim.u_env_x[-1, :] = 0
#         # bottom
#         sim.u_env_x[0, :] = 0
#
#         # left
#         sim.u_env_y[:, 0] = 0
#         # right
#         sim.u_env_y[:, -1] = 0
#         # top
#         sim.u_env_y[-1, :] = 0
#         # bottom
#         sim.u_env_y[0, :] = 0
#
#     # ---------------Flow through gap junction connected cells-------------------------------------------------------
#
#     # calculate the inverse viscocity for the cell collection, which is scaled by gj conductivity:
#     alpha_gj_o = p.gj_surface * sim.gjopen * (1 / p.mu_water)
#     # average to individual cells:
#     alpha_gj = np.dot(cells.M_sum_mems, alpha_gj_o) / cells.num_mems
#
#     if p.deform_electro is True:
#
#         Fe_cell_x = sim.F_electro_x
#         Fe_cell_y = sim.F_electro_y
#
#     else:
#         Fe_cell_x = np.zeros(len(cells.cell_i))
#         Fe_cell_y = np.zeros(len(cells.cell_i))
#
#     if p.deform_osmo is True:
#
#         F_osmo_x = sim.F_hydro_x
#         F_osmo_y = sim.F_hydro_y
#
#     else:
#
#         F_osmo_x = np.zeros(len(cells.cell_i))
#         F_osmo_y = np.zeros(len(cells.cell_i))
#
#     # net force is the sum of electrostatic and hydrostatic pressure induced body forces:
#     F_net_x = Fe_cell_x + F_osmo_x
#     F_net_y = Fe_cell_y + F_osmo_y
#
#     # integrate body forces:
#     # F_net_x = cells.integrator(F_net_x)
#     # F_net_y = cells.integrator(F_net_y)
#
#     # Calculate flow under body forces:
#     u_gj_xo = np.dot(cells.lapGJinv, -alpha_gj * F_net_x)
#     u_gj_yo = np.dot(cells.lapGJinv, -alpha_gj * F_net_y)
#
#     # calculate divergence of the flow field using general definition:
#     # first interpolate flow field at membrane midpoints:
#     ux_mem = interp.griddata((cells.cell_centres[:, 0], cells.cell_centres[:, 1]), u_gj_xo,
#         (cells.mem_mids_flat[:, 0], cells.mem_mids_flat[:, 1]), fill_value=0)
#
#     uy_mem = interp.griddata((cells.cell_centres[:, 0], cells.cell_centres[:, 1]), u_gj_yo,
#         (cells.mem_mids_flat[:, 0], cells.mem_mids_flat[:, 1]), fill_value=0)
#
#     # get the component of the velocity field normal to the membranes:
#     u_n = ux_mem * cells.mem_vects_flat[:, 2] + uy_mem * cells.mem_vects_flat[:, 3]
#
#     # calculate divergence as the sum of this vector x each surface area, divided by cell volume:
#     div_u = (np.dot(cells.M_sum_mems, u_n * cells.mem_sa) / cells.cell_vol)
#
#     # calculate the reaction pressure required to counter-balance the flow field:
#     P_react = np.dot(cells.lapGJ_P_inv, 2 * div_u)
#
#     # calculate its gradient:
#     gradP_react = (P_react[cells.cell_nn_i[:, 1]] - P_react[cells.cell_nn_i[:, 0]]) / (cells.nn_len)
#
#     gP_x = gradP_react * cells.cell_nn_tx
#     gP_y = gradP_react * cells.cell_nn_ty
#
#     # average the components of the reaction force field at cell centres and get boundary values:
#     gPx_cell = np.dot(cells.M_sum_mems, gP_x) / cells.num_mems
#     gPy_cell = np.dot(cells.M_sum_mems, gP_y) / cells.num_mems
#
#     sim.u_cells_x = u_gj_xo - gPx_cell
#     sim.u_cells_y = u_gj_yo - gPy_cell
#
#     # enforce the boundary conditions:
#     sim.u_cells_x[cells.bflags_cells] = 0
#     sim.u_cells_y[cells.bflags_cells] = 0