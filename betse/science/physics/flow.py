#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

import numpy as np
from scipy import interpolate as interp
from scipy.ndimage.filters import gaussian_filter
from betse.science import finitediff as fd
from scipy.sparse.linalg import lsmr


# FIXME!!! change over to HH decomposition!

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

            Fe_x = Qenv * sim.E_env_x
            Fe_y = Qenv * sim.E_env_y

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
        P = np.dot(cells.lapENV_P_inv, div_uo.ravel())
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
    P_react = np.dot(cells.lapGJ_P_inv, div_u)

    # calculate its gradient:
    gradP_react = (P_react[cells.cell_nn_i[:, 1]] - P_react[cells.cell_nn_i[:, 0]]) / (cells.nn_len)

    gP_x = gradP_react * cells.mem_vects_flat[:,2]
    gP_y = gradP_react * cells.mem_vects_flat[:,3]

    sim.u_gj_x = u_gj_xo - gP_x
    sim.u_gj_y = u_gj_yo - gP_y

    # enforce zero flow at outer boundary:
    sim.u_gj_x[cells.bflags_mems] = 0
    sim.u_gj_y[cells.bflags_mems] = 0

    # average the components at cell centres:
    sim.u_cells_x = np.dot(cells.M_sum_mems, sim.u_gj_x) / cells.num_mems
    sim.u_cells_y = np.dot(cells.M_sum_mems, sim.u_gj_y) / cells.num_mems


#--------WASTELANDS-----------------------------------------------------------------------------------------------
