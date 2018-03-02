#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

import numpy as np
from betse.science import sim_toolbox as stb
from scipy.ndimage.filters import gaussian_filter
from betse.science.math import finitediff as fd


def getFlow(sim, cells, p):
    """
    Calculate the electroosmotic fluid flow in the cell and extracellular
    networks using Stokes-Equation

    """

    # First do extracellular space electroosmotic flow--------------------------------------------------------------

    if p.is_ecm is True:

        # Use Helmholtz-Smoluchowski equation:
        muFx = -p.eo*p.er*(1 / p.mu_water) * sim.E_env_x * sim.v_env.reshape(cells.X.shape) * sim.D_env_weight
        muFy = -p.eo*p.er*(1 / p.mu_water) * sim.E_env_y * sim.v_env.reshape(cells.X.shape) * sim.D_env_weight

        uxo = np.dot(cells.lapENVinv, -muFx.ravel())
        uyo = np.dot(cells.lapENVinv, -muFy.ravel())

        _, sim.u_env_x, sim.u_env_y, _, _, _ = stb.HH_Decomp(uxo, uyo, cells)


    # -------Next do flow through gap junction connected cells-------------------------------------------------------

    # net force is the electrostatic body force on net volume charge in cells:
    Fxc = sim.E_cell_x*sim.rho_cells*(1/p.mu_water)*p.gj_surface
    Fyc = sim.E_cell_y*sim.rho_cells*(1/p.mu_water)*p.gj_surface

    # Calculate flow under body forces using Stokes flow:
    u_gj_xo = np.dot(cells.lapGJinv, -Fxc)
    u_gj_yo = np.dot(cells.lapGJinv, -Fyc)

    # Flow must be made divergence-free: use the Helmholtz-Hodge decomposition method:
    _, sim.u_cells_x, sim.u_cells_y, _, _, _ = cells.HH_cells(u_gj_xo, u_gj_yo, rot_only=True)


    #---Wastelands---------------------------------------------------------------------------------------------------
    # --Alternative method #1 calculates fluid flow in terms of slip velocity only:----------------------------------

    # scaleF =  (sim.ko_env*p.cell_height) # conversion concentrating charge density in the double layer
    #
    # # total charge in an env-grid square in Coulombs:
    # q_env = sim.rho_env.reshape(cells.X.shape)*(cells.delta**2)*(p.cell_height)
    #
    # # slip velocity forces:
    # # This assumes that the *screening* charge is pulled on by the self-generated electric field.
    # # Further investigation may be warranted in the future as to the sign of the slip velocity
    # muFx = -(1/p.mu_water)*sim.E_env_x*q_env*scaleF*sim.D_env_weight
    # muFy = -(1/p.mu_water)*sim.E_env_y*q_env*scaleF*sim.D_env_weight
    #
    # _, sim.u_env_x, sim.u_env_y, _, _, _ = stb.HH_Decomp(muFx, muFy, cells)


    # --Alternative method #2 calculates fluid flow in terms of slip velocity only using Helmholtz-Smoluchowski:-----

    # slip velocity forces:
    # This assumes that the *screening* charge is pulled on by the self-generated electric field.
    # Further investigation may be warranted in the future as to the sign of the slip velocity

    # vme = np.zeros(sim.edl)
    # vme[cells.map_mem2ecm] = sim.vm/2
    # vme = vme.reshape(cells.X.shape)
    #
    # muFx = (1/p.mu_water)*sim.E_env_x*(vme)*p.eo*p.er
    # muFy = (1/p.mu_water)*sim.E_env_y*(vme)*p.eo*p.er
    #
    # _, sim.u_env_x, sim.u_env_y, _, _, _ = stb.HH_Decomp(muFx, muFy, cells)

    # Alternative method #3 calculates flow in terms of curl component of current density:--------------------------
    # cc = sim.cc_env.mean(axis=0).reshape(cells.X.shape)
    # zz = sim.zs.mean()
    #
    # sim.u_env_x = -sim.J_env_x / (p.F * cc * zz)
    # sim.u_env_y = -sim.J_env_y / (p.F * cc * zz)

