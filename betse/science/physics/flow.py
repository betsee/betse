#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

import numpy as np
from betse.science import sim_toolbox as stb
from scipy.ndimage.filters import gaussian_filter


def getFlow(sim, cells, p):
    """
    Calculate the electroosmotic fluid flow in the cell and extracellular
    networks using Stokes-Equation

    """

    # First do extracellular space electroosmotic flow--------------------------------------------------------------

    if p.is_ecm is True:

        scaleF = (p.cell_space / cells.delta)

        muFx = (1 / p.mu_water)*sim.E_env_x*sim.rho_env.reshape(cells.X.shape)*scaleF*sim.D_env_weight
        muFy = (1 / p.mu_water)*sim.E_env_y*sim.rho_env.reshape(cells.X.shape)*scaleF*sim.D_env_weight

        uxo = np.dot(cells.lapENVinv, -muFx.ravel())
        uyo = np.dot(cells.lapENVinv, -muFy.ravel())

        _, sim.u_env_x, sim.u_env_y, _, _, _ = stb.HH_Decomp(uxo, uyo, cells)


    # -------Next do flow through gap junction connected cells-------------------------------------------------------

    # Charge density per unit volume:
    rho_cells = np.dot(sim.zs * p.F, sim.cc_cells) + sim.extra_rho_cells

    # net force is the electrostatic body force on net volume charge in cells:
    Fxc = sim.E_cell_x*rho_cells*(1/p.mu_water)
    Fyc = sim.E_cell_y*rho_cells*(1/p.mu_water)

    # Calculate flow under body forces using Stokes flow:
    u_gj_xo = np.dot(cells.lapGJinv, -Fxc)
    u_gj_yo = np.dot(cells.lapGJinv, -Fyc)

    # Flow must be made divergence-free: use the Helmholtz-Hodge decomposition method:
    _, sim.u_cells_x, sim.u_cells_y, _, _, _ = cells.HH_cells(u_gj_xo, u_gj_yo, rot_only=True)

