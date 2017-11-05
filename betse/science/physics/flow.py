#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

import numpy as np
from betse.science import sim_toolbox as stb
from scipy.ndimage.filters import gaussian_filter


def getFlow(sim, cells, p):
    """
    Calculate the electroosmotic fluid flow in the cell and extracellular
    networks using Hagen–Poiseuille "pipe flow" equation.

    """

    # First do extracellular space electroosmotic flow--------------------------------------------------------------

    if p.is_ecm is True:

        #  calculating velocity wrt the raw charge:-------------------------
        #conversion to consider only surface charge on cell membranes (rather than volume density):
        scaleF = cells.memSa_per_envSquare.reshape(cells.X.shape) / cells.delta

        # electrostatic body forces in environment:
        FFx = sim.rho_env.reshape(cells.X.shape)*sim.E_env_x*scaleF
        FFy = sim.rho_env.reshape(cells.X.shape)*sim.E_env_y*scaleF

        # volume forces scaled by water viscocity and environmental weight map defining TJ barrier:
        muFx = ((1/p.mu_water)*sim.D_env_weight)*FFx
        muFy = ((1/p.mu_water)*sim.D_env_weight)*FFy

        # Uxo = np.dot(cells.lapENVinv, -muFx.ravel())
        # Uyo = np.dot(cells.lapENVinv, -muFy.ravel())

        # Electroosmotic velocity in terms of slip velocity with charge at the screening layer:
        Uxo = muFx*(1/sim.ko_env)
        Uyo = muFy*(1/sim.ko_env)

        # Helmholtz-Hodge decomposition to obtain divergence-free projection of flow (zero at external boundary):
        _, Ux, Uy, _, _, _ = stb.HH_Decomp(Uxo.reshape(cells.X.shape),
                                            Uyo.reshape(cells.X.shape), cells)

        sim.u_env_x = Ux
        sim.u_env_y = Uy



    # -------Next do flow through gap junction connected cells-------------------------------------------------------

    # sigma = np.dot((((sim.zs ** 2) * p.q * p.F * sim.D_free) / (p.kb * p.T)), sim.cc_cells)

    # Charge density per unit volume:
    rho_cells = np.dot(sim.zs * p.F, sim.cc_cells) + sim.extra_rho_cells

    # net force is the electrostatic body force on net volume charge in cells:
    Fxc = sim.E_cell_x*rho_cells*(1/p.mu_water)
    Fyc = sim.E_cell_y*rho_cells*(1/p.mu_water)

    # Calculate flow under body forces using Stokes flow:
    u_gj_xo = np.dot(cells.lapGJinv, -Fxc)
    u_gj_yo = np.dot(cells.lapGJinv, -Fyc)

    # Electroosmotic velocity in terms of slip velocity with charge at the screening layer:
    # u_gj_xo = Fxc*(1/sim.ko_cell)
    # u_gj_yo = Fyc*(1/sim.ko_cell)

    # Flow must be made divergence-free: use the Helmholtz-Hodge decomposition method:
    _, sim.u_cells_x, sim.u_cells_y, _, _, _ = cells.HH_cells(u_gj_xo, u_gj_yo, rot_only=True)

