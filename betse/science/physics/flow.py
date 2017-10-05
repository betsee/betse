#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

import numpy as np
from betse.science import sim_toolbox as stb


def getFlow(sim, cells, p):
    """
    Calculate the electroosmotic fluid flow in the cell and extracellular
    networks using Hagen–Poiseuille "pipe flow" equation.

    """

    # First do extracellular space electroosmotic flow--------------------------------------------------------------

    if p.is_ecm is True:

        # electrostatic body forces in environment:
        FFx = sim.rho_env.reshape(cells.X.shape)*sim.E_env_x
        FFy = sim.rho_env.reshape(cells.X.shape)*sim.E_env_y

        # non-divergence free currents using Stokes flow equation:
        muFx = ((1/p.mu_water)*sim.D_env_weight)*FFx
        muFy = ((1/p.mu_water)*sim.D_env_weight)*FFy

        Uxo = np.dot(cells.lapENVinv, -muFx.ravel())
        Uyo = np.dot(cells.lapENVinv, -muFy.ravel())

        # Helmholtz-Hodge decomposition to obtain divergence-free projection of flow (zero at boundary):
        _, Ux, Uy, _, _, _ = stb.HH_Decomp(Uxo.reshape(cells.X.shape),
                                            Uyo.reshape(cells.X.shape), cells)

        sim.u_env_x = Ux
        sim.u_env_y = Uy


        # Gxo = FFx.ravel()[cells.map_mem2ecm]
        # Gyo = FFy.ravel()[cells.map_mem2ecm]
        #
        # Gx = (np.dot(cells.M_sum_mems, Gxo*cells.mem_sa)/cells.cell_sa)
        # Gy = (np.dot(cells.M_sum_mems, Gyo*cells.mem_sa) / cells.cell_sa)
        #
        # # Calculate flow under body forces using Stokes flow:
        # u_gj_xo = np.dot(cells.lapGJinv, -Gx*(1/p.mu_water))
        # u_gj_yo = np.dot(cells.lapGJinv, -Gy*(1/p.mu_water))
        #
        # # Flow must be made divergence-free: use the Helmholtz-Hodge decomposition method:
        # _, ux, uy, _, _, _ = cells.HH_cells(u_gj_xo, u_gj_yo, rot_only=True)
        #
        # u_env_x = np.zeros(sim.edl)
        # u_env_y = np.zeros(sim.edl)
        #
        # u_env_x[cells.map_cell2ecm] = ux*1
        # u_env_y[cells.map_cell2ecm] = uy*1
        #
        # sim.u_env_x = u_env_x.reshape(cells.X.shape)
        # sim.u_env_y = u_env_y.reshape(cells.X.shape)


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

    # Flow must be made divergence-free: use the Helmholtz-Hodge decomposition method:
    _, sim.u_cells_x, sim.u_cells_y, _, _, _ = cells.HH_cells(u_gj_xo, u_gj_yo, rot_only=True)

