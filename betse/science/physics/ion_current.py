#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

import numpy as np
from scipy.ndimage.filters import gaussian_filter
from betse.science import sim_toolbox as stb
from betse.science.math import finitediff as fd


def get_current(sim, cells, p):

    # calculate membrane current density (- as fluxes were defined into cell)
    sim.Jmem = -np.dot(sim.zs * p.F, sim.fluxes_mem) + sim.extra_J_mem

    # calculate current density across cell membranes via gap junctions:
    sim.Jgj = np.dot(sim.zs * p.F, sim.fluxes_gj)

    # add the free current sources together into a single transmembrane current:
    sim.Jn = sim.Jmem + sim.Jgj

    # multiply final result by membrane surface area to obtain current (direction into cell is +)
    sim.I_mem = -sim.Jn*cells.mem_sa

    # components of intracellular current:
    Jcx = sim.Jn * cells.mem_vects_flat[:,2]
    Jcy = sim.Jn * cells.mem_vects_flat[:,3]

    # average intracellular current to cell centres
    sim.J_cell_x = np.dot(cells.M_sum_mems, Jcx*cells.mem_sa) / cells.cell_sa
    sim.J_cell_y = np.dot(cells.M_sum_mems, Jcy*cells.mem_sa) / cells.cell_sa

    # normal component of J_cell at the membranes:
    sim.Jc = sim.J_cell_x[cells.mem_to_cells]*cells.mem_vects_flat[:,2] + sim.J_cell_y[cells.mem_to_cells]*cells.mem_vects_flat[:,3]

    # divergence of cell current (for use in calculations, not actually true divergence):
    sim.divJ_cell =  -np.dot(cells.M_sum_mems, sim.Jn * cells.mem_sa)/ cells.cell_sa

    # this approximation assumes that the intracellular space is slightly lower conductivity than free media:
    sim.E_cell_x = sim.J_cell_x*(1/(sim.sigma*0.1))
    sim.E_cell_y = sim.J_cell_y*(1/(sim.sigma*0.1))


    # Current in the environment --------------------------------------------------------------------------------------
    if p.is_ecm is True:

        # calculate the net charge in the environment:
        sim.rho_env = np.dot(sim.zs * p.F, sim.cc_env) + sim.extra_rho_env

        # diffusive component of current densities in the environment:
        J_env_x_o = np.dot(p.F*sim.zs, sim.fluxes_env_x)
        J_env_y_o = np.dot(p.F*sim.zs, sim.fluxes_env_y)

        # reshape the matrix:
        J_env_x_o = J_env_x_o.reshape(cells.X.shape)
        J_env_y_o = J_env_y_o.reshape(cells.X.shape)

        #---------------------------------------------

        #Helmholtz-Hodge decomposition to obtain divergence-free projection of actual currents (zero n_hat at boundary):
        _, sim.J_env_x, sim.J_env_y, BB, Jbx, Jby = stb.HH_Decomp(J_env_x_o, J_env_y_o, cells)

        # sim.v_env = np.dot(cells.lapENVinv, -(sim.rho_env.ravel() / (sim.ko_env * p.er * p.eo)))

        # The solution to the Screened Poisson Equation in the limit of large screening constant Ko, is simply
        # Phi = f/Ko2. This makes a perfect voltage estimate for the extracellular space.

        if p.smooth_level > 0.0:

            sim.v_env = (gaussian_filter(sim.rho_env.reshape(cells.X.shape), p.smooth_level).ravel())/((sim.ko_env**2)*p.er*p.eo)

        else:

            sim.v_env = (sim.rho_env.ravel()) / ((sim.ko_env ** 2) * p.er * p.eo)

        # gradient of the polarization voltage yields the electric field
        gVex, gVey = fd.gradient(sim.v_env.reshape(cells.X.shape), cells.delta)

        sim.E_env_x = -gVex
        sim.E_env_y = -gVey

        if p.smooth_level > 0.0:

            # smooth currents
            sim.J_env_x = gaussian_filter(sim.J_env_x, p.smooth_level)
            sim.J_env_y = gaussian_filter(sim.J_env_y, p.smooth_level)


        # Method #1 -----------------------------------------
        # Change in free current density at each membrane:
        sim.Jtx = sim.J_env_x + Jbx
        sim.Jty = sim.J_env_y + Jby

        # map current from extracellular space to membrane normal
        sim.Jme = (sim.Jtx.ravel()[cells.map_mem2ecm] * cells.mem_vects_flat[:, 2] +
               sim.Jty.ravel()[cells.map_mem2ecm] * cells.mem_vects_flat[:, 3])


        #-------------------------------------------------------

        # Add in any electric field resulting from an applied voltage:
        sim.E_env_x += -((sim.bound_V['R'] - sim.bound_V['L'])/(cells.xmax - cells.xmin))*np.ones(cells.X.shape)
        sim.E_env_y += -((sim.bound_V['T'] - sim.bound_V['B'])/(cells.ymax - cells.ymin))*np.ones(cells.X.shape)

        # sim.Jme = (sim.E_env_x.ravel()[cells.map_mem2ecm] * cells.mem_vects_flat[:, 2] +
        #        sim.E_env_y.ravel()[cells.map_mem2ecm] * cells.mem_vects_flat[:, 3])

    else:

        sigma = sim.sigma*sim.D_env_weight

        # add in the contribution from cell charge exchange across membranes:
        divJc = np.dot(cells.M_sum_mems, sim.Jn * cells.mem_sa) / cells.cell_vol

        vcell = np.dot(cells.lapGJ_P_inv, divJc)

        sim.v_env = np.zeros(len(cells.xypts))

        sim.v_env[cells.map_cell2ecm] = vcell

        # gradient of the polarization voltage yields the electric field
        gVex, gVey = fd.gradient(sim.v_env.reshape(cells.X.shape), cells.delta)

        sim.E_env_x = -gVex
        sim.E_env_y = -gVey

        # approximate environmental free currents:
        Jox = sim.E_env_x*sigma
        Joy = sim.E_env_y*sigma

        _, sim.J_env_x, sim.J_env_y, BB, _, _ = stb.HH_Decomp(Jox.reshape(cells.X.shape),
                                                              Joy.reshape(cells.X.shape), cells)

        # incrementally add in change to environmental voltage due to current flux in environment:
        # sim.v_env += BB*p.dt*(1/sim.cedl_env)


        # map current from extracellular space to membrane normal
        sim.Jme = (Jox.ravel()[cells.map_mem2ecm] * cells.mem_vects_flat[:, 2] +
               Joy.ravel()[cells.map_mem2ecm] * cells.mem_vects_flat[:, 3])

























