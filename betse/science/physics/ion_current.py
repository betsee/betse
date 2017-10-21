#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

import numpy as np
from scipy.ndimage.filters import gaussian_filter
from betse.science import sim_toolbox as stb
from betse.science.math import finitediff as fd


def get_current(sim, cells, p):

    # calculate the net charge in cells:
    sim.rho_cells = np.dot(sim.zs * p.F, sim.cc_cells) + sim.extra_rho_cells

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
    # sim.divJ_cell =  -np.dot(cells.M_sum_mems, sim.Jn * cells.mem_sa)/ cells.cell_sa
    # conductivity of cells is assumed to be 0.1x lower than that of free environment:
    sim.sigma_cell = np.asarray([((z ** 2) * p.q * p.F * cc * D * 0.1) / (p.kb * p.T) for (z, cc, D) in
                                 zip(sim.zs, sim.cc_cells, sim.D_free)]).mean(axis=0)

    # calculate electric field in cells using net intracellular current and cytosol conductivity:
    sim.E_cell_x = sim.J_cell_x*(1/(sim.sigma_cell))
    sim.E_cell_y = sim.J_cell_y*(1/(sim.sigma_cell))


    # Current in the environment --------------------------------------------------------------------------------------
    if p.is_ecm is True:

        # calculate the net charge in the environment:
        sim.rho_env = np.dot(sim.zs * p.F, sim.cc_env) + sim.extra_rho_env

        # calculate current densities in the environment from ion fluxes:
        J_env_x_o = np.dot(p.F*sim.zs, sim.fluxes_env_x) + sim.extra_Jenv_x
        J_env_y_o = np.dot(p.F*sim.zs, sim.fluxes_env_y) + sim.extra_Jenv_y

        # reshape the matrix:
        J_env_x_o = J_env_x_o.reshape(cells.X.shape)
        J_env_y_o = J_env_y_o.reshape(cells.X.shape)

        #---------------------------------------------

        #Helmholtz-Hodge decomposition to obtain divergence-free projection of actual currents (zero n_hat at boundary):
        _, sim.J_env_x, sim.J_env_y, BB, Jbx, Jby = stb.HH_Decomp(J_env_x_o, J_env_y_o, cells)


        # The solution to the Screened Poisson Equation in the limit of large screening constant Ko, is simply
        # Phi = f/Ko2. This makes a perfect voltage estimate for the extracellular space.

        if p.smooth_level > 0.0:

            v_env = (gaussian_filter(sim.rho_env.reshape(cells.X.shape), p.smooth_level).ravel())/((sim.ko_env**2)*p.er*p.eo)

        else:

            v_env = (sim.rho_env.ravel()) / ((sim.ko_env ** 2) * p.er * p.eo)


        v_env = v_env.reshape(cells.X.shape)

        # add in extra boundary conditions for the case of an externally-applied voltage event:
        v_env[:, -1] = sim.bound_V['R']
        v_env[:, 0] = sim.bound_V['L']
        v_env[-1, :] = sim.bound_V['T']
        v_env[0, :] = sim.bound_V['B']

        # gradient of the polarization voltage yields the electric field:
        gVex, gVey = fd.gradient(v_env, cells.delta)

        # keep the gradient of environmental voltage calculated from raw charge; this works in sim.update_ecm:
        sim.gVex = gVex
        sim.gVey = gVey

        # assign to environmental voltage array:
        sim.v_env = 1*v_env.ravel()

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

        # recalculate the conductivity as it will change every itteration:
        # calculate specific maps of conductivity in cells and environment

        # conductivity map for environment:
        sim.sigma_env = np.asarray(
            [((z ** 2) * p.q * p.F * cc * D) / (p.kb * p.T) for (z, cc, D) in
             zip(sim.zs, sim.cc_env, sim.D_env)]).mean(
            axis=0)

        if p.smooth_level_sigma_map > 0.0:
            # environmental conductivity matrix needs to be smoothed to assured simulation stability:
            sim.sigma_env = gaussian_filter(sim.sigma_env.reshape(cells.X.shape), p.smooth_level_sigma_map)


        # calculate electric field components; these are smoothed and more ideal for other environmental updates:
        sim.E_env_x = sim.Jtx*(1/sim.sigma_env.reshape(cells.X.shape))
        sim.E_env_y = sim.Jty*(1/sim.sigma_env.reshape(cells.X.shape))

        # # Add in any electric field resulting from an applied voltage:
        # sim.E_env_x += -((sim.bound_V['R'] - sim.bound_V['L'])/(cells.xmax - cells.xmin))*np.ones(cells.X.shape)
        # sim.E_env_y += -((sim.bound_V['T'] - sim.bound_V['B'])/(cells.ymax - cells.ymin))*np.ones(cells.X.shape)
        #
        # # Calculate the normal component of electric field at the cell
        #
        # sim.Eme = (sim.E_env_x.ravel()[cells.map_mem2ecm] * cells.mem_vects_flat[:, 2] +
        #        sim.E_env_y.ravel()[cells.map_mem2ecm] * cells.mem_vects_flat[:, 3])

    else:

        sigma = sim.sigma*sim.D_env_weight

        # divergence of transmembrane current represents rate of change of charge in space:
        divJc = np.dot(cells.M_sum_mems, sim.Jn * cells.mem_sa) / cells.cell_vol

        # solving the Poisson equation for divJc yields an approximation for the voltage change:
        venv = np.dot(cells.lapGJ_P_inv, divJc)

        sim.v_env[cells.map_cell2ecm] += venv*p.dt*(1/sim.cedl_cell)

        # gradient of the polarization voltage yields the electric field
        gVex, gVey = fd.gradient(sim.v_env.reshape(cells.X.shape), cells.delta)

        sim.E_env_x = -gVex
        sim.E_env_y = -gVey

        # work backwards to approximate environmental free currents:
        Jox = sim.E_env_x*sigma
        Joy = sim.E_env_y*sigma

        # obtain divergence-free component of current:
        _, sim.J_env_x, sim.J_env_y, _, _, _ = stb.HH_Decomp(Jox.reshape(cells.X.shape),
                                                              Joy.reshape(cells.X.shape), cells)

        # map current from extracellular space to membrane normal
        sim.Jme = (Jox.ravel()[cells.map_mem2ecm] * cells.mem_vects_flat[:, 2] +
               Joy.ravel()[cells.map_mem2ecm] * cells.mem_vects_flat[:, 3])

























