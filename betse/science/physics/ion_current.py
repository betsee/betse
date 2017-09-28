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

    # divergence base:
    divJo = np.dot(cells.M_sum_mems, sim.Jn * cells.mem_sa)

    # divergence of cell current (for use in calculations, not actually true divergence):
    sim.divJ_cell =  -divJo/ cells.cell_sa

    # Electric field calculation:
    divJc = divJo / cells.cell_vol

    # Change to long-range voltage in intracellular spaces:
    sim.v_cell += np.dot(cells.lapGJ_P_inv, divJc)*p.dt

    # Calculate electric field in intracellular spaces:
    Eco = (sim.v_cell[cells.cell_nn_i[:, 1]] - sim.v_cell[cells.cell_nn_i[:, 0]]) / (cells.nn_len)

    Exo = cells.mem_vects_flat[:,2]*Eco
    Eyo = cells.mem_vects_flat[:,3]*Eco

    # average electric field to cell centres
    sim.E_cell_x = np.dot(cells.M_sum_mems, Exo * cells.mem_sa) / cells.cell_sa
    sim.E_cell_y = np.dot(cells.M_sum_mems, Eyo * cells.mem_sa) / cells.cell_sa


    # Current in the environment --------------------------------------------------------------------------------------
    if p.is_ecm is True:

        # calculate the net charge in the environment:
        sim.rho_env = np.dot(sim.zs * p.F, sim.cc_env) + sim.extra_rho_env

        # diffusive component of current densities in the environment:
        J_env_x_o = np.dot(p.F*sim.zs, sim.fluxes_env_x)
        J_env_y_o = np.dot(p.F*sim.zs, sim.fluxes_env_y)

        # add in cell currents
        # J_env_x_o[cells.map_cell2ecm] += sim.J_cell_x
        # J_env_y_o[cells.map_cell2ecm] += sim.J_cell_y

        # reshape the matrix:
        J_env_x_o = J_env_x_o.reshape(cells.X.shape)
        J_env_y_o = J_env_y_o.reshape(cells.X.shape)

        # # Method #2-------------------------------
        #
        # # take the divergence:
        # divJ = fd.divergence(J_env_x_o, J_env_y_o, cells.delta, cells.delta)
        #
        # # calculate the rate of change of polarization (macroscopic) voltage:
        # dPhi = np.dot(cells.lapENV_P_inv, divJ.ravel())
        #
        # dPhi = dPhi.reshape(cells.X.shape)
        #
        # # gradient of the polarization voltage yields the electric field
        # gdPhix, gdPhiy = fd.gradient(dPhi, cells.delta)

        #---------------------------------------------

        #Helmholtz-Hodge decomposition to obtain divergence-free projection of actual currents (zero n_hat at boundary):
        _, sim.J_env_x, sim.J_env_y, BB, Jbx, Jby = stb.HH_Decomp(J_env_x_o, J_env_y_o, cells)

        # add to the environmental voltage:
        sim.v_env += BB*p.dt

        # gradient of the polarization voltage yields the electric field
        gVex, gVey = fd.gradient(sim.v_env.reshape(cells.X.shape), cells.delta)

        sim.E_env_x = -gVex
        sim.E_env_y = -gVey

        # smooth field:
        sim.E_env_x = gaussian_filter(sim.E_env_x, 2)
        sim.E_env_y = gaussian_filter(sim.E_env_y, 2)

        # smooth currents
        sim.J_env_x = gaussian_filter(sim.J_env_x, 2)
        sim.J_env_y = gaussian_filter(sim.J_env_y, 2)


        # Method #1 -----------------------------------------
        # Change in free current density at each membrane:
        sim.Jtx = sim.J_env_x + Jbx
        sim.Jty = sim.J_env_y + Jby

        # map current from extracellular space to membrane normal
        sim.Jme = (sim.Jtx.ravel()[cells.map_mem2ecm] * cells.mem_vects_flat[:, 2] +
               sim.Jty.ravel()[cells.map_mem2ecm] * cells.mem_vects_flat[:, 3])

        # Method #2 ------------------------------------------

        # dPx = -gdPhix.ravel()
        # dPy = -gdPhiy.ravel()
        #
        # # map current from extracellular space to membrane normal
        # sim.Jme = (dPx[cells.map_mem2ecm] * cells.mem_vects_flat[:, 2] +
        #        dPy[cells.map_mem2ecm] * cells.mem_vects_flat[:, 3])

        #-------------------------------------------------------

        # Add in any electric field resulting from an applied voltage:
        sim.E_env_x += -((sim.bound_V['R'] - sim.bound_V['L'])/(cells.xmax - cells.xmin))*np.ones(cells.X.shape)
        sim.E_env_y += -((sim.bound_V['T'] - sim.bound_V['B'])/(cells.ymax - cells.ymin))*np.ones(cells.X.shape)

    else:

        sigma = sim.sigma*sim.D_env_weight

        Jox = np.zeros(len(cells.xypts))
        Joy = np.zeros(len(cells.xypts))

        Jox[cells.map_cell2ecm] = sim.J_cell_x
        Joy[cells.map_cell2ecm] = sim.J_cell_y

        Jox = Jox.reshape(cells.X.shape)
        Joy = Joy.reshape(cells.Y.shape)

        # # recalculate a voltage and field based on a media conductivity that may vary in space:
        divJo = fd.divergence(Jox, Joy, cells.delta, cells.delta)
        dPhi = np.dot(cells.lapENV_P_inv, (divJo/sigma).ravel())

        gdpx, gdpy = fd.gradient(dPhi.reshape(cells.X.shape), cells.delta)

        # this is an approximation for environmental free currents:
        Jox = -gdpx
        Joy = -gdpy

        _, sim.J_env_x, sim.J_env_y, BB, _, _ = stb.HH_Decomp(Jox.reshape(cells.X.shape),
                                                              Joy.reshape(cells.X.shape), cells)

        # incrementally add in change to environmental voltage due to current flux:
        sim.v_env = BB*p.dt

        sim.J_cell_x = Jox.ravel()[cells.map_cell2ecm]
        sim.J_cell_y = Joy.ravel()[cells.map_cell2ecm]

        # dPx = -sim.Chi.ravel()*p.eo*gdpx.ravel()
        # dPy = -sim.Chi.ravel()*p.eo*gdpy.ravel()
        #
        # # map current from extracellular space to membrane normal
        # sim.Jme = (dPx[cells.map_mem2ecm] * cells.mem_vects_flat[:, 2] +
        #        dPy[cells.map_mem2ecm] * cells.mem_vects_flat[:, 3])

        # Jtx = Jax + Jbx
        # Jty = Jay + Jby

        # map current from extracellular space to membrane normal
        sim.Jme = (Jox.ravel()[cells.map_mem2ecm] * cells.mem_vects_flat[:, 2] +
               Joy.ravel()[cells.map_mem2ecm] * cells.mem_vects_flat[:, 3])

























