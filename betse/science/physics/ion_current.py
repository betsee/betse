#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

import numpy as np
from scipy.ndimage.filters import gaussian_filter

from betse.science import sim_toolbox as stb
from betse.science.math import finitediff as fd


def get_current(sim, cells, p):

    # nmx = cells.mem_vects_flat[:,2]
    # nmy = cells.mem_vects_flat[:,3]

    # calculate membrane current density (- as fluxes were defined into cell)
    sim.Jmem = -np.dot(sim.zs * p.F, sim.fluxes_mem)

    # calculate current density across cell membranes via gap junctions:
    sim.Jgj = np.dot(sim.zs * p.F, sim.fluxes_gj)

    # add the free current sources together into a single transmembrane current:
    sim.Jn = sim.Jmem + sim.Jgj + sim.extra_J_mem

    # multiply final result by membrane surface area to obtain current (direction into cell is +)
    sim.I_mem = -sim.Jn*cells.mem_sa

    # calculate current density across whole cells via cytoplasm:
    sim.Jc = np.dot(sim.zs*p.F, sim.fluxes_intra)

    # components of intracellular current:
    Jcx = sim.Jc * cells.mem_vects_flat[:,2]
    Jcy = sim.Jc * cells.mem_vects_flat[:,3]

    # average intracellular current to cell centres FIXME I believe the intracellular current should be curl only!:
    sim.J_cell_x = np.dot(cells.M_sum_mems, Jcx*cells.mem_sa) / cells.cell_sa
    sim.J_cell_y = np.dot(cells.M_sum_mems, Jcy*cells.mem_sa) / cells.cell_sa

    # update inverse electrical double layer in cells based on internal concentrations:
    sim.ko_cell = (
    np.sqrt(np.dot((p.NAv * (p.q ** 2) * sim.zs ** 2) / (p.eo * p.er * p.kb * p.T), sim.cc_cells))).mean()


    # Current in the environment --------------------------------------------------------------------------------------
    if p.sim_ECM is True:

        # calculate the net charge in the environment:
        sim.rho_env = np.dot(sim.zs * p.F, sim.cc_env) + sim.extra_rho_env

        # update inverse electrical double layer in environment based on external concentrations:
        sim.ko_env = (
            np.sqrt(np.dot((p.NAv * (p.q ** 2) * sim.zs ** 2) / (p.eo * p.er * p.kb * p.T), sim.cc_env))).mean()

        # currents in environment due to transmembrane fluxes:
        # J_tm_env_x = np.zeros(sim.edl)
        # J_tm_env_x[cells.map_mem2ecm] = sim.Jmem*cells.mem_vects_flat[:, 2]*(cells.mem_sa/cells.ecm_sa).mean()
        #
        # J_tm_env_y = np.zeros(sim.edl)
        # J_tm_env_y[cells.map_mem2ecm] = sim.Jmem*cells.mem_vects_flat[:, 3]*(cells.mem_sa/cells.ecm_sa).mean()

        # diffusive component of current densities in the environment:
        J_env_x_o = np.dot(p.F*sim.zs, sim.fluxes_env_x)
        J_env_y_o = np.dot(p.F*sim.zs, sim.fluxes_env_y)

        # add in component of current from transmembrane fluxes:
        # J_env_x_o += J_tm_env_x
        # J_env_y_o += J_tm_env_y

        # reshape the matrix:
        J_env_x_o = J_env_x_o.reshape(cells.X.shape)
        J_env_y_o = J_env_y_o.reshape(cells.X.shape)

        # conductivity in the media is modified by the environmental diffusion weight matrix:
        # sigma = np.dot((((sim.zs ** 2) * p.q * p.F) / (p.kb * p.T)), sim.cc_env*sim.D_env).reshape(cells.X.shape)

        #---Calculate divergences for concentration & transmembrane fluxes ---------------------------------------------
        # div_Jo = fd.divergence(J_env_x_o/sigma, J_env_y_o/sigma, cells.delta, cells.delta)

        # term describing source of environmental potential:
        source_term =  (sim.rho_env/((sim.ko_env)*p.eo*p.eedl)).reshape(cells.X.shape)

        # term producing vector spherical harmonics (becomes equivalent to Vector Helmholtz Equation):
        # source_term = div_Jo + (sim.rho_env / ((sim.ko_env**2) * p.eo * p.er)).reshape(cells.X.shape)

        # set boundary conditions
        source_term[:,0] = -sim.bound_V['L']*(1/cells.delta**2)
        source_term[:,-1] = -sim.bound_V['R']*(1/cells.delta**2)
        source_term[0,:] = -sim.bound_V['B']*(1/cells.delta**2)
        source_term[-1,:] = -sim.bound_V['T']*(1/cells.delta**2)

        # Calculate a voltage that resists the divergence:
        Phi = np.dot(cells.lapENVinv, -source_term.ravel())

        if p.smooth_level > 0.0:

            # smoothing of Phi:
            Phi = gaussian_filter(Phi.reshape(cells.X.shape), p.smooth_level, mode='constant')

        sim.v_env = Phi.ravel()

        #--------------------------------------------------------------------------------------------------------------
        # calculate the gradient of v_env:
        gPhix, gPhiy = fd.gradient(Phi.reshape(cells.X.shape), cells.delta)

        sim.E_env_x = -gPhix
        sim.E_env_y = -gPhiy

        #Helmholtz-Hodge decomposition to obtain divergence-free projection of actual currents (zero n_hat at boundary):
        _, sim.J_env_x, sim.J_env_y, _, _, _ = stb.HH_Decomp(J_env_x_o, J_env_y_o, cells)

        # sim.J_env_x = J_env_x_o
        # sim.J_env_y = J_env_y_o





















