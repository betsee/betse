#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

import numpy as np
from betse.science import finitediff as fd
from scipy.ndimage.filters import gaussian_filter
from betse.science import sim_toolbox as stb

def get_current(sim, cells, p):


    # membrane normal vectors short form:
    nx = cells.mem_vects_flat[:, 2]
    ny = cells.mem_vects_flat[:, 3]

    # calculate membrane current density (- as fluxes were defined into cell)
    Jn = -np.dot(sim.zs * p.F, sim.fluxes_mem)

    # calculate current density across cell membranes via gap junctions:
    Jgj = np.dot(sim.zs * p.F, sim.fluxes_gj)

    # # calculate the change in polarization from the displacement current:
    # sim.Pol_mem = sim.Pol_mem + p.dt*sim.dvm*p.cm*p.cell_polarizability

    # add the free current sources together into a single transmembrane current:
    sim.Jn = Jn + Jgj

    # multiply final result by membrane surface area to obtain current (direction into cell is +)
    sim.I_mem = -sim.Jn*cells.mem_sa

    # components of current:
    Jnx = sim.Jn * nx
    Jny = sim.Jn * ny

    # average to cell centres:
    sim.J_cell_x = np.dot(cells.M_sum_mems, Jnx) / cells.num_mems
    sim.J_cell_y = np.dot(cells.M_sum_mems, Jny) / cells.num_mems

    Jgj_cell_x = np.dot(cells.M_sum_mems, Jgj*nx) / cells.num_mems
    Jgj_cell_y = np.dot(cells.M_sum_mems, Jgj*ny) / cells.num_mems

    Jm_cell_x = np.dot(cells.M_sum_mems, Jn*nx) / cells.num_mems
    Jm_cell_y = np.dot(cells.M_sum_mems, Jn*ny) / cells.num_mems




    # Current in the environment --------------------------------------------------------------------------------------
    if p.sim_ECM is True:

        # non divergence free current densities in the environment:
        # J_env_x_o = np.dot(p.F*sim.zs, sim.fluxes_env_x)
        # J_env_y_o = np.dot(p.F*sim.zs, sim.fluxes_env_y)

        J_env_x_o = np.zeros(sim.edl)
        J_env_y_o = np.zeros(sim.edl)

        # # base current in environment is equivalent to transmembrane current of cells:
        # J_env_x_o[cells.map_cell2ecm] = J_env_x_o[cells.map_cell2ecm] + sim.J_cell_x*(
        #     cells.memSa_per_envSquare[cells.map_cell2ecm]/cells.delta**2)
        #
        # J_env_y_o[cells.map_cell2ecm] = J_env_y_o[cells.map_cell2ecm] + sim.J_cell_y*(
        #     cells.memSa_per_envSquare[cells.map_cell2ecm]/cells.delta**2)

        # J_env_x_o[cells.map_mem2ecm] = J_env_x_o[cells.map_mem2ecm] + Jnx*(cells.memSa_per_envSquare[cells.map_mem2ecm]/cells.delta**2)
        # J_env_y_o[cells.map_mem2ecm] = J_env_y_o[cells.map_mem2ecm] + Jny*(cells.memSa_per_envSquare[cells.map_mem2ecm]/cells.delta**2)

        J_env_x_o[cells.map_cell2ecm] = J_env_x_o[cells.map_cell2ecm] - Jgj_cell_x*(
            cells.memSa_per_envSquare[cells.map_cell2ecm]/cells.delta**2) + Jm_cell_x*(
            cells.memSa_per_envSquare[cells.map_cell2ecm]/cells.delta**2)

        J_env_y_o[cells.map_cell2ecm] = J_env_y_o[cells.map_cell2ecm] - Jgj_cell_y*(
            cells.memSa_per_envSquare[cells.map_cell2ecm]/cells.delta**2)  +  Jm_cell_y*(
            cells.memSa_per_envSquare[cells.map_cell2ecm]/cells.delta**2)



        # if p.smooth_level > 0.0:
        #     J_env_x_o = gaussian_filter(J_env_x_o.reshape(cells.X.shape), p.smooth_level).ravel()
        #     J_env_y_o = gaussian_filter(J_env_y_o.reshape(cells.X.shape), p.smooth_level).ravel()

        AA, Fx, Fy, BB, _, _ = stb.HH_Decomp(J_env_x_o.reshape(cells.X.shape), J_env_y_o.reshape(cells.X.shape), cells)


        # sim.v_env = -BB.ravel()/p.media_sigma

        sim.J_env_x = Fx.reshape(cells.X.shape)
        sim.J_env_y = Fy.reshape(cells.X.shape)

        sim.E_env_x = sim.J_env_x*p.media_sigma
        sim.E_env_y = sim.J_env_y*p.media_sigma





        # # determine corrected current densities assuming
        # # bulk electrolyte neutrality of entire system:
        #
        # # Next, calculate the divergence of the environmental current density:
        # div_J_env_o = fd.divergence(J_env_x_o.reshape(cells.X.shape), J_env_y_o.reshape(cells.X.shape),
        #     cells.delta, cells.delta)
        #
        #
        # # Next, calculate the divergence of the environmental current density divided by media conductivity:
        # # weight_factor = p.media_sigma*sim.D_env_weight
        # #
        # # div_J_env_o = fd.divergence((1/weight_factor)*J_env_x_o.reshape(cells.X.shape),
        # #                             (1/weight_factor)*J_env_y_o.reshape(cells.X.shape),
        # #                             cells.delta, cells.delta)
        #
        # # Find the value of the environmental electric potential:
        # Phi = np.dot(cells.lapENV_P_inv, -div_J_env_o.ravel())
        # Phi = Phi.reshape(cells.X.shape)
        #
        # sim.Phi_env = Phi
        #
        # sim.v_env = np.copy(sim.Phi_env.ravel())
        # # sim.v_env[cells.inds_env] = 0.0
        # # sim.v_env[cells.ecm_bound_k] = 0.0
        #
        # # Take the grid gradient of the scaled internal potential:
        # gPhix, gPhiy = fd.gradient(Phi, cells.delta)
        #
        # # subtract the potential term from the solution to yield the actual current density in the environment:
        # # sim.J_env_x = (J_env_x_o.reshape(cells.X.shape) + (weight_factor)*gPhix)
        # # sim.J_env_y = (J_env_y_o.reshape(cells.X.shape) + (weight_factor)*gPhiy)
        #
        # sim.J_env_x = (J_env_x_o.reshape(cells.X.shape) + gPhix)
        # sim.J_env_y = (J_env_y_o.reshape(cells.X.shape) + gPhiy)


        # Calculate a divergence-free net current, which assumes charge compensation of the bulk electrolyte:
        # This method uses the Helmholtz-Hodge decomposition for a 2D vector field:

        # Jxr = -J_env_y_o.reshape(cells.X.shape)
        # Jyr = J_env_x_o.reshape(cells.X.shape)
        #
        # divJr = fd.divergence(Jxr, Jyr, cells.delta, cells.delta)
        #
        # AA = np.dot(cells.lapENVinv, -divJr.ravel())
        #
        # gAx, gAy = fd.gradient(AA.reshape(cells.X.shape), cells.delta)
        #
        # sim.J_env_x = -gAy
        # sim.J_env_y = gAx

        # sim.E_env_x = p.media_sigma*sim.J_env_x
        # sim.E_env_y = p.media_sigma*sim.J_env_y







