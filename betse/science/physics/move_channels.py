#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

import numpy as np
from betse.science import sim_toolbox as stb


def eosmosis(sim, cells, p):
    """
    Electroosmosis of ion pumps and channels to potentially create directional fluxes in individual cells.

    This is presently simulated by calculating the Nernst-Planck concentration flux of a weighting
    agent, rho, which moves under its own concentration gradient and
    through the influence of the extracellular voltage gradient and fluid flows tangential to the membrane.

    """

    # components of fluid flow velocity at the membrane:
    if p.fluid_flow is True and p.sim_ECM is True:
        ux_mem = sim.u_env_x.ravel()[cells.map_mem2ecm]
        uy_mem = sim.u_env_y.ravel()[cells.map_mem2ecm]

    else:
        ux_mem = 0
        uy_mem = 0

    # get the gradient of rho concentration around each membrane:

    rho_pump = np.dot(cells.M_sum_mems, sim.rho_pump) / cells.num_mems
    rho_channel = np.dot(cells.M_sum_mems, sim.rho_channel) / cells.num_mems

    grad_c = (rho_pump[cells.cell_nn_i[:, 0]] - rho_pump[cells.cell_nn_i[:, 1]]) / cells.nn_len
    grad_c_ch = (rho_channel[cells.cell_nn_i[:, 0]] - rho_channel[cells.cell_nn_i[:, 1]]) / cells.nn_len

    # get the gradient components:
    gcx = grad_c * cells.cell_nn_tx
    gcy = grad_c * cells.cell_nn_ty

    gcx_ch = grad_c_ch * cells.cell_nn_tx
    gcy_ch = grad_c_ch * cells.cell_nn_ty

    # total average electric field at each membrane
    if p.sim_ECM is True:

        Ex = sim.E_env_x.ravel()[cells.map_mem2ecm]
        Ey = sim.E_env_y.ravel()[cells.map_mem2ecm]

        # Ex = self.E_env_x.ravel()[cells.map_mem2ecm]
        # Ey = self.E_env_y.ravel()[cells.map_mem2ecm]

    else:
        Ex = sim.E_gj_x
        Ey = sim.E_gj_y

    # calculate the total Nernst-Planck flux at each membrane for rho_pump factor:

    fx_pump, fy_pump = stb.nernst_planck_flux(sim.rho_pump, gcx, gcy, Ex, Ey, ux_mem, uy_mem, p.D_membrane, p.z_pump,
        sim.T, p)

    # component of total flux in direction of membrane
    # ftot = fx*cells.mem_vects_flat[:,4] + fy*cells.mem_vects_flat[:,5]

    # map the flux to cell centres:
    fx_pump_o = np.dot(cells.M_sum_mems, fx_pump) / cells.num_mems
    fy_pump_o = np.dot(cells.M_sum_mems, fy_pump) / cells.num_mems

    # divergence of the total flux:
    gfx_o = (fx_pump_o[cells.cell_nn_i[:, 1]] - fx_pump_o[cells.cell_nn_i[:, 0]]) / cells.nn_len

    fxx = gfx_o * cells.cell_nn_tx

    gfy_o = (fy_pump_o[cells.cell_nn_i[:, 1]] - fy_pump_o[cells.cell_nn_i[:, 0]]) / cells.nn_len
    fyy = gfy_o * cells.cell_nn_ty

    divF_pump = fxx + fyy

    sim.rho_pump = sim.rho_pump + divF_pump * p.dt

    # ------------------------------------------------
    # calculate the total Nernst-Planck flux at each membrane for rho_channel factor:

    fx_chan, fy_chan = stb.nernst_planck_flux(sim.rho_channel, gcx_ch, gcy_ch, Ex, Ey, ux_mem, uy_mem, p.D_membrane,
        p.z_channel, sim.T, p)

    # map the flux to cell centres:
    fx_chan_o = np.dot(cells.M_sum_mems, fx_chan) / cells.num_mems
    fy_chan_o = np.dot(cells.M_sum_mems, fy_chan) / cells.num_mems

    # divergence of the total flux:
    gfx_o = (fx_chan_o[cells.cell_nn_i[:, 1]] - fx_chan_o[cells.cell_nn_i[:, 0]]) / cells.nn_len

    fxx = gfx_o * cells.cell_nn_tx

    gfy_o = (fy_chan_o[cells.cell_nn_i[:, 1]] - fy_chan_o[cells.cell_nn_i[:, 0]]) / cells.nn_len
    fyy = gfy_o * cells.cell_nn_ty

    divF_chan = fxx + fyy

    sim.rho_channel = sim.rho_channel + divF_chan * p.dt

    if p.sim_ECM is False:
        # average to the cell centre:
        sim.rho_pump_o = np.dot(cells.M_sum_mems, sim.rho_pump) / cells.num_mems
        sim.rho_channel_o = np.dot(cells.M_sum_mems, sim.rho_channel) / cells.num_mems

    # ------------------------------------------------
    # make sure nothing is non-zero:
    fix_inds = (sim.rho_pump < 0).nonzero()
    sim.rho_pump[fix_inds] = 0

    fix_inds2 = (sim.rho_channel < 0).nonzero()
    sim.rho_channel[fix_inds2] = 0