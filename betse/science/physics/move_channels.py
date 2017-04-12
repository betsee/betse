#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

import numpy as np


class MoveChannel(object):

    """
    Movement of ion pumps and channels to potentially create directional fluxes in individual cells.

    This is presently simulated by calculating the Nernst-Planck concentration flux of a weighting
    agent, rho, which moves under its own concentration gradient (a homogeneity restoring influence)
    as well as under the influence of the extracellular voltage gradient and fluid flows tangential
    to the membrane.

    """

    def __init__(self, sim, cells, p):

        # transfer rate for entities in the membrane:
        self.alpha = p.D_membrane/(p.tm*p.cell_height*cells.num_mems[cells.mem_to_cells])

        # initialize pumps:
        self.pgrad_x = np.zeros(sim.mdl)
        self.pgrad_y = np.zeros(sim.mdl)

        # initialize channels:
        self.cgrad_x = np.zeros(sim.mdl)
        self.cgrad_y = np.zeros(sim.mdl)

    def run(self, sim, cells, p):

        self.update_channel(sim, cells, p)
        self.update_pump(sim, cells, p)

    def update_pump(self, sim, cells, p):

        # calculate the equillibrium gradient vector pointing from maximum to minimum density on the membrane:
        ceqm_x = ((p.z_pump * p.q) / (p.kb * p.T)) * 1.0 * sim.J_cell_x[cells.mem_to_cells] * p.media_rho
        ceqm_y = ((p.z_pump * p.q) / (p.kb * p.T)) * 1.0 * sim.J_cell_y[cells.mem_to_cells] * p.media_rho

        # calculate update to the actual gradient concentration:
        self.pgrad_x = self.pgrad_x - self.alpha*(self.pgrad_x - ceqm_x) * p.dt
        self.pgrad_y = self.pgrad_y - self.alpha*(self.pgrad_y - ceqm_y) * p.dt

        sim.rho_pump = 1.0 + self.pgrad_x * (cells.nx_rads/2) + self.pgrad_y * (cells.ny_rads/2)

    def update_channel(self, sim, cells, p):
        # calculate the equillibrium gradient vector pointing from maximum to minimum density on the membrane:
        ceqm_x = ((p.z_channel * p.q) / (p.kb * p.T)) * 1.0 * sim.J_cell_x[cells.mem_to_cells] * p.media_rho
        ceqm_y = ((p.z_channel * p.q) / (p.kb * p.T)) * 1.0 * sim.J_cell_y[cells.mem_to_cells] * p.media_rho

        # calculate update to the actual gradient concentration:
        self.cgrad_x = self.cgrad_x - self.alpha * (self.cgrad_x - ceqm_x) * p.dt
        self.cgrad_y = self.cgrad_y - self.alpha * (self.cgrad_y - ceqm_y) * p.dt

        sim.rho_channel = 1.0 + self.cgrad_x * (cells.nx_rads/2) + self.cgrad_y * (cells.ny_rads/2)


    #------Old Method----------------------------------------------

    # # x and y components of membrane tangent unit vectors
    # tx = cells.mem_vects_flat[:, 4]
    # ty = cells.mem_vects_flat[:, 5]
    #
    # # map rho pump and rho channel to cell vertices:
    # pump_at_verts = np.dot(sim.rho_pump, cells.matrixMap2Verts)
    # channel_at_verts = np.dot(sim.rho_channel, cells.matrixMap2Verts)
    #
    # # get the gradient of rho concentration around each membrane:
    # grad_c_p =  np.dot(cells.gradMem, pump_at_verts)
    # grad_c_ch = np.dot(cells.gradMem, channel_at_verts)
    #
    # if p.sim_ECM is True:
    #
    #     Ex = sim.E_env_x.ravel()[cells.map_cell2ecm]
    #     Ey = sim.E_env_y.ravel()[cells.map_cell2ecm]
    #
    #     E_tang = (Ex[cells.mem_to_cells] * tx + Ey[cells.mem_to_cells] * ty)
    #
    # else:
    #
    #     Ex = p.media_rho*sim.J_cell_x[cells.mem_to_cells]
    #     Ey = p.media_rho*sim.J_cell_y[cells.mem_to_cells]
    #
    #
    #     E_tang = (Ex * tx + Ey * ty)
    #
    #
    # # calculate the total Nernst-Planck flux at each membrane for rho_pump factor:
    #
    # flux_pump = -p.D_membrane*grad_c_p  + \
    #             ((p.z_pump*p.D_membrane*p.F)/(p.R*p.T))*sim.rho_pump*E_tang
    #
    # flux_chan = -p.D_membrane * grad_c_ch  * sim.rho_channel + \
    #             ((p.z_channel * p.D_membrane * p.F) / (p.R * p.T)) * sim.rho_channel * E_tang
    #
    #
    # # divergence of the total flux via dfx/dx + dfy/dy:
    # # components of flux in x and y directions, mapped to cell verts:
    # # pumps:
    # gfx_p = np.dot(-flux_pump * tx, cells.matrixMap2Verts)
    # gfy_p = np.dot(-flux_pump * ty, cells.matrixMap2Verts)
    #
    # ddfx_p = np.dot(cells.gradMem, gfx_p) * tx
    # ddfy_p = np.dot(cells.gradMem, gfy_p) * ty
    #
    # divF_pump = ddfx_p + ddfy_p
    #
    # gfx_ch = np.dot(-flux_chan * tx, cells.matrixMap2Verts)
    # gfy_ch = np.dot(-flux_chan * ty, cells.matrixMap2Verts)
    #
    # ddfx_ch = np.dot(cells.gradMem, gfx_ch) * tx
    # ddfy_ch = np.dot(cells.gradMem, gfy_ch) * ty
    #
    # divF_chan = ddfx_ch + ddfy_ch
    #
    #
    #
    # sim.rho_pump = sim.rho_pump - divF_pump * p.dt
    # sim.rho_channel = sim.rho_channel - divF_chan * p.dt

    # ------------------------------------------------
    # make sure nothing is non-zero:
    # fix_inds = (sim.rho_pump < 0).nonzero()
    # sim.rho_pump[fix_inds] = 0
    #
    # fix_inds2 = (sim.rho_channel < 0).nonzero()
    # sim.rho_channel[fix_inds2] = 0

    def remove_data(self, targets_mem):

        # remove items from the mem concentration lists:
        pgrad_x = np.delete(self.pgrad_x, targets_mem)
        # reassign the new data vector to the object:
        self.pgrad_x = pgrad_x*1

        # remove items from the mem concentration lists:
        pgrad_y = np.delete(self.pgrad_y, targets_mem)
        # reassign the new data vector to the object:
        self.pgrad_y = pgrad_y * 1

        # remove items from the mem concentration lists:
        cgrad_x = np.delete(self.cgrad_x, targets_mem)
        # reassign the new data vector to the object:
        self.cgrad_x = cgrad_x*1

        # remove items from the mem concentration lists:
        cgrad_y = np.delete(self.cgrad_y, targets_mem)
        # reassign the new data vector to the object:
        self.cgrad_y = cgrad_y * 1

        # remove items from the mem concentration lists:
        alpha = np.delete(self.alpha, targets_mem)
        # reassign the new data vector to the object:
        self.alpha = alpha * 1

