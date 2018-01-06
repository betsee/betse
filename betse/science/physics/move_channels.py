#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

import numpy as np
from betse.science import sim_toolbox as stb


class MoveChannel(object):

    """
    Movement of ion pumps and channels to potentially create directional fluxes in individual cells.

    This is presently simulated by calculating the Nernst-Planck concentration flux of a weighting
    agent, rho, which moves under its own concentration gradient (a homogeneity restoring influence)
    as well as under the influence of the extracellular voltage gradient and fluid flows tangential
    to the membrane.

    """

    def __init__(self, sim, cells, p):

        sim.rho_pump = np.ones(sim.mdl)
        sim.rho_channel = np.ones(sim.mdl)

        # transfer rate for entities in the membrane:
        # self.alpha = p.D_membrane/(p.tm*p.cell_height*cells.num_mems[cells.mem_to_cells])
        #
        # # initialize pumps:
        # self.pgrad_x = np.zeros(sim.mdl)
        # self.pgrad_y = np.zeros(sim.mdl)
        #
        # # initialize channels:
        # self.cgrad_x = np.zeros(sim.mdl)
        # self.cgrad_y = np.zeros(sim.mdl)

    def run(self, sim, cells, p):

        self.update_channel(sim, cells, p)
        self.update_pump(sim, cells, p)

    def update_pump(self, sim, cells, p):

        cav = 1.0  # concentration at cell centre
        cpi = sim.rho_pump  # concentration at membrane
        z = p.z_pump  # charge of ion
        Do = p.D_membrane  # diffusion constant of ion

        cap = (cav + cpi) / 2  # concentration at midpoint between cell centre and membrane
        cgp = (cpi - cav) / cells.R_rads  # concentration gradients

        # Nernst-Planck equation to calculate flux:
        cfluxpo = -Do * cgp + ((Do * p.q * cap * z) / (p.kb * sim.T)) * sim.Ec

        # as no net mass must leave this intracellular movement, make the flux divergence-free:
        cfluxp = stb.single_cell_div_free(cfluxpo, cells)

        # calculate the actual concentration at membranes by unpacking to concentration vectors:
        sim.rho_pump = cpi + cfluxp * (cells.mem_sa / cells.mem_vol) * p.dt

        # calculate the equillibrium gradient vector pointing from maximum to minimum density on the membrane:
        # ceqm_x = ((p.z_pump * p.q) / (p.kb * p.T)) * 1.0 * sim.J_cell_x[cells.mem_to_cells] * (1/sim.sigma)
        # ceqm_y = ((p.z_pump * p.q) / (p.kb * p.T)) * 1.0 * sim.J_cell_y[cells.mem_to_cells] * (1/sim.sigma)
        #
        # # calculate update to the actual gradient concentration:
        # self.pgrad_x = self.pgrad_x - self.alpha*(self.pgrad_x - ceqm_x) * p.dt
        # self.pgrad_y = self.pgrad_y - self.alpha*(self.pgrad_y - ceqm_y) * p.dt
        #
        # sim.rho_pump = 1.0 + self.pgrad_x * (cells.nx_rads/2) + self.pgrad_y * (cells.ny_rads/2)

    def update_channel(self, sim, cells, p):

        cav = 1.0  # concentration at cell centre
        cpi = sim.rho_channel  # concentration at membrane
        z = p.z_channel  # charge of ion
        Do = p.D_membrane  # diffusion constant of ion

        cap = (cav + cpi) / 2  # concentration at midpoint between cell centre and membrane
        cgp = (cpi - cav) / cells.R_rads  # concentration gradients

        # Nernst-Planck equation to calculate flux:
        cfluxco = -Do * cgp + ((Do * p.q * cap * z) / (p.kb * sim.T)) * sim.Ec

        # as no net mass must leave this intracellular movement, make the flux divergence-free:
        cfluxc = stb.single_cell_div_free(cfluxco, cells)

        # calculate the actual concentration at membranes by unpacking to concentration vectors:
        sim.rho_channel = cpi + cfluxc * (cells.mem_sa / cells.mem_vol) * p.dt

        # calculate the equillibrium gradient vector pointing from maximum to minimum density on the membrane:
        # ceqm_x = ((p.z_channel * p.q) / (p.kb * p.T)) * 1.0 * sim.J_cell_x[cells.mem_to_cells] * (1/sim.sigma)
        # ceqm_y = ((p.z_channel * p.q) / (p.kb * p.T)) * 1.0 * sim.J_cell_y[cells.mem_to_cells] * (1/sim.sigma)
        #
        # # calculate update to the actual gradient concentration:
        # self.cgrad_x = self.cgrad_x - self.alpha * (self.cgrad_x - ceqm_x) * p.dt
        # self.cgrad_y = self.cgrad_y - self.alpha * (self.cgrad_y - ceqm_y) * p.dt
        #
        # sim.rho_channel = 1.0 + self.cgrad_x * (cells.nx_rads/2) + self.cgrad_y * (cells.ny_rads/2)

    def remove_data(self, targets_cell):

        pass

        # # remove items from the mem concentration lists:
        # rpo2 = np.delete(self.rpo, targets_cell)
        # # reassign the new data vector to the object:
        # self.rpo = rpo2*1
        #
        # # remove items from the mem concentration lists:
        # rco2 = np.delete(self.rco, targets_cell)
        # # reassign the new data vector to the object:
        # self.rco = rco2*1


