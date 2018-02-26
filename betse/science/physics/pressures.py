#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.


import numpy as np
from scipy.ndimage.filters import gaussian_filter


def osmotic_P(sim, cells, p):

    # initialize osmotic pressures in cells and env
    sim.osmo_P_cell = np.sum(p.R * sim.T * sim.cc_cells, axis=0)
    sim.osmo_P_env = np.sum(p.R * sim.T * sim.cc_env, axis=0)

    if p.is_ecm is False:

        op_env = np.dot(cells.M_sum_mems, sim.osmo_P_env)/cells.num_mems

        sim.osmo_P_delta = op_env - sim.osmo_P_cell

    else:

        sim.osmo_P_delta = sim.osmo_P_env[cells.map_cell2ecm] - sim.osmo_P_cell


    # Calculate the transmembrane flow of water due to osmotic pressure.
    # High, positive osmotic pressure leads to water flow into the cell. Existing pressure in the cell
    # resists the degree of osmotic influx. The effect also depends on aquaporin fraction in membrane:
    # sim.u_osmo = (sim.osmo_P_delta - sim.P_cells) * (p.aquaporins / p.mu_water)*cells.cell_sa
    sim.u_osmo = (sim.osmo_P_delta - sim.P_cells)*(p.aquaporins/p.mu_water)*(3e-10**2)*(1/p.tm)

    # obtain the divergence of the flow -- this is a strain rate:
    sim.div_u_osmo = sim.u_osmo*(cells.cell_sa/cells.cell_vol)

    # calculate pressure in whole network resulting from divergence :
    sim.PP = np.dot(cells.lapGJ_P_inv, -sim.div_u_osmo * p.rho * p.dt)

    # ------------------------------------------------------------------------------------------------------------
    # actual volume change is amount of flow over the cell surface area per unit time:
    sim.delta_vol = p.dt*sim.u_osmo*cells.cell_sa

    vol_ratio = cells.cell_vol / (cells.cell_vol - sim.delta_vol)

    # new concentrations in cells from C1*V1 = C2*V2 expression:
    sim.cc_cells = sim.cc_cells * vol_ratio

    # reassign cell volume:
    cells.cell_vol = cells.cell_vol - sim.delta_vol

    # reassign mem volume:
    cells.mem_vol = cells.mem_vol/vol_ratio[cells.mem_to_cells]

    # if p.is_ecm is True:
    #     vo_ecm = cells.ecm_vol[cells.map_cell2ecm]
    #     v1_ecm = (1 - sim.delta_vol) * vo_ecm
    #
    #     for i, cc_array in enumerate(sim.cc_env):
    #         sim.cc_env[i][cells.map_cell2ecm] = cc_array[cells.map_cell2ecm] * (vo_ecm / v1_ecm)
    #
    #     cells.ecm_vol[cells.map_cell2ecm] = v1_ecm
    #
    # if p.voltage_dye is True:
    #     sim.cDye_cell = sim.cDye_cell * vol_ratio
    #
    # if p.scheduled_options['IP3'] != 0 or p.Ca_dyn is True:
    #     sim.cIP3 = sim.cIP3 * (vo / v1)

# def get_mass_flux(sim, cells, p):
#     """
#     Sum up individual trans-membrane and
#     trans-gap junction ion fluxes to obtain the
#     net flow of mass into a cell.
#
#     Assumes that each ion travels with a hydration
#     shell of 6 water molecules.
#
#     """
#
#     # calculate mass flux across cell membranes:
#     mass_flux = np.zeros(len(cells.mem_i))
#
#     for flux_array, mm in zip(sim.fluxes_mem, sim.molar_mass):
#         m_flx = flux_array * (mm + 6 * 18.01e-3)  # flux x molar mass of ion x 6 water molecules at 18e-3 kg/mol
#
#         mass_flux = mass_flux + m_flx
#
#     return mass_flux

    # # total mass change in cell
    # mass_change = self.mass_flux*p.dt*cells.mem_sa
    # # sum the change over the membranes to get the total mass change of salts:
    # self.delta_m_salts = np.dot(cells.M_sum_mems,mass_change)