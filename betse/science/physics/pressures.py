#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.


import numpy as np
from scipy.ndimage.filters import gaussian_filter


def osmotic_P(sim, cells, p):
    # initialize osmotic pressures in cells and env

    sim.osmo_P_cell = np.zeros(len(sim.cc_cells[0]))
    sim.osmo_P_env = np.zeros(len(sim.cc_env[0]))

    # calculate osmotic pressure in cells based on total molarity:
    for c_ion in sim.cc_cells:
        sim.osmo_P_cell = c_ion * p.R * sim.T + sim.osmo_P_cell

    # calculate osmotic pressure in environment based on total molarity:

    for c_ion_env in sim.cc_env:
        sim.osmo_P_env = c_ion_env * p.R * sim.T + sim.osmo_P_env

    if p.sim_ECM is False:

        op_env = np.dot(cells.M_sum_mems, sim.osmo_P_env)

        sim.osmo_P_delta = sim.osmo_P_cell - op_env

    else:
        # smooth out the environmental osmotic pressure:
        sim.osmo_P_env = gaussian_filter(sim.osmo_P_env.reshape(cells.X.shape),p.smooth_level).ravel()

        sim.osmo_P_delta = sim.osmo_P_cell - sim.osmo_P_env[cells.map_cell2ecm]

    # Calculate the transmembrane flow of water due to osmotic pressure.
    # High, positive osmotic pressure leads to water flow into the cell. Existing pressure in the cell
    # resists the degree of osmotic influx. The effect also depends on aquaporin fraction in membrane:
    u_osmo_o = (sim.osmo_P_delta - sim.P_cells) * (p.aquaporins / p.mu_water)

    # map osmotic influx to membranes:
    u_osmo = u_osmo_o[cells.mem_to_cells]

    # u_osmo = interp.griddata((cells.cell_centres[:, 0], cells.cell_centres[:, 1]), sim.d_cells_x,
    #     (cells.mem_mids_flat[:, 0], cells.mem_mids_flat[:, 1]), fill_value=0)
    #
    # uy_mem = interp.griddata((cells.cell_centres[:, 0], cells.cell_centres[:, 1]), sim.d_cells_y,
    #     (cells.mem_mids_flat[:, 0], cells.mem_mids_flat[:, 1]), fill_value=0)


    # calculate the flow due to net mass flux into the cell due to selective active/passive ion transport:
    u_mass_flux = (1 / p.rho) * get_mass_flux(sim, cells, p)

    sim.u_net = u_osmo + u_mass_flux

    # obtain the divergence of the flow:
    div_u_osmo = np.dot(cells.M_sum_mems, sim.u_net * cells.mem_sa) / cells.cell_vol

    # the divergence is sequential/accumulative with each time step [???]:
    sim.div_u_osmo = div_u_osmo

    # pressure developing in the cell depends on how much the volume can change:
    P_react = (1 - (1 / p.youngMod)) * sim.div_u_osmo

    # the inflow of mass adds to base pressure in cells
    # (this format is used to avoid conflict with pressure channels):
    if p.run_sim is True:
        sim.P_base = sim.P_base + P_react[:]

    else:
        sim.P_cells = sim.P_cells + P_react[:]

    # ------------------------------------------------------------------------------------------------------------
    # actual volume change depends on the mechanical properties (young's modulus) of the
    # tissue and the ability for cells to deform:

    sim.delta_vol = (1 / p.youngMod) * div_u_osmo * p.dt

    # update concentrations and volume in the cell:
    vo = cells.cell_vol[:]

    v1 = (1 + sim.delta_vol) * vo

    vol_ratio = (vo/v1)

    # sim.cc_mems = sim.cc_mems * vol_ratio[cells.mem_to_cells]
    sim.cc_cells = sim.cc_cells * vol_ratio

    # reassign cell volume:
    cells.cell_vol = v1[:]

    # reassign mem volume:
    cells.mem_vol = cells.mem_vol/vol_ratio[cells.mem_to_cells]

    # if p.sim_ECM is True:
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

def get_mass_flux(sim, cells, p):
    """
    Sum up individual trans-membrane and
    trans-gap junction ion fluxes to obtain the
    net flow of mass into a cell.

    Assumes that each ion travels with a hydration
    shell of 6 water molecules.

    """

    # calculate mass flux across cell membranes:
    mass_flux = np.zeros(len(cells.mem_i))

    for flux_array, mm in zip(sim.fluxes_mem, sim.molar_mass):
        m_flx = flux_array * (mm + 6 * 18.01e-3)  # flux x molar mass of ion x 6 water molecules at 18e-3 kg/mol

        mass_flux = mass_flux + m_flx

    return mass_flux

    # # total mass change in cell
    # mass_change = self.mass_flux*p.dt*cells.mem_sa
    # # sum the change over the membranes to get the total mass change of salts:
    # self.delta_m_salts = np.dot(cells.M_sum_mems,mass_change)