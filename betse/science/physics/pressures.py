#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.


import numpy as np
from betse.science import finitediff as fd


def electro_F(sim, cells, p):
    """
    Calculates electrostatic body force between gap junctions
     of a networked cell collective.

    """

    # # map charge in cell to the membrane, averaging between two cells:
    # Q_mem = (sim.rho_cells[cells.cell_nn_i[:, 1]] + sim.rho_cells[cells.cell_nn_i[:, 0]]) / 2

    # calculate force at each membrane:
    sim.F_gj_x = sim.rho_cells[cells.mem_to_cells] * sim.E_gj_x
    sim.F_gj_y = sim.rho_cells[cells.mem_to_cells] * sim.E_gj_y

    # calculate a shear electrostatic body force at the cell centre:
    sim.F_electro_x = np.dot(cells.M_sum_mems, sim.F_gj_x) / cells.num_mems
    sim.F_electro_y = np.dot(cells.M_sum_mems, sim.F_gj_y) / cells.num_mems

    sim.F_electro = np.sqrt(sim.F_electro_x ** 2 + sim.F_electro_y ** 2)

    # define this in terms of pressure (force per unit area)
    P_x = (sim.F_electro_x * cells.cell_vol) / cells.cell_sa
    P_y = (sim.F_electro_y * cells.cell_vol) / cells.cell_sa

    sim.P_electro = np.sqrt(P_x ** 2 + P_y ** 2)

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

        sim.osmo_P_delta = sim.osmo_P_cell - sim.osmo_P_env

    else:
        # smooth out the environmental osmotic pressure:
        sim.osmo_P_env = sim.osmo_P_env.reshape(cells.X.shape)
        sim.osmo_P_env = fd.integrator(sim.osmo_P_env)
        sim.osmo_P_env = sim.osmo_P_env.ravel()

        sim.osmo_P_delta = sim.osmo_P_cell[cells.mem_to_cells] - sim.osmo_P_env[cells.map_mem2ecm]

        # average the pressure to individual cells
        sim.osmo_P_delta = np.dot(cells.M_sum_mems, sim.osmo_P_delta) / cells.num_mems

    # Calculate the transmembrane flow of water due to osmotic pressure.
    # High, positive osmotic pressure leads to water flow into the cell. Existing pressure in the cell
    # resists the degree of osmotic influx. The effect also depends on aquaporin fraction in membrane:
    u_osmo_o = (sim.osmo_P_delta - sim.P_cells) * (p.aquaporins / p.mu_water)

    # map osmotic influx to membranes:
    u_osmo = u_osmo_o[cells.mem_to_cells]

    # calculate the flow due to net mass flux into the cell due to selective active/passive ion transport:
    u_mass_flux = (1 / p.rho) * sim.get_mass_flux(cells, p)

    u_net = u_osmo + u_mass_flux

    # obtain the divergence of the flow in a timestep, which yields the fractional volume change:
    div_u_osmo = p.dt * np.dot(cells.M_sum_mems, u_net * cells.mem_sa) / cells.cell_vol

    # pressure developing in the cell depends on how much the volume can change:
    P_react = (1 - (1 / p.youngMod)) * div_u_osmo

    # the inflow of mass adds to base pressure in cells
    # (this format is used to avoid conflict with pressure channels):
    if p.run_sim is True:
        sim.P_base = sim.P_base + P_react[:]

    else:
        sim.P_cells = sim.P_cells + P_react[:]

    # ------------------------------------------------------------------------------------------------------------
    # actual volume change depends on the mechanical properties (young's modulus) of the
    # tissue:

    sim.delta_vol = (1 / p.youngMod) * div_u_osmo

    # update concentrations and volume in the cell:
    vo = cells.cell_vol[:]

    v1 = (1 + sim.delta_vol) * vo

    sim.cc_cells = sim.cc_cells * (vo / v1)

    # reassign cell volume:
    cells.cell_vol = v1[:]

    # if p.sim_ECM is True:   FIXME environmental changes need to be updated as well....
    #     vo_ecm = cells.ecm_vol[cells.map_cell2ecm]
    #     v1_ecm = (1 - sim.delta_vol) * vo_ecm
    #
    #     for i, cc_array in enumerate(sim.cc_env):
    #         sim.cc_env[i][cells.map_cell2ecm] = cc_array[cells.map_cell2ecm] * (vo_ecm / v1_ecm)
    #
    #     cells.ecm_vol[cells.map_cell2ecm] = v1_ecm

    if p.voltage_dye is True:
        sim.cDye_cell = sim.cDye_cell * (vo / v1)

    if p.scheduled_options['IP3'] != 0 or p.Ca_dyn is True:
        sim.cIP3 = sim.cIP3 * (vo / v1)

def getHydroF(sim, cells, p):
    # ----Calculate body forces due to hydrostatic pressure gradients---------------------------------------------

    # determine body force due to hydrostatic pressure gradient between cells:

    gPcells = -(sim.P_cells[cells.cell_nn_i[:, 1]] - sim.P_cells[cells.cell_nn_i[:, 0]]) / cells.nn_len

    sim.F_hydro_x_gj = gPcells * cells.cell_nn_tx
    sim.F_hydro_y_gj = gPcells * cells.cell_nn_ty

    # calculate a shear electrostatic body force at the cell centre:
    sim.F_hydro_x = np.dot(cells.M_sum_mems, sim.F_hydro_x_gj) / cells.num_mems
    sim.F_hydro_y = np.dot(cells.M_sum_mems, sim.F_hydro_y_gj) / cells.num_mems

    sim.F_hydro = np.sqrt(sim.F_hydro_x ** 2 + sim.F_hydro_y ** 2)

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