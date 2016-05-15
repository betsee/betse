
#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

import numpy as np
from scipy import interpolate as interp
from scipy.ndimage.filters import gaussian_filter
from betse.exceptions import BetseExceptionSimulation
from betse.science import sim_toolbox as stb
from betse.util.io.log import logs


def getDeformation(sim, cells, t, p):
    """
    Calculates the deformation of the cell cluster under the action
    of intracellular forces and pressures, assuming steady-state
    (slow) changes.

    The method assumes that material is entirely incompressible.
    First, the equation of linear elastic motion is used to calculate
    deformation assuming full compressibility. Then, the divergence
    is calculated, an internal reaction pressure is calculated from
    the divergence, and the gradient of the reaction pressure is subtracted
    from the initial solution to create a divergence-free deformation field.

    """

    # Determine action forces ------------------------------------------------
    # # FIXME: Update this when focusing on forces at a later time
    # # assume action forces for deformation are pseudo forces, resulting from cell galvanotropism,
    # # but are in relation to the global electric field (proportional to environmental current)
    # F_x = p.galvanotropism * sim.J_env_x.ravel()
    # F_y = p.galvanotropism * sim.J_env_y.ravel()
    #
    # # now interpolate these forces at the cell centres:
    # # interpolate charge from environmental grid to the ecm_mids:
    # F_cell_x = interp.griddata((cells.xypts[:, 0], cells.xypts[:, 1]),
    #     F_x, (cells.cell_centres[:, 0], cells.cell_centres[:, 1]), method='nearest', fill_value=0)
    #
    # F_cell_y = interp.griddata((cells.xypts[:, 0], cells.xypts[:, 1]),
    #     F_y, (cells.cell_centres[:, 0], cells.cell_centres[:, 1]), method='nearest', fill_value=0)

    #---------------------------------------------------------------------------------

    # body forces from hydrostatic pressure
    F_hydro_x = sim.F_hydro_x
    F_hydro_y = sim.F_hydro_y

    # first determine body force components due to electrostatics, if desired:
    if p.deform_electro is True:
        F_electro_x = sim.F_electro_x
        F_electro_y = sim.F_electro_y

    else:

        F_electro_x = np.zeros(len(cells.cell_i))
        F_electro_y = np.zeros(len(cells.cell_i))

    # Take the total component of pressure from all contributions:
    F_cell_x = F_electro_x + F_hydro_x
    F_cell_y = F_electro_y + F_hydro_y


    # integrate the forces, as is mandated by finite volume methods:
    F_cell_x = cells.integrator(F_cell_x)
    F_cell_y = cells.integrator(F_cell_y)

    # --calculate displacement field for incompressible medium------------------------------------------------

    # calculate the initial displacement field (not divergence free!) for the forces using the linear elasticity
    # equation:

    if p.fixed_cluster_bound is True:

        u_x_o = np.dot(cells.lapGJinv, -(1 / p.lame_mu) * (F_cell_x))
        u_y_o = np.dot(cells.lapGJinv, -(1 / p.lame_mu) * (F_cell_y))

        # enforce boundary conditions on u:
        if p.fixed_cluster_bound is True:
            u_x_o[cells.bflags_cells] = 0
            u_y_o[cells.bflags_cells] = 0
            u_x_o[cells.nn_bound] = 0
            u_y_o[cells.nn_bound] = 0

    else:

        u_x_o = np.dot(cells.lapGJ_P_inv, -(1 / p.lame_mu) * (F_cell_x))
        u_y_o = np.dot(cells.lapGJ_P_inv, -(1 / p.lame_mu) * (F_cell_y))

        # first interpolate displacement field at membrane midpoints:
    ux_mem = interp.griddata((cells.cell_centres[:, 0], cells.cell_centres[:, 1]), u_x_o,
        (cells.mem_mids_flat[:, 0], cells.mem_mids_flat[:, 1]), fill_value=0)

    uy_mem = interp.griddata((cells.cell_centres[:, 0], cells.cell_centres[:, 1]), u_y_o,
        (cells.mem_mids_flat[:, 0], cells.mem_mids_flat[:, 1]), fill_value=0)

    # get the component of the displacement field normal to the membranes:
    u_n = ux_mem * cells.mem_vects_flat[:, 2] + uy_mem * cells.mem_vects_flat[:, 3]

    # calculate divergence as the sum of this vector x each surface area, divided by cell volume:
    div_u = (np.dot(cells.M_sum_mems, u_n * cells.mem_sa) / cells.cell_vol)

    # calculate the reaction pressure required to counter-balance the flow field:

    # if p.fixed_cluster_bound is True:

    P_react = np.dot(cells.lapGJ_P_inv, 2 * div_u)  # trial and error suggest this needs to be 2x the divergence?

    # calculate its gradient:
    gradP_react = (P_react[cells.cell_nn_i[:, 1]] - P_react[cells.cell_nn_i[:, 0]]) / (cells.nn_len)

    gP_x = gradP_react * cells.cell_nn_tx
    gP_y = gradP_react * cells.cell_nn_ty

    # average the components of the reaction force field at cell centres and get boundary values:
    gPx_cell = np.dot(cells.M_sum_mems, gP_x) / cells.num_mems
    gPy_cell = np.dot(cells.M_sum_mems, gP_y) / cells.num_mems

    # calculate the displacement of cell centres under the applied force under incompressible conditions:
    sim.d_cells_x = u_x_o - gPx_cell
    sim.d_cells_y = u_y_o - gPy_cell

    # # enforce boundary conditions:
    if p.fixed_cluster_bound is True:
        sim.d_cells_x[cells.bflags_cells] = 0
        sim.d_cells_y[cells.bflags_cells] = 0
        # self.d_cells_x[cells.nn_bound] = 0
        # self.d_cells_y[cells.nn_bound] = 0

def timeDeform(sim, cells, t, p):
    """
    Calculates the deformation of the cell cluster under the action
    of intracellular pressure, considering the full time-dependent
    linear elasticity equation for an incompressible medium.

    The solution method for this equation is similar to the
    steady-state method of deformation(). First the displacement
    field is calculated assuming compressibility,
    a reaction pressure is calculated from the divergence of the
    initial field, and the gradient of the internal pressure is
    subtracted from the initial field to produce a divergence
    free solution.

    This method is working much better than the timeDeform_o()
    so is presently in active use.

    """

    # Check for the adequacy of the time step:
    step_check = (p.dt / (2 * p.rc)) * np.sqrt(p.lame_mu / 1000)

    if step_check > 1.0:
        new_ts = (0.9 * 2 * p.rc) / (np.sqrt(p.lame_mu / 1000))

        raise BetseExceptionSimulation(
            'Time dependent deformation is tricky business, requiring a small time step! '
            'The time step you are using is too large to bother going further with. '
            'Please set your time step to ' + str(new_ts) + ' and try again.')

    k_const = (p.dt ** 2) * (p.lame_mu / 1000)

    # # Determine action forces ------------------------------------------------

    # body force from hydrostatic pressure:
    F_hydro_x = sim.F_hydro_x
    F_hydro_y = sim.F_hydro_y

    # first determine body force components due to electrostatics, if desired:
    if p.deform_electro is True:
        F_electro_x = sim.F_electro_x
        F_electro_y = sim.F_electro_y

    else:

        F_electro_x = np.zeros(len(cells.cell_i))
        F_electro_y = np.zeros(len(cells.cell_i))

    # Take the total component of pressure from all contributions:
    F_cell_x = F_electro_x + F_hydro_x
    F_cell_y = F_electro_y + F_hydro_y

    # integrate the forces:
    F_cell_x = cells.integrator(F_cell_x)
    F_cell_y = cells.integrator(F_cell_y)

    # -------------------------------------------------------------------------------------------------

    sim.dx_time.append(sim.d_cells_x[:])  # append the solution to the time-save vector
    sim.dy_time.append(sim.d_cells_y[:])

    # Initial value solution--------------------------------------------------------------------------------
    if t == 0.0:

        wave_speed = np.sqrt(p.lame_mu / 1000)
        wave_speed = np.float(wave_speed)
        wave_speed = np.round(wave_speed, 2)

        logs.log_info(
            'Your wave speed is approximately: ' +
            str(wave_speed) + ' m/s '
        )

        logs.log_info('Try a world size of at least: ' + str(round((5 / 3) * (wave_speed / 500) * 1e6))
                      + ' um for resonance.')

        if p.fixed_cluster_bound is True:

            sim.d_cells_x = k_const * np.dot(cells.lapGJ, sim.dx_time[-1]) + (k_const / p.lame_mu) * F_cell_x + \
                            sim.dx_time[-1]
            sim.d_cells_y = k_const * np.dot(cells.lapGJ, sim.dy_time[-1]) + (k_const / p.lame_mu) * F_cell_y + \
                            sim.dy_time[-1]

        else:

            sim.d_cells_x = k_const * np.dot(cells.lapGJ_P, sim.dx_time[-1]) + (k_const / p.lame_mu) * F_cell_x + \
                            sim.dx_time[-1]
            sim.d_cells_y = k_const * np.dot(cells.lapGJ_P, sim.dy_time[-1]) + (k_const / p.lame_mu) * F_cell_y + \
                            sim.dy_time[-1]



    elif t > 0.0:

        # do the non-initial value, standard solution iteration:

        # calculate the velocity for viscous damping:
        d_ux_dt = (sim.dx_time[-1] - sim.dx_time[-2]) / (p.dt)
        d_uy_dt = (sim.dy_time[-1] - sim.dy_time[-2]) / (p.dt)

        gamma = ((p.dt ** 2) * (p.mu_tissue * p.lame_mu)) / (1000 * (2 * p.rc))

        if p.fixed_cluster_bound is True:

            sim.d_cells_x = k_const * np.dot(cells.lapGJ, sim.dx_time[-1]) - gamma * d_ux_dt + \
                             (k_const / p.lame_mu) * F_cell_x + 2 * sim.dx_time[-1] - sim.dx_time[-2]

            sim.d_cells_y = k_const * np.dot(cells.lapGJ, sim.dy_time[-1]) - gamma * d_uy_dt + \
                             (k_const / p.lame_mu) * F_cell_y + 2 * sim.dy_time[-1] - sim.dy_time[-2]

        else:

            sim.d_cells_x = k_const * np.dot(cells.lapGJ_P, sim.dx_time[-1]) - gamma * d_ux_dt + \
                             (k_const / p.lame_mu) * F_cell_x + 2 * sim.dx_time[-1] - sim.dx_time[-2]

            sim.d_cells_y = k_const * np.dot(cells.lapGJ_P, sim.dy_time[-1]) - gamma * d_uy_dt + \
                             (k_const / p.lame_mu) * F_cell_y + 2 * sim.dy_time[-1] - sim.dy_time[-2]


            # calculate divergence of u  -----------------------------------------------------------------------

            # first interpolate displacement field at membrane midpoints:
    ux_mem = interp.griddata((cells.cell_centres[:, 0], cells.cell_centres[:, 1]), sim.d_cells_x,
        (cells.mem_mids_flat[:, 0], cells.mem_mids_flat[:, 1]), fill_value=0)

    uy_mem = interp.griddata((cells.cell_centres[:, 0], cells.cell_centres[:, 1]), sim.d_cells_y,
        (cells.mem_mids_flat[:, 0], cells.mem_mids_flat[:, 1]), fill_value=0)

    # get the component of the displacement field normal to the membranes:
    u_n = ux_mem * cells.mem_vects_flat[:, 2] + uy_mem * cells.mem_vects_flat[:, 3]

    # calculate divergence as the sum of this vector x each surface area, divided by cell volume:
    div_u = (np.dot(cells.M_sum_mems, u_n * cells.mem_sa) / cells.cell_vol)

    # calculate the reaction pressure required to counter-balance the flow field:

    P_react = np.dot(cells.lapGJ_P_inv, 2 * div_u)  # FIXME why a two?

    # self.P_cells = (p.lame_mu/k_const)*P_react[:]

    # calculate its gradient:
    gradP_react = (P_react[cells.cell_nn_i[:, 1]] - P_react[cells.cell_nn_i[:, 0]]) / (cells.nn_len)

    gP_x = gradP_react * cells.cell_nn_tx
    gP_y = gradP_react * cells.cell_nn_ty

    # average the components of the reaction force field at cell centres and get boundary values:
    gPx_cell = np.dot(cells.M_sum_mems, gP_x) / cells.num_mems
    gPy_cell = np.dot(cells.M_sum_mems, gP_y) / cells.num_mems

    # calculate the displacement of cell centres under the applied force under incompressible conditions:
    sim.d_cells_x = sim.d_cells_x - gPx_cell
    sim.d_cells_y = sim.d_cells_y - gPy_cell

    if p.fixed_cluster_bound is True:  # enforce zero displacement boundary condition:

        sim.d_cells_x[cells.bflags_cells] = 0
        sim.d_cells_y[cells.bflags_cells] = 0

        sim.d_cells_x[cells.nn_bound] = 0
        sim.d_cells_y[cells.nn_bound] = 0

    # check the displacement for NANs:
    stb.check_v(sim.d_cells_x)

def implement_deform_timestep(sim, cells, t, p):
    # Map individual cell deformations to their membranes. In this case,
    # this is better than interpolation.
    ux_at_mem = sim.d_cells_x[cells.mem_to_cells]
    uy_at_mem = sim.d_cells_y[cells.mem_to_cells]

    ux_at_ecm = np.dot(cells.M_sum_mem_to_ecm, ux_at_mem)
    uy_at_ecm = np.dot(cells.M_sum_mem_to_ecm, uy_at_mem)

    # get new ecm verts:
    new_ecm_verts_x = sim.ecm_verts_unique_to[:, 0] + np.dot(cells.deforM, ux_at_ecm)
    new_ecm_verts_y = sim.ecm_verts_unique_to[:, 1] + np.dot(cells.deforM, uy_at_ecm)

    ecm_new = np.column_stack((new_ecm_verts_x, new_ecm_verts_y))

    # set the voronoi points originally tagged to the ecm to the value of these new points
    cells.voronoi_grid[cells.map_voronoi2ecm] = ecm_new[:]

    # recreate ecm_verts_unique:
    cells.ecm_verts_unique = ecm_new[:]

    # Repackage ecm verts so that the World module can do its magic:
    ecm_new_flat = ecm_new[cells.ecmInds]  # first expand it to a flattened form (include duplictes)

    # Repackage the structure to include individual cell data.
    cells.ecm_verts = []  # null the original ecm verts data structure...

    # Convert region to a numpy array so it can be sorted.
    for i in range(0, len(cells.cell_to_mems)):
        ecm_nest = ecm_new_flat[cells.cell_to_mems[i]]
        ecm_nest = np.asarray(ecm_nest)
        cells.ecm_verts.append(ecm_nest)

    # Voila! Deformed ecm_verts!
    cells.ecm_verts = np.asarray(cells.ecm_verts)
    cells.deformWorld(p)

    sim.cell_centres_time.append(cells.cell_centres[:])
    sim.mem_mids_time.append(cells.mem_mids_flat[:])
    sim.maskM_time.append(cells.maskM[:])
    sim.mem_edges_time.append(cells.mem_edges_flat[:])
    sim.cell_verts_time.append(cells.cell_verts[:])

    sim.dx_cell_time.append(sim.d_cells_x[:])
    sim.dy_cell_time.append(sim.d_cells_y[:])