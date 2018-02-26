#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

import numpy as np
from scipy.ndimage.filters import gaussian_filter
from betse.science import sim_toolbox as stb
from betse.science.math import finitediff as fd


def get_current(sim, cells, p):

    # calculate the net charge in cells:
    sim.rho_cells = np.dot(sim.zs * p.F, sim.cc_cells) + sim.extra_rho_cells

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

    # conductivity of cells is assumed to be 0.1x lower than that of free environment:
    sim.sigma_cell = np.asarray([((z ** 2) * p.q * p.F * cc * D * 0.1) / (p.kb * p.T) for (z, cc, D) in
                                 zip(sim.zs, sim.cc_cells, sim.D_free)]).mean(axis=0)

    # calculate electric field in cells using net intracellular current and cytosol conductivity:
    sim.E_cell_x = sim.J_cell_x*(1/(sim.sigma_cell))
    sim.E_cell_y = sim.J_cell_y*(1/(sim.sigma_cell))

    sim.Emc = (sim.E_cell_x[cells.mem_to_cells] * cells.mem_vects_flat[:, 2] +
               sim.E_cell_y[cells.mem_to_cells] * cells.mem_vects_flat[:, 3])


    # Current in the environment --------------------------------------------------------------------------------------
    if p.is_ecm is True:

        # calculate current densities in the environment from ion fluxes:
        J_env_x_o = np.dot(p.F*sim.zs, sim.fluxes_env_x) + sim.extra_Jenv_x
        J_env_y_o = np.dot(p.F*sim.zs, sim.fluxes_env_y) + sim.extra_Jenv_y

        # reshape the matrix:
        J_env_x_o = J_env_x_o.reshape(cells.X.shape)
        J_env_y_o = J_env_y_o.reshape(cells.X.shape)

        #---------------------------------------------

        #Helmholtz-Hodge decomposition to obtain divergence-free projection of actual currents (zero n_hat at boundary):
        _, sim.J_env_x, sim.J_env_y, BB, Jbx, Jby = stb.HH_Decomp(J_env_x_o, J_env_y_o, cells)

        # free current density at each membrane, "smoothed" using Helmholtz-Hodge (and without transmem contribution):
        sim.Jtx = sim.J_env_x + Jbx
        sim.Jty = sim.J_env_y + Jby

        # calculate the net charge in the environment:
        sim.rho_env = np.dot(sim.zs * p.F, sim.cc_env) + sim.extra_rho_env
        # sim.rho_env = fd.integrator(sim.rho_env.reshape(cells.X.shape), sharp = 0.5).ravel()

        # Method 6-------------------------------------------------------------------------------------------------
        # Calculate voltage from charge in environment; electric field is divergence-free assuming Laplace Eqn holds

        v_env = (1/(2*p.cm))*(sim.rho_env*cells.delta)

        v_env = v_env.reshape(cells.X.shape)

        v_env = fd.integrator(v_env, 0.5)

        # assign to environmental voltage array:
        sim.v_env = 1*v_env.ravel()

        # Voltage generating electric fields:
        Phi_env = ((sim.rho_env) / ((sim.ko_env ** 2) * p.eo * p.er)).reshape(cells.X.shape)

        # add in extra boundary conditions for the case of an externally-applied voltage event:
        Phi_env[:, -1] = sim.bound_V['R']
        Phi_env[:, 0] = sim.bound_V['L']
        Phi_env[-1, :] = sim.bound_V['T']
        Phi_env[0, :] = sim.bound_V['B']

        # gradient of the polarization voltage yields the electric field:
        gVex, gVey = fd.gradient(Phi_env, cells.delta)

        # use Hodgkin-Huxley decomposition and recomposition to "smooth" the electric field:
        gVex, gVey = stb.smooth_flux(gVex.reshape(cells.X.shape), gVey.reshape(cells.X.shape), cells)

        # assign to electric field of the system:
        sim.E_env_x = -gVex
        sim.E_env_y = -gVey

        sim.Eme = (sim.E_env_x.ravel()[cells.map_mem2ecm] * cells.mem_vects_flat[:, 2] +
               sim.E_env_y.ravel()[cells.map_mem2ecm] * cells.mem_vects_flat[:, 3])

    else:

        vc = sim.vm/2

        vce = np.zeros(len(cells.xypts))
        vce[cells.map_mem2ecm] = vc[cells.mem_to_cells]

        # Local field potential:
        Phi = -vce

        # Smooth it:
        Phi = fd.integrator(Phi.reshape(cells.X.shape), 0.5)

        # Calculate the gradient:
        gVex, gVey = fd.gradient(cells.delta*Phi.reshape(cells.X.shape), cells.delta)

        # Smooth it:
        gVex, gVey = stb.smooth_flux(gVex.reshape(cells.X.shape), gVey.reshape(cells.X.shape), cells)

        # Assign electric field:
        sim.E_env_x = -gVex
        sim.E_env_y = -gVey

        # Environmental current:
        Jxo = sim.sigma * sim.E_env_x
        Jyo = sim.sigma * sim.E_env_y

        _, sim.J_env_x, sim.J_env_y, _, _, _ = stb.HH_Decomp(Jxo, Jyo, cells)

        # map current from extracellular space to membrane normal
        sim.Eme = (sim.E_env_x.ravel()[cells.map_mem2ecm] * cells.mem_vects_flat[:, 2] +
               sim.E_env_y.ravel()[cells.map_mem2ecm] * cells.mem_vects_flat[:, 3])

        # assign environmental voltage:
        sim.v_env = Phi.ravel()*1


# WASTELANDS (Options)--------------------------------------------------------------------------------------------------

        # Method 2: Integrated charge calculation from currents:-------------------------------------------------------

        #divergence of the environmental current from fluxes in the environment:
        # div_Je = fd.divergence(sim.Jtx, sim.Jty, cells.delta, cells.delta)
        #
        # #divergence of trans-membrane fluxes:
        # div_Je_fromcells = stb.div_env(-sim.Jmem, cells, p).reshape(cells.X.shape)
        # div_Je += div_Je_fromcells
        #
        # #calculate mapped current component from transmembrane fluxes (which is always curl-free):
        # #Important for 100% biophysical correctness, but we can skip adding these in for efficiency
        # Phi_cells = np.dot(cells.lapENVinv, -div_Je_fromcells.ravel())
        # Jex_cells, Jey_cells = fd.gradient(Phi_cells.reshape(cells.X.shape), cells.delta)
        #
        # sim.Jtx += Jex_cells
        # sim.Jty += Jey_cells
        #
        # # #the negative divergence of the total environmental current is the change in charge density:
        # sim.rho_env += -div_Je.ravel()*p.dt


    # ---Method # 1----------------------------------------------------------------------------------------
    # Handle the electric field in the environment using the Screened Poisson equation for currents:

    # divJ = fd.divergence(sim.Jtx, sim.Jty, cells.delta, cells.delta)
    #
    # # lambda_screen = 1e-9
    # v_env = sim.v_env.reshape(cells.X.shape)
    # v_env += -(divJ /((sim.ko_env**2)*p.er*p.eo))*p.dt
    #
    # # if p.sharpness < 1.0:
    # v_env = fd.integrator(v_env, 0.5)
    #
    # # add in extra boundary conditions for the case of an externally-applied voltage event:
    # v_env[:, -1] = sim.bound_V['R']
    # v_env[:, 0] = sim.bound_V['L']
    # v_env[-1, :] = sim.bound_V['T']
    # v_env[0, :] = sim.bound_V['B']
    #
    # gVex, gVey = fd.gradient(v_env, cells.delta)
    #
    # # assign to environmental voltage array:
    # sim.v_env = 1*v_env.ravel()
    #
    # # use Hodgkin-Huxley decomposition and re-composition to "smooth" the electric field:
    # gVex, gVey = stb.smooth_flux(gVex.reshape(cells.X.shape), gVey.reshape(cells.X.shape), cells)
    #
    # # assign to electric field of the system:
    # sim.E_env_x = -gVex
    # sim.E_env_y = -gVey
    #
    # sim.Eme = (sim.E_env_x.ravel()[cells.map_mem2ecm] * cells.mem_vects_flat[:, 2] +
    #            sim.E_env_y.ravel()[cells.map_mem2ecm] * cells.mem_vects_flat[:, 3])

    # --Method #2 ---------------------------------------------------------------

    # The solution to the Screened Poisson Equation in the limit of large screening constant Ko, is simply
    # Phi = +f/Ko2. This makes a perfect voltage estimate for the extracellular space.
    # (the Screened Poisson Equation is Lap(Phi) - ko2 Phi = -rho/(eta)
    # Note that the relative permittivity of the double layer is known to be 6 rather than 80 of pure water
    # (see Srinivasan 2006).

    # v_env = ((sim.rho_env) / ((sim.ko_env ** 2) * p.eo * p.er))
    # v_env = v_env.reshape(cells.X.shape)
    #
    # v_env = fd.integrator(v_env, 0.5)
    #
    # # add in extra boundary conditions for the case of an externally-applied voltage event:
    # v_env[:, -1] = sim.bound_V['R']
    # v_env[:, 0] = sim.bound_V['L']
    # v_env[-1, :] = sim.bound_V['T']
    # v_env[0, :] = sim.bound_V['B']
    #
    # # gradient of the polarization voltage yields the electric field:
    # gVex, gVey = fd.gradient(v_env, cells.delta)
    #
    # # assign to environmental voltage array:
    # sim.v_env = 1*v_env.ravel()
    #
    # # use Hodgkin-Huxley decomposition and recomposition to "smooth" the electric field:
    # gVex, gVey = stb.smooth_flux(gVex.reshape(cells.X.shape), gVey.reshape(cells.X.shape), cells)
    #
    # # assign to electric field of the system:
    # sim.E_env_x = -gVex
    # sim.E_env_y = -gVey
    #
    #
    # sim.Eme = (sim.E_env_x.ravel()[cells.map_mem2ecm] * cells.mem_vects_flat[:, 2] +
    #        sim.E_env_y.ravel()[cells.map_mem2ecm] * cells.mem_vects_flat[:, 3])

    # --Method 3: Local field potentials--------------------------------------------------------------------------

    # env_sig = sim.sigma * sim.D_env_weight # environmental conductivity map
    #
    # # env_sig = sim.sigma_env.reshape(cells.X.shape) # environmental conductivity map
    # # divergence of environmental current, scaled by the conductivity map
    # divJo = fd.divergence(sim.Jtx/env_sig, sim.Jty/env_sig, cells.delta, cells.delta)
    #
    # # set boundary conditions for any applied voltages:
    # divJo[:, -1] = -sim.bound_V['R'] / (cells.delta ** 2)
    # divJo[:, 0] = -sim.bound_V['L'] / (cells.delta ** 2)
    # divJo[-1, :] = -sim.bound_V['T'] / (cells.delta ** 2)
    # divJo[0, :] = -sim.bound_V['B'] / (cells.delta ** 2)
    #
    # # divJo = fd.integrator(divJo.reshape(cells.X.shape), 0.5)
    #
    # # environmental local field potential:
    # Phi = np.dot(cells.lapENVinv, -divJo.ravel())
    #
    # # smooth it:
    # Phi = fd.integrator(Phi.reshape(cells.X.shape), 0.5)
    #
    # # take the gradient:
    # gVex, gVey = fd.gradient(Phi.reshape(cells.X.shape), cells.delta)
    #
    # # use Hodgkin-Huxley decomposition and recomposition to "smooth" the electric field:
    # gVex, gVey = stb.smooth_flux(gVex.reshape(cells.X.shape), gVey.reshape(cells.X.shape), cells)
    #
    # # assign to electric field of the system:
    # sim.E_env_x = -gVex
    # sim.E_env_y = -gVey
    #
    # sim.Eme = (sim.E_env_x.ravel()[cells.map_mem2ecm] * cells.mem_vects_flat[:, 2] +
    #        sim.E_env_y.ravel()[cells.map_mem2ecm] * cells.mem_vects_flat[:, 3])
    #
    # # assign environmental voltage:
    # sim.v_env = Phi.ravel() * 1

    # # Method 4------------------------------------------------------------------------------------------------
    # # Charge in the environmental space is surface charge that screens the cell charge.
    #
    # # Voltage seen in the extracellular region:
    # vce = np.zeros(sim.edl)
    # vce[cells.map_mem2ecm] = -sim.vm/2
    #
    # # env_sig = sim.sigma * sim.D_env_weight # environmental conductivity map
    # # # divergence of environmental current, scaled by the conductivity map
    # # divJo = fd.divergence(sim.Jtx/env_sig, sim.Jty/env_sig, cells.delta, cells.delta)
    # #
    # # # set boundary conditions for any applied voltages:
    # # divJo[:, -1] = -sim.bound_V['R'] / (cells.delta ** 2)
    # # divJo[:, 0] = -sim.bound_V['L'] / (cells.delta ** 2)
    # # divJo[-1, :] = -sim.bound_V['T'] / (cells.delta ** 2)
    # # divJo[0, :] = -sim.bound_V['B'] / (cells.delta ** 2)
    # #
    # # # environmental local field potential:
    # # lfp = np.dot(cells.lapENVinv, -divJo.ravel())
    #
    # # v_env = vce + ((sim.rho_env) / ((sim.ko_env ** 2) * p.eo * p.er))
    #
    # v_env = vce
    #
    # v_env = v_env.reshape(cells.X.shape)
    #
    # v_env = fd.integrator(v_env, 0.5)
    #
    # # assign to environmental voltage array:
    # # sim.v_env = 1*v_env.ravel()
    #
    # # gradient of the polarization voltage yields the electric field:
    # gVex, gVey = fd.gradient(v_env.reshape(cells.X.shape), cells.delta)
    #
    # # use Hodgkin-Huxley decomposition and recomposition to "smooth" the electric field:
    # # gVex, gVey = stb.smooth_flux(gVex.reshape(cells.X.shape), gVey.reshape(cells.X.shape), cells)
    # sim.v_env, gVex, gVey, _, _, _ = stb.HH_Decomp(gVex, gVey, cells)
    #
    # # assign to electric field of the system:
    # sim.E_env_x = -gVex
    # sim.E_env_y = -gVey
    #
    # sim.Eme = (sim.E_env_x.ravel()[cells.map_mem2ecm] * cells.mem_vects_flat[:, 2] +
    #        sim.E_env_y.ravel()[cells.map_mem2ecm] * cells.mem_vects_flat[:, 3])

    # Method 6-------------------------------------------------------------------------------------------------
    # # Calculate voltage from divergence of individual current components:
    #
    # # from transmembrane flux:
    # divJ_tm = np.dot(cells.M_sum_mems, sim.Jmem * cells.mem_sa) / cells.cell_vol
    #
    # # from lateral field:
    # Jc = np.sqrt(sim.J_cell_x ** 2 + sim.J_cell_y ** 2)
    #
    # # divergence of lateral current flow in cells:
    # divJ_lat = np.dot(cells.M_sum_mems, Jc[cells.mem_to_cells] * cells.mem_sa) / cells.cell_vol
    #
    # # total divergence from cell current flow:
    # divJ = divJ_tm + divJ_lat
    #
    # Phi_cells = np.dot(cells.lapGJ_P_inv, -divJ / sim.sigma)
    #
    # Phi_envo = np.zeros(sim.edl)
    # Phi_envo[cells.map_mem2ecm] = -Phi_cells[cells.mem_to_cells]
    #
    # # sum of all contributions:
    # v_env = Phi_envo + ((sim.rho_env) / ((sim.ko_env ** 2) * p.eo * p.er))
    #
    # # v_env = Phi_envo
    #
    # v_env = v_env.reshape(cells.X.shape)
    #
    # v_env = fd.integrator(v_env, 0.5)
    #
    # # assign to environmental voltage array:
    # sim.v_env = 1 * v_env.ravel()
    #
    # # gradient of the polarization voltage yields the electric field:
    # gVex, gVey = fd.gradient(v_env.reshape(cells.X.shape), cells.delta)
    #
    # # use Hodgkin-Huxley decomposition and recomposition to "smooth" the electric field:
    # gVex, gVey = stb.smooth_flux(gVex.reshape(cells.X.shape), gVey.reshape(cells.X.shape), cells)
    # # sim.v_env, gVex, gVey, _, _, _ = stb.HH_Decomp(gVex, gVey, cells)
    #
    # # assign to electric field of the system:
    # sim.E_env_x = -gVex
    # sim.E_env_y = -gVey
    #
    # sim.Eme = (sim.E_env_x.ravel()[cells.map_mem2ecm] * cells.mem_vects_flat[:, 2] +
    #            sim.E_env_y.ravel()[cells.map_mem2ecm] * cells.mem_vects_flat[:, 3])


























