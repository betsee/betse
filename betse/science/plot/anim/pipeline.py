#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level facilities for displaying and/or saving all enabled animations.
'''

#FIXME: This module would be a *GREAT* candidate for testing out Python 3.5-
#based asynchronicity and parallelization. Ideally, we'd be able to segregate
#the generation of each animation to its own Python process. Verdant shimmers!

# ....................{ IMPORTS                            }....................
import numpy as np
from betse.science.plot.anim.anim import (
    AnimCellsTimeSeries,
    AnimCurrent,
    AnimateDeformation,
    AnimEnvTimeSeries,
    AnimGapJuncTimeSeries,
    AnimMembraneTimeSeries,
    AnimMorphogenTimeSeries,
    AnimFieldIntracellular,
    AnimFieldExtracellular,
    AnimVelocityIntracellular,
    AnimVelocityExtracellular,
)
from betse.science.plot import plot
from betse.util.type import types

# ....................{ PIPELINES                          }....................
#FIXME: To improve memory consumption, split this monolithic function up into
#one function per top-level animation if conditional. Doing so will prevent
#large arrays specific to such conditionals from being retained for the duration
#of the animation pipeline. Against, this really warrants refactoring into a
#class (e.g., to avoid repassing "sim", "cells", and "p" everywhere). Summers!
def anim_all(sim: 'Simulator', cells: 'Cells', p: 'Parameters') -> None:
    '''
    Serially (i.e., in series) display and/or save all enabled animations for
    the current simulation phase if animations are enabled _or_ noop otherwise.

    Parameters
    ----------------------------
    sim : Simulator
        Current simulation.
    cells : Cells
        Current cell cluster.
    p : Parameters
        Current simulation configuration.
    plot_type : str
        String constant corresponding to the current simulation phase. Valid
        values include:
        * `init`, for plotting simulation initialization results.
        * `sim`, for plotting simulation run results.
    '''
    assert types.is_simulator(sim), types.assert_not_simulator(sim)
    assert types.is_cells(cells), types.assert_not_parameters(cells)
    assert types.is_parameters(p), types.assert_not_parameters(p)

    # If animations are disabled, noop.
    if not p.createAnimations:
       return

    # If animating IP3 calcium dynamics, do so.
    if p.ani_ip32d is True and p.Ca_dyn is True:
        AnimCellsTimeSeries(
            sim=sim, cells=cells, p=p,
            time_series=np.array(sim.cIP3_time) * 1000,
            type='IP3',
            figure_title='IP3 concentration',
            colorbar_title='Concentration [umol/L]',
            is_color_autoscaled=p.autoscale_IP3_ani,
            color_min=p.IP3_ani_min_clr,
            color_max=p.IP3_ani_max_clr,
        )

    # If animating voltage-sensitive dye concentration, do so.
    if p.ani_dye2d is True and p.voltage_dye is True:
        # Cell voltage-sensitive dye concentration as a function of time.
        cell_morphogen_time_series = np.array(sim.cDye_time) * 1000

        AnimCellsTimeSeries(
            sim=sim, cells=cells, p=p,
            time_series=cell_morphogen_time_series,
            type='Morph_cell',
            figure_title='Cellular Morphogen Concentration',
            colorbar_title='Concentration [umol/L]',
            is_color_autoscaled=p.autoscale_Dye_ani,
            color_min=p.Dye_ani_min_clr,
            color_max=p.Dye_ani_max_clr,
        )

        if p.sim_ECM is True:
            AnimMorphogenTimeSeries(
                sim=sim, cells=cells, p=p,
                cell_time_series=cell_morphogen_time_series,
                env_time_series=np.array(sim.cDye_env_time)*1000,
                type='Morph_all',
                figure_title='Total Morphogen Concentration',
                colorbar_title='Concentration [umol/L]',
                is_color_autoscaled=p.autoscale_Dye_ani,
                color_min=p.Dye_ani_min_clr,
                color_max=p.Dye_ani_max_clr,
            )


    if p.ani_ca2d is True and p.ions_dict['Ca'] == 1:
        AnimCellsTimeSeries(
            sim=sim, cells=cells, p=p,
            time_series=[1e6*arr[sim.iCa] for arr in sim.cc_time],
            type='Ca',
            figure_title='Cytosolic Ca2+',
            colorbar_title='Concentration [nmol/L]',
            is_color_autoscaled=p.autoscale_Ca_ani,
            color_min=p.Ca_ani_min_clr,
            color_max=p.Ca_ani_max_clr,
        )

    if p.ani_pH2d is True and p.ions_dict['H'] == 1:
        AnimCellsTimeSeries(
            sim=sim, cells=cells, p=p,
            time_series=[-np.log10(1.0e-3*arr[sim.iH]) for arr in sim.cc_time],
            type='pH',
            figure_title='Cytosolic pH',
            colorbar_title='pH',
            is_color_autoscaled=p.autoscale_Ca_ani,
            color_min=p.Ca_ani_min_clr,
            color_max=p.Ca_ani_max_clr,
        )

    if p.ani_vm2d is True:
        if p.sim_ECM is False:
            vmplt = [1000*arr for arr in sim.vm_time[1:]]
            scale_v = [1000*arr for arr in sim.vm_time[1:]]
        else:
            vmplt = [1000*arr for arr in sim.vm_Matrix[1:]]
            scale_v = [1000*arr for arr in sim.vm_time[1:]]

        AnimCellsTimeSeries(
            sim=sim, cells=cells, p=p,
            time_series=vmplt,
            scaling_series=scale_v,
            is_ecm_ignored=False,
            type='Vmem',
            figure_title='Cell Vmem',
            colorbar_title='Voltage [mV]',
            is_color_autoscaled=p.autoscale_Vmem_ani,
            color_min=p.Vmem_ani_min_clr,
            color_max=p.Vmem_ani_max_clr,
        )

    # Animate the gap junction state over cell membrane voltage if desired.
    if p.ani_vmgj2d is True:
        AnimGapJuncTimeSeries(
            sim=sim, cells=cells, p=p,
            cell_time_series=_get_vmem_time_series(sim, p),
            gapjunc_time_series=sim.gjopen_time,
            type='Vmem_gj',
            figure_title='Gap Junction State over Vmem',
            colorbar_title='Voltage [mV]',
            is_color_autoscaled=p.autoscale_Vgj_ani,
            color_min=p.Vgj_ani_min_clr,
            color_max=p.Vgj_ani_max_clr,
        )

    if p.ani_vcell is True and p.sim_ECM is True:
        AnimCellsTimeSeries(
            sim=sim, cells=cells, p=p,
            time_series=[1000*arr for arr in sim.vcell_time],
            type='vcell',

            figure_title='Cell voltage',
            colorbar_title='Voltage [mV]',
            is_color_autoscaled=p.autoscale_vcell_ani,
            color_min=p.vcell_ani_min_clr,
            color_max=p.vcell_ani_max_clr,
        )

    if p.ani_I is True:
        # Always animate the gap junction current.

        AnimCurrent(
            sim=sim, cells=cells, p=p,
            is_overlaying_current_gj_only=True,
            type='current_gj',
            figure_title='Intercellular Current',
            colorbar_title='Current Density [uA/cm2]',
            is_color_autoscaled=p.autoscale_I_ani,
            color_min=p.I_ani_min_clr,
            color_max=p.I_ani_max_clr,
        )

        # Also animate the extracellular spaces current if desired.
        if p.sim_ECM is True:
            AnimCurrent(
                sim=sim, cells=cells, p=p,
                is_overlaying_current_gj_only=False,
                type='current_ecm',
                figure_title='Total Current',
                colorbar_title='Current Density [uA/cm2]',
                is_color_autoscaled=p.autoscale_I_ani,
                color_min=p.I_ani_min_clr,
                color_max=p.I_ani_max_clr,
            )

    if p.ani_Efield is True:
        # Always animate the gap junction electric field.
        AnimFieldIntracellular(
            sim=sim, cells=cells, p=p,
            x_time_series=sim.efield_gj_x_time,
            y_time_series=sim.efield_gj_y_time,
            type='Efield_gj',
            figure_title='Intercellular E Field',
            colorbar_title='Electric Field [V/m]',
            is_color_autoscaled=p.autoscale_Efield_ani,
            color_min=p.Efield_ani_min_clr,
            color_max=p.Efield_ani_max_clr,
        )

        # Also animate the extracellular spaces electric field if desired.
        if p.sim_ECM is True:
            AnimFieldExtracellular(
                sim=sim, cells=cells, p=p,
                x_time_series=sim.efield_ecm_x_time,
                y_time_series=sim.efield_ecm_y_time,
                type='Efield_ecm',
                figure_title='Extracellular E Field',
                colorbar_title='Electric Field [V/m]',
                is_color_autoscaled=p.autoscale_Efield_ani,
                color_min=p.Efield_ani_min_clr,
                color_max=p.Efield_ani_max_clr,
            )

    if p.ani_Pcell is True and p.deform_osmo is True:
        AnimCellsTimeSeries(
            sim=sim, cells=cells, p=p,
            time_series=sim.P_cells_time,
            type='Pcell',
            figure_title='Hydrostatic Pressure in Cells',
            colorbar_title='Pressure [Pa]',
            is_color_autoscaled=p.autoscale_Pcell_ani,
            color_min=p.Pcell_ani_min_clr,
            color_max=p.Pcell_ani_max_clr,
        )

    if p.ani_Pcell is True and p.deform_osmo is True:
        AnimCellsTimeSeries(
            sim=sim, cells=cells, p=p,
            time_series=sim.osmo_P_delta_time,
            type='OsmoP',
            figure_title='Osmotic Pressure in Cells',
            colorbar_title='Pressure [Pa]',
            is_color_autoscaled=p.autoscale_Pcell_ani,
            color_min=p.Pcell_ani_min_clr,
            color_max=p.Pcell_ani_max_clr,
        )

    if p.ani_force is True:
        if p.deform_electro is True:
            AnimFieldIntracellular(
                sim=sim, cells=cells, p=p,
                x_time_series=[(1/p.um)*arr for arr in sim.F_electro_x_time],
                y_time_series=[(1/p.um)*arr for arr in sim.F_electro_y_time],
                type='ElectrostaticFfield',
                figure_title='Electrostatic Body Force',
                colorbar_title='Force [N/cm3]',
                is_color_autoscaled=p.autoscale_force_ani,
                color_min=p.force_ani_min_clr,
                color_max=p.force_ani_max_clr,
            )

        if p.deform_osmo is True:
            AnimFieldIntracellular(
                sim=sim, cells=cells, p=p,
                x_time_series=[(1/p.um)*arr for arr in sim.F_hydro_x_time],
                y_time_series=[(1/p.um)*arr for arr in sim.F_hydro_y_time],
                type='HydroFfield',
                figure_title='Hydrostatic Body Force',
                colorbar_title='Force [N/cm3]',
                is_color_autoscaled=p.autoscale_force_ani,
                color_min=p.force_ani_min_clr,
                color_max=p.force_ani_max_clr,
            )

    # Animate environment voltage if requested.
    if p.ani_venv is True and p.sim_ECM is True:
        # List of environment voltages, indexed by time step.
        venv_time_series = [
            venv.reshape(cells.X.shape)*1000 for venv in sim.venv_time]
        AnimEnvTimeSeries(
            sim=sim, cells=cells, p=p,
            time_series=venv_time_series,
            type='Venv',
            figure_title='Environmental Voltage',
            colorbar_title='Voltage [V]',
            is_color_autoscaled=p.autoscale_venv_ani,
            color_min=p.venv_min_clr,
            color_max=p.venv_max_clr,
        )

    # Display and/or save animations specific to the "sim" simulation phase.
    anim_sim(sim, cells, p)


def anim_sim(sim: 'Simulator', cells: 'Cells', p: 'Parameters') -> None:
    '''
    Serially (i.e., in series) display and/or save all enabled animations if
    the current simulation phase is `sim` _or_ noop otherwise.

    Parameters
    ----------------------------
    sim : Simulator
        Current simulation.
    cells : Cells
        Current cell cluster.
    p : Parameters
        Current simulation configuration.
    '''
    assert types.is_simulator(sim), types.assert_not_simulator(sim)
    assert types.is_cells(cells), types.assert_not_parameters(cells)
    assert types.is_parameters(p), types.assert_not_parameters(p)

    # If the current simulation phase is *NOT* "sim", noop.
    if not sim.run_sim:
       return

    if (p.ani_Velocity is True and p.fluid_flow is True and
        p.deform_electro is True):
        # Always animate the gap junction fluid velocity.
        AnimVelocityIntracellular(
            sim=sim, cells=cells, p=p,
            type='Velocity_gj',
            figure_title='Intercellular Fluid Velocity',
            colorbar_title='Fluid Velocity [nm/s]',
            is_color_autoscaled=p.autoscale_Velocity_ani,
            color_min=p.Velocity_ani_min_clr,
            color_max=p.Velocity_ani_max_clr,
        )

        # Also animate the extracellular spaces fluid velocity if desired.
        if p.sim_ECM is True:
            AnimVelocityExtracellular(
                sim=sim, cells=cells, p=p,
                type='Velocity_ecm',
                figure_title='Extracellular Fluid Velocity',
                colorbar_title='Fluid Velocity [nm/s]',
                is_color_autoscaled=p.autoscale_Velocity_ani,
                color_min=p.Velocity_ani_min_clr,
                color_max=p.Velocity_ani_max_clr,
            )

    # Animate if desired.
    if p.ani_Deformation is True and p.deformation is True:
        AnimateDeformation(
            sim, cells, p,
            ani_repeat=True,
            save=p.saveAnimations,
        )

        # if p.ani_Deformation_type == 'Displacement':
        #     displacement_time_series = [
        #         np.sqrt(cell_dx_series**2 + cell_dy_series**2) * self.p.um
        #         for cell_dx_series, cell_dy_series in zip(
        #            self.sim.dx_cell_time, self.sim.dy_cell_time)]
        #     AnimDeformTimeSeries(
        #         sim=sim, cells=cells, p=p,
        #         cell_time_series=displacement_time_series,
        #         type='Deform_dxdy',
        #         figure_title='Displacement Field and Deformation',
        #         colorbar_title='Displacement [um]',
        #         is_color_autoscaled=p.autoscale_Deformation_ani,
        #         color_min=p.Deformation_ani_min_clr,
        #         color_max=p.Deformation_ani_max_clr,
        #         colormap=p.background_cm,
        #     )
        # elif p.ani_Deformation_type == 'Vmem':
        #     AnimDeformTimeSeries(
        #         sim=sim, cells=cells, p=p,
        #         cell_time_series=_get_vmem_time_series(sim, p),
        #         type='Deform_Vmem',
        #         figure_title='Cell Vmem and Deformation',
        #         colorbar_title='Voltage [mV]',
        #         is_color_autoscaled=p.autoscale_Deformation_ani,
        #         color_min=p.Deformation_ani_min_clr,
        #         color_max=p.Deformation_ani_max_clr,
        #         colormap=p.default_cm,
        #     )

    # Animate the cell membrane pump density factor as a function of time.
    if p.ani_mem is True and p.sim_eosmosis is True:
        AnimMembraneTimeSeries(
            sim=sim, cells=cells, p=p,
            time_series=sim.rho_pump_time,
            type='rhoPump',
            figure_title='Pump Density Factor',
            colorbar_title='mol fraction/m2',
            is_color_autoscaled=p.autoscale_mem_ani,
            color_min=p.mem_ani_min_clr,
            color_max=p.mem_ani_max_clr,
        )

# ....................{ PRIVATE ~ getters                  }....................
#FIXME: Use everywhere above. Since recomputing this is heavy, we probably want
#to refactor this module's functions into class methods. Fair dandylion hair!
def _get_vmem_time_series(sim: 'Simulator', p: 'Parameters') -> list:
    '''
    Get the membrane voltage time series for the current simulation, upscaled
    for use in animations.
    '''
    assert types.is_simulator(sim), types.assert_not_simulator(sim)
    assert types.is_parameters(p), types.assert_not_parameters(p)

    # Scaled membrane voltage time series.
    if p.sim_ECM is False:
        return plot.upscale_data(sim.vm_time)
    else:
        #FIXME: What's the difference between "sim.vcell_time" and
        #"sim.vm_Matrix"? Both the "p.ani_vm2d" and "AnimCellsWhileSolving"
        #animations leverage the latter for extracellular spaces, whereas most
        #animations leverage the former.
        #FIXME: It would seem that "sim.vm_Matrix" is used where continuous
        #plots (e.g., streamplots) are required and "sim.vcell_time" where
        #discrete plots suffice, suggesting we probably want two variants of
        #this method:
        #
        #* _get_vmem_time_series_continuous(), returning "sim.vm_Matrix" for
        #  ECM and "sim.vm_time" for non-ECM.
        #* _get_vmem_time_series_discontinuous(), returning "sim.vcell_time" for
        #  ECM and "sim.vm_time" for non-ECM.
        return plot.upscale_data(sim.vcell_time)
