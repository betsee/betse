#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level facilities for displaying and/or saving all enabled animations.
'''

#FIXME: This module would be a *GREAT* candidate for testing out Python 3.5-
#based asynchronicity and parallelization. Ideally, we'd be able to segregate
#the generation of each animation to its own Python process. Verdant shimmers!

# ....................{ IMPORTS                            }....................
import numpy as np
from betse.science.visual.anim.anim import (
    AnimCellsMembranesData,
    AnimCurrent,
    # AnimDeformTimeSeries,
    AnimateDeformation,
    AnimGapJuncTimeSeries,
    AnimMembraneTimeSeries,
    AnimFieldIntracellular,
    AnimFieldExtracellular,
    AnimVelocityIntracellular,
    AnimVelocityExtracellular,
    AnimFlatCellsTimeSeries,
    AnimEnvTimeSeries
)
from betse.science.visual import visuals
from betse.util.io.log import logs
from betse.util.type import types
from betse.util.type.types import type_check

# ....................{ PIPELINES                          }....................
@type_check
def pipeline_anims(
    sim: 'betse.science.sim.Simulator',
    cells: 'betse.science.cells.Cells',
    p: 'betse.science.parameters.Parameters',
) -> None:
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

    # If post-simulation animations are disabled, noop.
    if not p.anim.is_after_sim:
       return

    # Log animation creation.
    logs.log_info('Creating animations...')



    if p.ani_ca2d is True and p.ions_dict['Ca'] == 1:
        AnimFlatCellsTimeSeries(
            sim=sim, cells=cells, p=p,
            time_series=[1e6*arr[sim.iCa] for arr in sim.cc_time],
            label='Ca',
            figure_title='Cytosolic Ca2+',
            colorbar_title='Concentration [nmol/L]',
            is_color_autoscaled=p.autoscale_Ca_ani,
            color_min=p.Ca_ani_min_clr,
            color_max=p.Ca_ani_max_clr,
        )

    if p.ani_pH2d is True and p.ions_dict['H'] == 1:
        AnimFlatCellsTimeSeries(
            sim=sim, cells=cells, p=p,
            time_series=[-np.log10(1.0e-3*arr[sim.iH]) for arr in sim.cc_time],
            label='pH',
            figure_title='Cytosolic pH',
            colorbar_title='pH',
            is_color_autoscaled=p.autoscale_Ca_ani,
            color_min=p.Ca_ani_min_clr,
            color_max=p.Ca_ani_max_clr,
        )

    if p.ani_vm2d is True:
        AnimCellsMembranesData(
            sim=sim, cells=cells, p=p,
            times_membranes_midpoint_data=visuals.upscale_cell_data(
                sim.vm_time),
            label='Vmem',
            figure_title='Transmembrane Voltage',
            colorbar_title='Voltage [mV]',
            is_color_autoscaled=p.autoscale_Vmem_ani,
            color_min=p.Vmem_ani_min_clr,
            color_max=p.Vmem_ani_max_clr,
        )

    # Animate the gap junction state over cell membrane voltage if desired.
    if p.ani_vmgj2d is True:
        AnimGapJuncTimeSeries(
            sim=sim, cells=cells, p=p,
            time_series=sim.gjopen_time,
            label='Vmem_gj',
            figure_title='Gap Junction State over Vmem',
            colorbar_title='Voltage [mV]',
            is_color_autoscaled=p.autoscale_Vgj_ani,
            color_min=p.Vgj_ani_min_clr,
            color_max=p.Vgj_ani_max_clr,
        )


    if p.ani_I is True:
        # Always animate the gap junction current.
        AnimCurrent(
            sim=sim, cells=cells, p=p,
            is_current_overlay_only_gj=True,
            label='current_gj',
            figure_title='Intracellular Current',
            colorbar_title='Current Density [uA/cm2]',
            is_color_autoscaled=p.autoscale_I_ani,
            color_min=p.I_ani_min_clr,
            color_max=p.I_ani_max_clr,
        )

        # Animate the extracellular spaces current if desired as well.
        if p.sim_ECM is True:
            AnimCurrent(
                sim=sim, cells=cells, p=p,
                is_current_overlay_only_gj=False,
                label='current_ecm',
                figure_title='Extracellular Current',
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
            label='Efield_gj',
            figure_title='Intracellular E Field',
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
                label='Efield_ecm',
                figure_title='Extracellular E Field',
                colorbar_title='Electric Field [V/m]',
                is_color_autoscaled=p.autoscale_Efield_ani,
                color_min=p.Efield_ani_min_clr,
                color_max=p.Efield_ani_max_clr,
            )

    # if np.mean(sim.P_cells_time) != 0.0:

    if p.ani_Pcell is True and np.mean(sim.P_cells_time) != 0.0:
        AnimFlatCellsTimeSeries(
            sim=sim, cells=cells, p=p,
            time_series=sim.P_cells_time,
            label='Pcell',
            figure_title='Pressure in Cells',
            colorbar_title='Pressure [Pa]',
            is_color_autoscaled=p.autoscale_Pcell_ani,
            color_min=p.Pcell_ani_min_clr,
            color_max=p.Pcell_ani_max_clr,
        )


    if p.ani_Pcell is True and p.deform_osmo is True:
        AnimFlatCellsTimeSeries(
            sim=sim, cells=cells, p=p,
            time_series=sim.osmo_P_delta_time,
            label='Osmotic Pcell',
            figure_title='Osmotic Pressure in Cells',
            colorbar_title='Pressure [Pa]',
            is_color_autoscaled=p.autoscale_Pcell_ani,
            color_min=p.Pcell_ani_min_clr,
            color_max=p.Pcell_ani_max_clr,
        )


    # Animate environment voltage if requested.
    if p.ani_venv is True and p.sim_ECM is True:
        # List of environment voltages, indexed by time step.
        venv_time_series = [
            venv.reshape(cells.X.shape)*1000 for venv in sim.venv_time]
        AnimEnvTimeSeries(
            sim=sim, cells=cells, p=p,
            time_series=venv_time_series,
            label='Venv',
            figure_title='Environmental Voltage',
            colorbar_title='Voltage [mV]',
            is_color_autoscaled=p.autoscale_venv_ani,
            color_min=p.venv_ani_min_clr,
            color_max=p.venv_ani_max_clr,
        )

    # Display and/or save animations specific to the "sim" simulation phase.
    if (p.ani_Velocity is True and p.fluid_flow is True):
        # Always animate the gap junction fluid velocity.
        AnimVelocityIntracellular(
            sim=sim, cells=cells, p=p,
            label='Velocity_gj',
            figure_title='Intracellular Fluid Velocity',
            colorbar_title='Fluid Velocity [nm/s]',
            is_color_autoscaled=p.autoscale_Velocity_ani,
            color_min=p.Velocity_ani_min_clr,
            color_max=p.Velocity_ani_max_clr,
        )

        # Also animate the extracellular spaces fluid velocity if desired.
        if p.sim_ECM is True:
            AnimVelocityExtracellular(
                sim=sim, cells=cells, p=p,
                label='Velocity_ecm',
                figure_title='Extracellular Fluid Velocity',
                colorbar_title='Fluid Velocity [um/s]',
                is_color_autoscaled=p.autoscale_Velocity_ani,
                color_min=p.Velocity_ani_min_clr,
                color_max=p.Velocity_ani_max_clr,
            )

    # Animate deformation if desired.
    if p.ani_Deformation is True and p.deformation is True:
        AnimateDeformation(
            sim, cells, p,
            ani_repeat=True,
            save=p.anim.is_after_sim_save,
        )

    # if p.ani_Deformation is True and p.deformation is True:
    #
    #     if p.ani_Deformation_data == 'Displacement':
    #
    #         displacement_time_series = [
    #             np.sqrt(cell_dx_series**2 + cell_dy_series**2) * p.um
    #             for cell_dx_series, cell_dy_series in zip(
    #                sim.dx_cell_time, sim.dy_cell_time)]
    #
    #         AnimDeformTimeSeries(
    #             sim=sim, cells=cells, p=p,
    #             cell_time_series=displacement_time_series,
    #             label='Deform_dxdy',
    #             figure_title='Displacement and Deformation',
    #             colorbar_title='Displacement [um]',
    #             is_color_autoscaled=p.autoscale_Deformation_ani,
    #             color_min=p.Deformation_ani_min_clr,
    #             color_max=p.Deformation_ani_max_clr,
    #             colormap=p.background_cm,
    #         )
    #     elif p.ani_Deformation_data == 'Vmem':
    #
    #         AnimDeformTimeSeries(
    #             sim=sim, cells=cells, p=p,
    #             cell_time_series=_get_vmem_time_series(sim, p),
    #             label='Deform_Vmem',
    #             figure_title='Vmem and Deformation',
    #             colorbar_title='Voltage [mV]',
    #             is_color_autoscaled=p.autoscale_Deformation_ani,
    #             color_min=p.Deformation_ani_min_clr,
    #             color_max=p.Deformation_ani_max_clr,
    #             colormap=p.default_cm,
    #         )



        # AnimDeformTimeSeries(
        #     sim=sim,
        #     cells=cells,
        #     p=p,
        #     ani_repeat=True,
        #     save=p.anim.is_after_sim_save,
        # )


    # Animate the cell membrane pump density factor as a function of time.
    if p.ani_mem is True and p.sim_eosmosis is True:
        AnimMembraneTimeSeries(
            sim=sim, cells=cells, p=p,
            time_series=sim.rho_pump_time,
            label='rhoPump',
            figure_title='Pump Density Factor',
            colorbar_title='mol fraction/m2',
            is_color_autoscaled=p.autoscale_mem_ani,
            color_min=p.mem_ani_min_clr,
            color_max=p.mem_ani_max_clr,
        )


    #
    # anim_sim(sim, cells, p)

# FIXME: these are no longer sim-specific, content moved above, this can be deleted
# def anim_sim(sim: 'Simulator', cells: 'Cells', p: 'Parameters') -> None:
#     '''
#     Serially (i.e., in series) display and/or save all enabled animations if
#     the current simulation phase is `sim` _or_ noop otherwise.
#
#     Parameters
#     ----------------------------
#     sim : Simulator
#         Current simulation.
#     cells : Cells
#         Current cell cluster.
#     p : Parameters
#         Current simulation configuration.
#     '''
#     assert types.is_simulator(sim), types.assert_not_simulator(sim)
#     assert types.is_cells(cells), types.assert_not_parameters(cells)
#     assert types.is_parameters(p), types.assert_not_parameters(p)
#
#     # If the current simulation phase is *NOT* "sim", noop.
#     if not sim.run_sim:
#        return
#
#     if (p.ani_Velocity is True and p.fluid_flow is True):
#         # Always animate the gap junction fluid velocity.
#         AnimVelocityIntracellular(
#             sim=sim, cells=cells, p=p,
#             label='Velocity_gj',
#             figure_title='Intracellular Fluid Velocity',
#             colorbar_title='Fluid Velocity [nm/s]',
#             is_color_autoscaled=p.autoscale_Velocity_ani,
#             color_min=p.Velocity_ani_min_clr,
#             color_max=p.Velocity_ani_max_clr,
#         )
#
#         # Also animate the extracellular spaces fluid velocity if desired.
#         if p.sim_ECM is True:
#             AnimVelocityExtracellular(
#                 sim=sim, cells=cells, p=p,
#                 label='Velocity_ecm',
#                 figure_title='Extracellular Fluid Velocity',
#                 colorbar_title='Fluid Velocity [nm/s]',
#                 is_color_autoscaled=p.autoscale_Velocity_ani,
#                 color_min=p.Velocity_ani_min_clr,
#                 color_max=p.Velocity_ani_max_clr,
#             )
#
#     # Animate if desired.
#     if p.ani_Deformation is True and p.deformation is True:
#         AnimateDeformation(
#             sim, cells, p,
#             ani_repeat=True,
#             save=p.anim.is_after_sim_save,
#         )
#
#         # if p.ani_Deformation_type == 'Displacement':
#         #     displacement_time_series = [
#         #         np.sqrt(cell_dx_series**2 + cell_dy_series**2) * self.p.um
#         #         for cell_dx_series, cell_dy_series in zip(
#         #            self.sim.dx_cell_time, self.sim.dy_cell_time)]
#         #     AnimDeformTimeSeries(
#         #         sim=sim, cells=cells, p=p,
#         #         cell_time_series=displacement_time_series,
#         #         label='Deform_dxdy',
#         #         figure_title='Displacement Field and Deformation',
#         #         colorbar_title='Displacement [um]',
#         #         is_color_autoscaled=p.autoscale_Deformation_ani,
#         #         color_min=p.Deformation_ani_min_clr,
#         #         color_max=p.Deformation_ani_max_clr,
#         #         colormap=p.background_cm,
#         #     )
#         # elif p.ani_Deformation_type == 'Vmem':
#         #     AnimDeformTimeSeries(
#         #         sim=sim, cells=cells, p=p,
#         #         cell_time_series=_get_vmem_time_series(sim, p),
#         #         label='Deform_Vmem',
#         #         figure_title='Cell Vmem and Deformation',
#         #         colorbar_title='Voltage [mV]',
#         #         is_color_autoscaled=p.autoscale_Deformation_ani,
#         #         color_min=p.Deformation_ani_min_clr,
#         #         color_max=p.Deformation_ani_max_clr,
#         #         colormap=p.default_cm,
#         #     )
#
#     # Animate the cell membrane pump density factor as a function of time.
#     if p.ani_mem is True and p.sim_eosmosis is True:
#         AnimMembraneTimeSeries(
#             sim=sim, cells=cells, p=p,
#             time_series=sim.rho_pump_time,
#             label='rhoPump',
#             figure_title='Pump Density Factor',
#             colorbar_title='mol fraction/m2',
#             is_color_autoscaled=p.autoscale_mem_ani,
#             color_min=p.mem_ani_min_clr,
#             color_max=p.mem_ani_max_clr,
#         )

# ....................{ PRIVATE ~ getters                  }....................
#FIXME: Use everywhere above. Since recomputing this is heavy, we probably want
#to refactor this module's functions into class methods. Fair dandylion hair!
@type_check
def _get_vmem_time_series(
    sim: 'betse.science.sim.Simulator',
    p: 'betse.science.parameters.Parameters',
) -> list:
    '''
    Get the membrane voltage time series for the current simulation, upscaled
    for use in animations.
    '''

    # Scaled membrane voltage time series.
    if p.sim_ECM is False:
        return visuals.upscale_cell_data(sim.vm_time)
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
        return visuals.upscale_cell_data(sim.vcell_time)
