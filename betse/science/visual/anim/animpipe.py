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
from betse.science.simulate.simphase import SimPhaseWeak
from betse.science.simulate.simpipeabc import SimPipelayerABC
from betse.science.vector import vectormake
from betse.science.vector.field import fieldmake
from betse.science.vector.vectorcls import VectorCells
# from betse.science.visual import visualutil
from betse.science.visual.anim.anim import (
    AnimCurrent,
    AnimateDeformation,
    AnimGapJuncTimeSeries,
    AnimMembraneTimeSeries,
    AnimFieldExtracellular,
    AnimVelocityIntracellular,
    AnimVelocityExtracellular,
    AnimFlatCellsTimeSeries,
    AnimEnvTimeSeries
)
from betse.science.visual.anim.animafter import AnimCellsAfterSolvingLayered
from betse.science.visual.layer.field.layerfieldquiver import (
    LayerCellsFieldQuiver)
from betse.science.visual.layer.vector import layervectorsurface
from betse.science.visual.layer.vector.layervectorsurface import (
    LayerCellsVectorSurfaceContinuous)
from betse.util.type.call.memoizers import property_cached
from betse.util.type.types import type_check, SequenceTypes

# ....................{ CLASSES                            }....................
class AnimCellsPipelayer(SimPipelayerABC):
    '''
    **Post-simulation animation pipeline** (i.e., class iteratively creating all
    such animations requested by the current simulation configuration).
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self, *args, **kwargs) -> None:

        # Initialize our superclass with all passed parameters.
        super().__init__(
            *args,
            label_singular='animation',
            label_plural='animations',
            **kwargs)

    # ..................{ SUPERCLASS                         }..................
    @property_cached
    def runner_names(self) -> SequenceTypes:
        '''
        Sequence of the names of all post-simulation animations enabled by this
        simulation configuration.
        '''

        return tuple(anim.name for anim in self._phase.p.anim.postsim_pipeline)

    # ..................{ RUNNERS ~ current                  }..................
    def run_current_intra(self) -> None:
        '''
        Animate the intracellular current density for all time steps.
        '''

        # Log this animation attempt.
        self._die_unless_intra()

        # Animate this animation.
        AnimCurrent(
            sim=self._phase.sim, cells=self._phase.cells, p=self._phase.p,
            is_current_overlay_only_gj=True,
            label='current_gj',
            figure_title='Intracellular Current',
            colorbar_title='Current Density [uA/cm2]',
            is_color_autoscaled=self._phase.p.autoscale_I_ani,
            color_min=self._phase.p.I_ani_min_clr,
            color_max=self._phase.p.I_ani_max_clr,
        )


    def run_current_total(self) -> None:
        '''
        Animate the total current density (i.e., both intra- and extracellular)
        for all time steps.
        '''

        # Raise an exception unless extracellular spaces are enabled.
        self._die_unless_extra()

        # Animate this animation.
        AnimCurrent(
            sim=self._phase.sim, cells=self._phase.cells, p=self._phase.p,
            is_current_overlay_only_gj=False,
            label='current_ecm',
            figure_title='Extracellular Current',
            colorbar_title='Current Density [uA/cm2]',
            is_color_autoscaled=self._phase.p.autoscale_I_ani,
            color_min=self._phase.p.I_ani_min_clr,
            color_max=self._phase.p.I_ani_max_clr,
        )

    # ..................{ RUNNERS ~ deformation              }..................
    def run_deform(self) -> None:
        '''
        Animate physical cellular deformations for all time steps.
        '''

        # Raise an exception unless deformation is enabled.
        self._die_unless(
            is_satisfied=self._phase.p.deformation,
            exception_reason='deformation disabled')

        # Animate this animation.
        AnimateDeformation(
            sim=self._phase.sim, cells=self._phase.cells, p=self._phase.p,
            ani_repeat=True,
            save=self._phase.p.anim.is_after_sim_save,
        )

    # ..................{ RUNNERS ~ electric                 }..................
    def run_electric_intra(self) -> None:
        '''
        Animate the intracellular electric field for all time steps.
        '''

        # Log this animation attempt.
        self._die_unless_intra()

        # Vector field cache of the intracellular electric field for all time steps.
        field = fieldmake.make_electric_intra(
            sim=self._phase.sim, cells=self._phase.cells, p=self._phase.p)

        # Vector of all intracellular electric field magnitudes for all time steps,
        # spatially situated at cell centres.
        field_magnitudes = VectorCells(
            cells=self._phase.cells, p=self._phase.p,
            times_cells_centre=field.times_cells_centre.magnitudes)

        # Sequence of layers consisting of...
        layers = (
            # A lower layer animating these magnitudes.
            LayerCellsVectorSurfaceContinuous(vector=field_magnitudes),

            # A higher layer animating this field.
            LayerCellsFieldQuiver(field=field),
        )

        # Animate these layers.
        AnimCellsAfterSolvingLayered(
            sim=self._phase.sim, cells=self._phase.cells, p=self._phase.p,
            layers=layers,
            label='Efield_gj',
            figure_title='Intracellular E Field',
            colorbar_title='Electric Field [V/m]',
            is_color_autoscaled=self._phase.p.autoscale_Efield_ani,
            color_min=self._phase.p.Efield_ani_min_clr,
            color_max=self._phase.p.Efield_ani_max_clr,

            # Prefer an alternative colormap.
            colormap=self._phase.p.background_cm,
        )


    def run_electric_total(self) -> None:
        '''
        Animate the total electric field (i.e., both intra- and extracellular)
        for all time steps.
        '''

        # Raise an exception unless extracellular spaces are enabled.
        self._die_unless_extra()

        # Animate this animation.
        AnimFieldExtracellular(
            sim=self._phase.sim, cells=self._phase.cells, p=self._phase.p,
            x_time_series=self._phase.sim.efield_ecm_x_time,
            y_time_series=self._phase.sim.efield_ecm_y_time,
            label='Efield_ecm',
            figure_title='Extracellular E Field',
            colorbar_title='Electric Field [V/m]',
            is_color_autoscaled=self._phase.p.autoscale_Efield_ani,
            color_min=self._phase.p.Efield_ani_min_clr,
            color_max=self._phase.p.Efield_ani_max_clr,
        )

    # ..................{ RUNNERS ~ fluid                    }..................
    def run_fluid_intra(self) -> None:
        '''
        Animate the intracellular fluid field for all time steps.
        '''

        # Raise an exception unless fluid flow is enabled.
        self._die_unless(
            is_satisfied=self._phase.p.fluid_flow,
            exception_reason='fluid flow disabled')

        # Animate this animation.
        AnimVelocityIntracellular(
            sim=self._phase.sim, cells=self._phase.cells, p=self._phase.p,
            label='Velocity_gj',
            figure_title='Intracellular Fluid Velocity',
            colorbar_title='Fluid Velocity [nm/s]',
            is_color_autoscaled=self._phase.p.autoscale_Velocity_ani,
            color_min=self._phase.p.Velocity_ani_min_clr,
            color_max=self._phase.p.Velocity_ani_max_clr,
        )


    def run_fluid_total(self) -> None:
        '''
        Animate the total fluid field (i.e., both intra- and extracellular)
        for all time steps.
        '''

        # Raise an exception unless fluid flow and extracellular spaces are
        # enabled.
        self._die_unless(
            is_satisfied=(
                self._phase.p.fluid_flow and self._phase.p.sim_ECM),
            exception_reason=(
                'fluid flow and/or extracellular spaces disabled'))

        # Animate this animation.
        AnimVelocityExtracellular(
            sim=self._phase.sim, cells=self._phase.cells, p=self._phase.p,
            label='Velocity_ecm',
            figure_title='Extracellular Fluid Velocity',
            colorbar_title='Fluid Velocity [um/s]',
            is_color_autoscaled=self._phase.p.autoscale_Velocity_ani,
            color_min=self._phase.p.Velocity_ani_min_clr,
            color_max=self._phase.p.Velocity_ani_max_clr,
        )

    # ..................{ RUNNERS ~ gap junction             }..................
    def run_gap_junction(self) -> None:
        '''
        Animate all intracellular voltages overlayed by gap junction connection
        states for all time steps.
        '''

        # Log this animation attempt.
        self._die_unless_intra()

        # Animate this animation.
        AnimGapJuncTimeSeries(
            sim=self._phase.sim, cells=self._phase.cells, p=self._phase.p,
            time_series=self._phase.sim.gjopen_time,
            label='Vmem_gj',
            figure_title='Gap Junction State over Vmem',
            colorbar_title='Voltage [mV]',
            is_color_autoscaled=self._phase.p.autoscale_Vgj_ani,
            color_min=self._phase.p.Vgj_ani_min_clr,
            color_max=self._phase.p.Vgj_ani_max_clr,
        )

    # ..................{ RUNNERS ~ ion                      }..................
    def run_ion_calcium(self) -> None:
        '''
        Animate all calcium (i.e., Ca2+) ion concentrations for all time steps.
        '''

        # Raise an exception unless the calcium ion is enabled.
        self._die_unless_ion('Ca')

        # Array of all upscaled calcium ion concentrations.
        time_series = [
            1e6*arr[self._phase.sim.iCa] for arr in self._phase.sim.cc_time]

        # Animate this animation.
        AnimFlatCellsTimeSeries(
            sim=self._phase.sim, cells=self._phase.cells, p=self._phase.p,
            time_series=time_series,
            label='Ca',
            figure_title='Cytosolic Ca2+',
            colorbar_title='Concentration [nmol/L]',
            is_color_autoscaled=self._phase.p.autoscale_Ca_ani,
            color_min=self._phase.p.Ca_ani_min_clr,
            color_max=self._phase.p.Ca_ani_max_clr,
        )


    def run_ion_hydrogen(self) -> None:
        '''
        Animate all hydrogen (i.e., H+) ion concentrations for all time steps,
        scaled to correspond exactly to pH.
        '''

        # Raise an exception unless the hydrogen ion is enabled.
        self._die_unless_ion('H')

        # Array of all upscaled calcium ion concentrations.
        time_series = [
            -np.log10(1.0e-3 * arr[self._phase.sim.iH])
            for arr in self._phase.sim.cc_time
        ]

        # Animate this animation.
        AnimFlatCellsTimeSeries(
            sim=self._phase.sim, cells=self._phase.cells, p=self._phase.p,
            time_series=time_series,
            label='pH',
            figure_title='Cytosolic pH',
            colorbar_title='pH',
            is_color_autoscaled=self._phase.p.autoscale_Ca_ani,
            color_min=self._phase.p.Ca_ani_min_clr,
            color_max=self._phase.p.Ca_ani_max_clr,
        )

    # ..................{ RUNNERS ~ pressure                 }..................
    def run_pressure_mechanical(self) -> None:
        '''
        Animate the cellular mechanical pressure for all time steps.
        '''

        # Raise an exception unless mechanical pressure is enabled.
        self._die_unless(
            is_satisfied=self._phase.p.scheduled_options['pressure'] != 0,
            exception_reason='mechanical pressure event disabled')

        # Animate this animation.
        AnimFlatCellsTimeSeries(
            sim=self._phase.sim, cells=self._phase.cells, p=self._phase.p,
            time_series=self._phase.sim.P_cells_time,
            label='Pcell',
            figure_title='Pressure in Cells',
            colorbar_title='Pressure [Pa]',
            is_color_autoscaled=self._phase.p.autoscale_Pcell_ani,
            color_min=self._phase.p.Pcell_ani_min_clr,
            color_max=self._phase.p.Pcell_ani_max_clr,
        )


    def run_pressure_osmotic(self) -> None:
        '''
        Animate the cellular osmotic pressure for all time steps.
        '''

        # Raise an exception unless osmotic pressure is enabled.
        self._die_unless(
            is_satisfied=self._phase.p.deform_osmo,
            exception_reason='osmotic pressure disabled')

        # Animate this animation.
        AnimFlatCellsTimeSeries(
            sim=self._phase.sim, cells=self._phase.cells, p=self._phase.p,
            time_series=self._phase.sim.osmo_P_delta_time,
            label='Osmotic Pcell',
            figure_title='Osmotic Pressure in Cells',
            colorbar_title='Pressure [Pa]',
            is_color_autoscaled=self._phase.p.autoscale_Pcell_ani,
            color_min=self._phase.p.Pcell_ani_min_clr,
            color_max=self._phase.p.Pcell_ani_max_clr,
        )

    # ..................{ RUNNERS ~ pump                     }..................
    def run_pump_density(self) -> None:
        '''
        Animate the cellular membrane pump density factor for all time steps.
        '''

        # Raise an exception unless channel electroosmosis is enabled.
        self._die_unless(
            is_satisfied=self._phase.p.sim_eosmosis,
            exception_reason='channel electroosmosis disabled')

        # Animate this animation.
        AnimMembraneTimeSeries(
            sim=self._phase.sim, cells=self._phase.cells, p=self._phase.p,
            time_series=self._phase.sim.rho_pump_time,
            label='rhoPump',
            figure_title='Pump Density Factor',
            colorbar_title='mol fraction/m2',
            is_color_autoscaled=self._phase.p.autoscale_mem_ani,
            color_min=self._phase.p.mem_ani_min_clr,
            color_max=self._phase.p.mem_ani_max_clr,
        )

    # ..................{ RUNNERS ~ voltage                  }..................
    def run_voltage_intra(self) -> None:
        '''
        Animate all intracellular voltages for all time steps.
        '''

        # Log this animation attempt.
        self._die_unless_intra()

        # Vector of all cell membrane voltages for all time steps.
        vector = vectormake.make_voltages_intra(
            sim=self._phase.sim, cells=self._phase.cells, p=self._phase.p)

        # Sequence of layers, consisting of only one layer animating these voltages
        # as a Gouraud-shaded surface.
        layers = (layervectorsurface.make(p=self._phase.p, vector=vector),)

        # Animate these layers.
        AnimCellsAfterSolvingLayered(
            sim=self._phase.sim, cells=self._phase.cells, p=self._phase.p,
            layers=layers,
            label='Vmem',
            figure_title='Transmembrane Voltage',
            colorbar_title='Voltage [mV]',
            is_color_autoscaled=self._phase.p.autoscale_Vmem_ani,
            color_min=self._phase.p.Vmem_ani_min_clr,
            color_max=self._phase.p.Vmem_ani_max_clr,
        )


    def run_voltage_total(self) -> None:
        '''
        Animate all voltages (i.e., both intra- and extracellular) for all time
        steps.
        '''

        # Raise an exception unless extracellular spaces are enabled.
        self._die_unless_extra()

        # List of environment voltages, indexed by time step.
        venv_time_series = [
            venv.reshape(self._phase.cells.X.shape)*1000
            for venv in self._phase.sim.venv_time
        ]

        # Animate this animation.
        AnimEnvTimeSeries(
            sim=self._phase.sim, cells=self._phase.cells, p=self._phase.p,
            time_series=venv_time_series,
            label='Venv',
            figure_title='Environmental Voltage',
            colorbar_title='Voltage [mV]',
            is_color_autoscaled=self._phase.p.autoscale_venv_ani,
            color_min=self._phase.p.venv_ani_min_clr,
            color_max=self._phase.p.venv_ani_max_clr,
        )

# ....................{ OBSOLETE                           }....................
#FIXME: Replace *ALL* functionality defined below with the "AnimCellsPipelayer"
#class defined above.
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

    # Post-simulation animation pipeline producing all such animations.
    pipelayer = AnimCellsPipelayer(
        phase=SimPhaseWeak(sim=sim, cells=cells, p=p))

    #FIXME: Replace *ALL* logic below with the following single call:
    #    pipelayer.run()
    #FIXME: When doing so, note that *ALL* uses of hardcoded animation-specific
    #parameter options (e.g., "self._phase.p.I_ani_min_clr") will need to be
    #refactored to use the general-purpose settings for the current animation.
    #FIXME: Likewise, refactor tests to exercise the new dynamic pipeline schema
    #rather than the obsolete hardcoded schema.

    if p.ani_ca2d is True and p.ions_dict['Ca'] == 1:
        pipelayer.run_ion_calcium()

    if p.ani_pH2d is True and p.ions_dict['H'] == 1:
        pipelayer.run_ion_hydrogen()

    # If animating cell membrane voltage, do so.
    if p.ani_vm2d:
        pipelayer.run_voltage_intra()

    # Animate environment voltage if requested.
    if p.ani_venv and p.sim_ECM:
        pipelayer.run_voltage_total()

    # If animating gap junction states, do so.
    if p.ani_vmgj2d:
        pipelayer.run_gap_junction()

    # If animating current density, do so.
    if p.ani_I:
        # Always animate intracellular current density.
        pipelayer.run_current_intra()

        # Animate extracellular spaces current if desired as well.
        if p.sim_ECM:
            pipelayer.run_current_total()

    if p.ani_Efield:
        # Always animate the gap junction electric field.
        pipelayer.run_electric_intra()

        # Also animate the extracellular spaces electric field if desired.
        if p.sim_ECM:
            pipelayer.run_electric_total()

    if p.ani_Pcell:
        if p.scheduled_options['pressure'] != 0:
            pipelayer.run_pressure_mechanical()

        if p.deform_osmo:
            pipelayer.run_pressure_osmotic()

    # Display and/or save animations specific to the "sim" simulation phase.
    if p.ani_Velocity and p.fluid_flow:
        # Always animate the gap junction fluid velocity.
        pipelayer.run_fluid_intra()

        # Also animate the extracellular spaces fluid velocity if desired.
        if p.sim_ECM:
            pipelayer.run_fluid_total()

    # Animate deformation if desired.
    if p.ani_Deformation and p.deformation:
        pipelayer.run_deform()

    # Animate the cell membrane pump density factor as a function of time.
    if p.ani_mem and p.sim_eosmosis:
        pipelayer.run_pump_density()
