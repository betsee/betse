#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level facilities for **pipelining** (i.e., iteratively displaying and/or
exporting) post-simulation animations.
'''

#FIXME: Force all optional "conf: SimConfVisualListable" parameters
#below to be mandatory.
#FIXME: Refactor the "VisualCellsABC" superclass to accept a single
#"SimConfVisualListable" parameter in lieu of the three current
#"is_color_autoscaled", "color_min", and "color_max" parameters.

#FIXME: This module would be a *GREAT* candidate for testing out Python 3.5-
#based asynchronicity and parallelization. Ideally, we'd be able to segregate
#the generation of each animation to its own Python process. Verdant shimmers!

# ....................{ IMPORTS                            }....................
import numpy as np
from betse.science.config.export.confvisabc import SimConfVisualListable
from betse.science.simulate.simpipeabc import (
    SimPipelinerExportABC, exporter_metadata)
from betse.science.vector import vectormake
from betse.science.vector.field import fieldmake
from betse.science.vector.vectorcls import VectorCells
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
from betse.util.type.types import type_check, IterableTypes

# ....................{ SUBCLASSES                         }....................
class AnimCellsPipeliner(SimPipelinerExportABC):
    '''
    **Post-simulation animation pipeline** (i.e., class iteratively creating all
    post-simulation animations requested by the current simulation
    configuration).
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self, *args, **kwargs) -> None:

        # Initialize our superclass with all passed parameters.
        super().__init__(*args, label_singular='animation', **kwargs)

    # ..................{ SUPERCLASS                         }..................
    @property
    def is_enabled(self) -> bool:
        return self._phase.p.anim.is_after_sim

    @property
    def _runners_conf_enabled(self) -> IterableTypes:
        return self._phase.p.anim.after_sim_pipeline

    # ..................{ EXPORTERS ~ current                }..................
    @exporter_metadata(categories=('Current Density', 'Intracellular'))
    def export_current_intra(self, conf: SimConfVisualListable) -> None:
        '''
        Animate the intracellular current density for all time steps.
        '''

        # Log this animation attempt.
        self._die_unless_intra()

        # Animate this animation.
        AnimCurrent(
            phase=self._phase,
            conf=conf,
            is_current_overlay_only_gj=True,
            label='current_gj',
            figure_title='Intracellular Current',
            colorbar_title='Current Density [uA/cm2]',
        )


    @exporter_metadata(categories=('Current Density', 'Total'))
    def export_current_total(self, conf: SimConfVisualListable) -> None:
        '''
        Animate the total current density (i.e., both intra- and extracellular)
        for all time steps.
        '''

        # Raise an exception unless extracellular spaces are enabled.
        self._die_unless_extra()

        # Animate this animation.
        AnimCurrent(
            phase=self._phase,
            conf=conf,
            is_current_overlay_only_gj=False,
            label='current_ecm',
            figure_title='Extracellular Current',
            colorbar_title='Current Density [uA/cm2]',
        )

    # ..................{ EXPORTERS ~ deform                 }..................
    @exporter_metadata(categories=('Cellular Deformation', 'Physical'))
    def export_deform(self, conf: SimConfVisualListable) -> None:
        '''
        Animate physical cellular deformations for all time steps.
        '''

        # Raise an exception unless deformation is enabled.
        self._die_unless(
            is_satisfied=self._phase.p.deformation,
            exception_reason='deformation disabled')

        # Animate this animation.
        AnimateDeformation(
            phase=self._phase,
            conf=conf,
            ani_repeat=True,
            save=self._phase.p.anim.is_after_sim_save,
        )

    # ..................{ EXPORTERS ~ electric               }..................
    @exporter_metadata(categories=('Electric Field', 'Intracellular'))
    def export_electric_intra(self, conf: SimConfVisualListable) -> None:
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
            phase=self._phase,
            conf=conf,
            layers=layers,
            label='Efield_gj',
            figure_title='Intracellular E Field',
            colorbar_title='Electric Field [V/m]',

            # Prefer an alternative colormap.
            colormap=self._phase.p.background_cm,
        )


    @exporter_metadata(categories=('Electric Field', 'Total'))
    def export_electric_total(self, conf: SimConfVisualListable) -> None:
        '''
        Animate the total electric field (i.e., both intra- and extracellular)
        for all time steps.
        '''

        # Raise an exception unless extracellular spaces are enabled.
        self._die_unless_extra()

        # Animate this animation.
        AnimFieldExtracellular(
            phase=self._phase,
            conf=conf,
            x_time_series=self._phase.sim.efield_ecm_x_time,
            y_time_series=self._phase.sim.efield_ecm_y_time,
            label='Efield_ecm',
            figure_title='Extracellular E Field',
            colorbar_title='Electric Field [V/m]',
        )

    # ..................{ EXPORTERS ~ fluid                  }..................
    @exporter_metadata(categories=('Fluid Flow', 'Intracellular'))
    def export_fluid_intra(self, conf: SimConfVisualListable) -> None:
        '''
        Animate the intracellular fluid flow field for all time steps.
        '''

        # Raise an exception unless fluid flow is enabled.
        self._die_unless(
            is_satisfied=self._phase.p.fluid_flow,
            exception_reason='fluid flow disabled')

        # Animate this animation.
        AnimVelocityIntracellular(
            phase=self._phase,
            conf=conf,
            label='Velocity_gj',
            figure_title='Intracellular Fluid Velocity',
            colorbar_title='Fluid Velocity [nm/s]',
        )


    @exporter_metadata(categories=('Fluid Flow', 'Total'))
    def export_fluid_total(self, conf: SimConfVisualListable) -> None:
        '''
        Animate the total fluid flow field (i.e., both intra- and extracellular)
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
            phase=self._phase,
            conf=conf,
            label='Velocity_ecm',
            figure_title='Extracellular Fluid Velocity',
            colorbar_title='Fluid Velocity [um/s]',
        )

    # ..................{ EXPORTERS ~ ion                    }..................
    @exporter_metadata(categories=('Ion Concentration', 'Calcium'))
    def export_ion_calcium(self, conf: SimConfVisualListable) -> None:
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
            phase=self._phase,
            conf=conf,
            time_series=time_series,
            label='Ca',
            figure_title='Cytosolic Ca2+',
            colorbar_title='Concentration [nmol/L]',
        )


    @exporter_metadata(categories=('Ion Concentration', 'Hydrogen'))
    def export_ion_hydrogen(self, conf: SimConfVisualListable) -> None:
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
            phase=self._phase,
            conf=conf,
            time_series=time_series,
            label='pH',
            figure_title='Cytosolic pH',
            colorbar_title='pH',
        )

    # ..................{ EXPORTERS ~ membrane               }..................
    @exporter_metadata(categories=('Cellular Membrane', 'Gap Junctions'))
    def export_membrane_gap_junction(self, conf: SimConfVisualListable) -> None:
        '''
        Animate all gap junction connectivity states for all time steps.
        '''

        # Log this animation attempt.
        self._die_unless_intra()

        # Animate this animation.
        AnimGapJuncTimeSeries(
            phase=self._phase,
            conf=conf,
            time_series=self._phase.sim.gjopen_time,
            label='Vmem_gj',
            figure_title='Gap Junction State over Vmem',
            colorbar_title='Voltage [mV]',
        )


    @exporter_metadata(categories=('Cellular Membrane', 'Pump Density'))
    def export_membrane_pump_density(self, conf: SimConfVisualListable) -> None:
        '''
        Animate all cellular membrane pump density factors for all time steps.
        '''

        # Raise an exception unless channel electroosmosis is enabled.
        self._die_unless(
            is_satisfied=self._phase.p.sim_eosmosis,
            exception_reason='channel electroosmosis disabled')

        # Animate this animation.
        AnimMembraneTimeSeries(
            phase=self._phase,
            conf=conf,
            time_series=self._phase.sim.rho_pump_time,
            label='rhoPump',
            figure_title='Pump Density Factor',
            colorbar_title='mol fraction/m2',
        )

    # ..................{ EXPORTERS ~ pressure               }..................
    @exporter_metadata(categories=('Cellular Pressure', 'Mechanical'))
    def export_pressure_mechanical(self, conf: SimConfVisualListable) -> None:
        '''
        Animate the cellular mechanical pressure for all time steps.
        '''

        # Raise an exception unless mechanical pressure is enabled.
        self._die_unless(
            is_satisfied=self._phase.p.scheduled_options['pressure'] != 0,
            exception_reason='mechanical pressure event disabled')

        # Animate this animation.
        AnimFlatCellsTimeSeries(
            phase=self._phase,
            conf=conf,
            time_series=self._phase.sim.P_cells_time,
            label='Pcell',
            figure_title='Pressure in Cells',
            colorbar_title='Pressure [Pa]',
        )


    @exporter_metadata(categories=('Cellular Pressure', 'Osmotic'))
    def export_pressure_osmotic(self, conf: SimConfVisualListable) -> None:
        '''
        Animate the cellular osmotic pressure for all time steps.
        '''

        # Raise an exception unless osmotic pressure is enabled.
        self._die_unless(
            is_satisfied=self._phase.p.deform_osmo,
            exception_reason='osmotic pressure disabled')

        # Animate this animation.
        AnimFlatCellsTimeSeries(
            phase=self._phase,
            conf=conf,
            time_series=self._phase.sim.osmo_P_delta_time,
            label='Osmotic Pcell',
            figure_title='Osmotic Pressure in Cells',
            colorbar_title='Pressure [Pa]',
        )

    # ..................{ EXPORTERS ~ voltage                }..................
    @exporter_metadata(categories=('Voltage', 'Intracellular'))
    def export_voltage_intra(self, conf: SimConfVisualListable) -> None:
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
            phase=self._phase,
            conf=conf,
            layers=layers,
            label='Vmem',
            figure_title='Transmembrane Voltage',
            colorbar_title='Voltage [mV]',
        )


    @exporter_metadata(categories=('Voltage', 'Total'))
    def export_voltage_total(self, conf: SimConfVisualListable) -> None:
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
            phase=self._phase,
            conf=conf,
            time_series=venv_time_series,
            label='Venv',
            figure_title='Environmental Voltage',
            colorbar_title='Voltage [mV]',
        )
