#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level facilities for **pipelining** (i.e., iteratively displaying and/or
exporting) post-simulation animations.
'''

#FIXME: This module would be a *GREAT* candidate for testing out Python 3.5-
#based asynchronicity and parallelization. Ideally, we'd be able to segregate
#the generation of each animation to its own Python process. Verdant shimmers!

# ....................{ IMPORTS                            }....................
import numpy as np
from betse.science.config.export.confvisabc import SimConfVisualListable
from betse.science.math.vector import vectormake
from betse.science.math.vector.field import fieldmake
from betse.science.math.vector.vectorcls import VectorCells
from betse.science.simulate.pipe import piperunreq
from betse.science.simulate.pipe.pipeabc import SimPipelinerExportABC
from betse.science.simulate.pipe.piperun import piperunner
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
    def _runners_conf(self) -> IterableTypes:
        return self._phase.p.anim.after_sim_pipeline

    # ..................{ EXPORTERS ~ current                }..................
    @piperunner(categories=('Current Density', 'Intracellular',))
    def export_current_intra(self, conf: SimConfVisualListable) -> None:
        '''
        Animate all intracellular current densities for the cell cluster for all
        time steps.
        '''

        # Animate this animation.
        AnimCurrent(
            phase=self._phase,
            conf=conf,
            is_current_overlay_only_gj=True,
            label='current_gj',
            figure_title='Intracellular Current',
            colorbar_title='Current Density [uA/cm2]',
        )


    @piperunner(
        categories=('Current Density', 'Extracellular',),
        requirements={piperunreq.ECM,},
    )
    def export_current_extra(self, conf: SimConfVisualListable) -> None:
        '''
        Animate all extracellular current densities for the cell cluster
        environment over all time steps.
        '''

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
    @piperunner(
        categories=('Physical Deformation',),
        requirements={piperunreq.DEFORM,},
    )
    def export_deform(self, conf: SimConfVisualListable) -> None:
        '''
        Animate all physical deformations for the cell cluster over all time
        steps.
        '''

        # Animate this animation.
        AnimateDeformation(
            phase=self._phase,
            conf=conf,
            ani_repeat=True,
            save=self._phase.p.anim.is_after_sim_save,
        )

    # ..................{ EXPORTERS ~ electric               }..................
    @piperunner(categories=('Electric Field', 'Intracellular',))
    def export_electric_intra(self, conf: SimConfVisualListable) -> None:
        '''
        Animate all intracellular electric field lines for the cell cluster for
        all time steps.
        '''

        # Vector field cache of the intracellular electric field over all time steps.
        field = fieldmake.make_electric_intra(phase=self._phase)

        # Vector of all intracellular electric field magnitudes over all time steps.
        # spatially situated at cell centres.
        field_magnitudes = VectorCells(
            phase=self._phase,
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


    @piperunner(
        categories=('Electric Field', 'Extracellular',),
        requirements={piperunreq.ECM,},
    )
    def export_electric_extra(self, conf: SimConfVisualListable) -> None:
        '''
        Animate all extracellular electric field lines for the cell cluster
        environment over all time steps.
        '''

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
    @piperunner(
        categories=('Fluid Flow', 'Intracellular',),
        requirements={piperunreq.FLUID,},
    )
    def export_fluid_intra(self, conf: SimConfVisualListable) -> None:
        '''
        Animate all intracellular fluid flow field lines for the cell cluster
        over all time steps.
        '''

        # Animate this animation.
        AnimVelocityIntracellular(
            phase=self._phase,
            conf=conf,
            label='Velocity_gj',
            figure_title='Intracellular Fluid Velocity',
            colorbar_title='Fluid Velocity [nm/s]',
        )


    @piperunner(
        categories=('Fluid Flow', 'Extracellular',),
        requirements={piperunreq.FLUID, piperunreq.ECM,},
    )
    def export_fluid_extra(self, conf: SimConfVisualListable) -> None:
        '''
        Animate all extracellular fluid flow field lines for the cell cluster
        environment over all time steps.
        '''

        # Animate this animation.
        AnimVelocityExtracellular(
            phase=self._phase,
            conf=conf,
            label='Velocity_ecm',
            figure_title='Extracellular Fluid Velocity',
            colorbar_title='Fluid Velocity [um/s]',
        )

    # ..................{ EXPORTERS ~ ion                    }..................
    @piperunner(
        categories=('Ion Concentration', 'Calcium',),
        requirements={piperunreq.ION_CALCIUM, },
    )
    def export_ion_calcium(self, conf: SimConfVisualListable) -> None:
        '''
        Animate all calcium (i.e., Ca2+) ion concentrations over all time steps.
        '''

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


    @piperunner(
        categories=('Ion Concentration', 'Hydrogen',),
        requirements={piperunreq.ION_HYDROGEN, },
    )
    def export_ion_hydrogen(self, conf: SimConfVisualListable) -> None:
        '''
        Animate all hydrogen (i.e., H+) ion concentrations over all time steps.
        scaled to correspond exactly to pH.
        '''

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

    # ..................{ EXPORTERS ~ junction               }..................
    @piperunner(categories=('Gap Junction', 'Connectivity State',))
    def export_junction_state(self, conf: SimConfVisualListable) -> None:
        '''
        Animate all **gap junction connectivity states** (i.e., relative
        permeabilities of the gap junctions connecting all cell membranes) for
        the cell cluster over all time steps.
        '''

        # Animate this animation.
        AnimGapJuncTimeSeries(
            phase=self._phase,
            conf=conf,
            time_series=self._phase.sim.gjopen_time,
            label='Vmem_gj',
            figure_title='Gap Junction State over Vmem',
            colorbar_title='Voltage [mV]',
        )

    # ..................{ EXPORTERS ~ pressure               }..................
    @piperunner(
        categories=('Cellular Pressure', 'Total',),
        requirements={piperunreq.PRESSURE_TOTAL,},
    )
    def export_pressure_total(self, conf: SimConfVisualListable) -> None:
        '''
        Animate all **cellular pressure totals** (i.e., summations of all
        cellular mechanical and osmotic pressures) for the cell cluster over all
        time steps.
        '''

        # Animate this animation.
        AnimFlatCellsTimeSeries(
            phase=self._phase,
            conf=conf,
            time_series=self._phase.sim.P_cells_time,
            label='Pcell',
            figure_title='Pressure in Cells',
            colorbar_title='Pressure [Pa]',
        )


    @piperunner(
        categories=('Cellular Pressure', 'Osmotic',),
        requirements={piperunreq.PRESSURE_OSMOTIC,},
    )
    def export_pressure_osmotic(self, conf: SimConfVisualListable) -> None:
        '''
        Animate the cellular osmotic pressure over all time steps.
        '''

        # Animate this animation.
        AnimFlatCellsTimeSeries(
            phase=self._phase,
            conf=conf,
            time_series=self._phase.sim.osmo_P_delta_time,
            label='Osmotic Pcell',
            figure_title='Osmotic Pressure in Cells',
            colorbar_title='Pressure [Pa]',
        )

    # ..................{ EXPORTERS ~ pump                   }..................
    @piperunner(
        categories=('Ion Pump', 'Density Factor',),
        requirements={piperunreq.ELECTROOSMOSIS,},
    )
    def export_pump_density(self, conf: SimConfVisualListable) -> None:
        '''
        Animate all cell membrane ion pump density factors for the cell cluster
        over all time steps.
        '''

        # Animate this animation.
        AnimMembraneTimeSeries(
            phase=self._phase,
            conf=conf,
            time_series=self._phase.sim.rho_pump_time,
            label='rhoPump',
            figure_title='Pump Density Factor',
            colorbar_title='mol fraction/m2',
        )

    # ..................{ EXPORTERS ~ voltage                }..................
    @piperunner(categories=('Voltage', 'Transmembrane',))
    def export_voltage_membrane(self, conf: SimConfVisualListable) -> None:
        '''
        Animate all transmembrane voltages (Vmem) over all time steps.
        '''

        # Vector of all cell membrane voltages over all time steps.
        vector = vectormake.make_voltages_membrane(phase=self._phase)

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


    @piperunner(
        categories=('Voltage', 'Extracellular',),
        requirements={piperunreq.ECM,},
    )
    def export_voltage_extra(self, conf: SimConfVisualListable) -> None:
        '''
        Animate all extracellular voltages for the cell cluster environment over
        all time steps.
        '''

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
