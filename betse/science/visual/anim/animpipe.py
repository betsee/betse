#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level facilities for **pipelining** (i.e., iteratively displaying and/or
exporting) post-simulation animations.
'''

#FIXME: This module would be a *GREAT* candidate for testing out Python 3.5-
#based asynchronicity and parallelization. Ideally, we'd be able to segregate
#the generation of each animation to its own Python process. Verdant shimmers!

# ....................{ IMPORTS                            }....................
from betse.science.config.export.visual.confvisabc import (
    SimConfVisualCellsListItem)
from betse.science.math.vector.veccls import VectorCellsCache
from betse.science.phase.pipe.pipeabc import SimPipeExportABC
from betse.science.phase.pipe.piperun import piperunner
from betse.science.phase.require import phasereqs
from betse.science.visual.anim.anim import (
    AnimGapJuncTimeSeries,
    # AnimMembraneTimeSeries,
    AnimVelocityIntracellular,
    AnimVelocityExtracellular,
    AnimFlatCellsTimeSeries,
)
from betse.science.visual.anim.animafter import AnimCellsAfterSolvingLayered
from betse.science.visual.layer.vector import lyrvecabc
from betse.science.visual.layer.vector.lyrvecdiscrete import (
    LayerCellsVectorDiscreteMembranesDeformed)
from betse.science.visual.layer.vector.lyrvecsmooth import (
    LayerCellsVectorSmoothGrids, LayerCellsVectorSmoothRegions)
from betse.science.visual.layer.vectorfield.lyrvecfldquiver import (
    LayerCellsFieldQuiverCells,
    LayerCellsFieldQuiverGrids,
    LayerCellsFieldQuiverMembranes,
)
from betse.science.visual.layer.vectorfield.lyrvecfldstream import (
    LayerCellsFieldStream)
from betse.util.type.types import type_check, IterableTypes

# ....................{ SUBCLASSES                         }....................
#FIXME: Rename to "SimPipeExportAnimsCells" for disambiguity.
class AnimCellsPipe(SimPipeExportABC):
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
        return self._phase.p.anim.anims_after_sim

    # ..................{ EXPORTERS ~ current                }..................
    @piperunner(
        categories=('Current Density', 'Intracellular',),
        requirements=phasereqs.ELECTRIC_CURRENT,
    )
    def export_currents_intra(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Animate all intracellular current densities for the cell cluster for all
        time steps.
        '''

        # Intracellular current density field.
        field = self._phase.cache.vector_field.currents_intra

        # Vector cache of all intracellular current density field magnitudes
        # over all time steps, spatially situated at cell centres.
        field_magnitudes = VectorCellsCache(
            phase=self._phase,
            times_cells_centre=field.times_cells_centre.magnitudes)

        # Layer sequence containing...
        layers = (
            # A lower layer animating these magnitudes.
            LayerCellsVectorSmoothRegions(vector=field_magnitudes),

            # A higher layer animating this field.
            LayerCellsFieldStream(field=field),
        )

        # Animate these layers.
        AnimCellsAfterSolvingLayered(
            phase=self._phase,
            conf=conf,
            layers=layers,
            figure_title='Intracellular Current',
            colorbar_title='Current Density [uA/cm2]',

            # Prefer an alternative colormap.
            colormap=self._phase.p.background_cm,
        )


    @piperunner(
        categories=('Current Density', 'Extracellular',),
        requirements=phasereqs.ELECTRIC_CURRENT_EXTRA,
    )
    def export_currents_extra(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Animate all extracellular current densities for the cell cluster
        environment over all sampled time steps.
        '''

        # Extracellular current density field over all sampled time steps.
        field = self._phase.cache.vector_field.currents_extra

        # Vector of all extracellular current density magnitudes over all time
        # steps, spatially situated at environmental grid space centres.
        field_magnitudes = VectorCellsCache(
            phase=self._phase,
            times_grids_centre=field.times_grids_centre.magnitudes)

        # Layer sequence containing...
        layers = (
            # A lower layer animating these magnitudes.
            LayerCellsVectorSmoothGrids(vector=field_magnitudes),

            # A higher layer animating this field.
            LayerCellsFieldStream(field=field),
        )

        # Animate these layers.
        AnimCellsAfterSolvingLayered(
            phase=self._phase,
            conf=conf,
            layers=layers,
            figure_title='Extracellular Current',
            colorbar_title='Current Density [uA/cm2]',

            # Prefer an alternative colormap.
            colormap=self._phase.p.background_cm,
        )

    # ..................{ EXPORTERS ~ deform                 }..................
    @piperunner(
        categories=('Deformation', 'Total',),
        requirements=phasereqs.DEFORM,
    )
    def export_deform_total(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Animate all **total cellular displacements** (i.e., summations of all
        cellular deformations due to galvanotropic and osmotic pressure body
        forces) for the cell cluster over all sampled time steps.
        '''

        # Total cellular displacements.
        field = self._phase.cache.vector_field.deform_total

        # Vector cache of all total cellular displacement magnitudes over all
        # time steps, spatially situated at cell centres.
        field_magnitudes = VectorCellsCache(
            phase=self._phase,
            times_cells_centre=field.times_cells_centre.magnitudes,
        )

        # Layer sequence containing...
        layers = (
            # A lower layer animating all cell displacement magnitudes.
            LayerCellsVectorDiscreteMembranesDeformed(vector=field_magnitudes),

            # A higher layer animating all cell displacement directionalities.
            LayerCellsFieldQuiverCells(field=field),
        )

        # Animate these layers.
        AnimCellsAfterSolvingLayered(
            phase=self._phase,
            conf=conf,
            layers=layers,
            figure_title='Deformation',
            colorbar_title='Displacement [um]',
        )

    # ..................{ EXPORTERS ~ electric               }..................
    @piperunner(
        categories=('Electric Field', 'Intracellular',),
        requirements=phasereqs.ELECTRIC_FIELD,
    )
    def export_electric_intra(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Animate all intracellular electric field lines for the cell cluster over
        all sampled time steps.
        '''

        # Intracellular electric field.
        field = self._phase.cache.vector_field.electric_intra

        # Vector cache of all intracellular electric field magnitudes over all
        # time steps, spatially situated at cell centres.
        field_magnitudes = VectorCellsCache(
            phase=self._phase,
            times_cells_centre=field.times_cells_centre.magnitudes)

        # Layer sequence containing...
        layers = (
            # A lower layer animating these magnitudes.
            LayerCellsVectorSmoothRegions(vector=field_magnitudes),

            # A higher layer animating this field.
            LayerCellsFieldQuiverCells(field=field),
        )

        # Animate these layers.
        AnimCellsAfterSolvingLayered(
            phase=self._phase,
            conf=conf,
            layers=layers,
            figure_title='Intracellular E Field',
            colorbar_title='Electric Field [V/m]',

            # Prefer an alternative colormap.
            colormap=self._phase.p.background_cm,
        )


    @piperunner(
        categories=('Electric Field', 'Extracellular',),
        requirements=phasereqs.ELECTRIC_FIELD_EXTRA,
    )
    def export_electric_extra(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Animate all extracellular electric field lines for the cell cluster
        environment over all sampled time steps.
        '''

        # Extracellular electric field over all sampled time steps.
        field = self._phase.cache.vector_field.electric_extra

        # Vector of all extracellular electric field magnitudes over all time
        # steps, spatially situated at environmental grid space centres.
        field_magnitudes = VectorCellsCache(
            phase=self._phase,
            times_grids_centre=field.times_grids_centre.magnitudes)

        # Layer sequence containing...
        layers = (
            # A lower layer animating these magnitudes.
            LayerCellsVectorSmoothGrids(vector=field_magnitudes),

            # A higher layer animating this field.
            LayerCellsFieldQuiverGrids(field=field),
        )

        # Animate these layers.
        AnimCellsAfterSolvingLayered(
            phase=self._phase,
            conf=conf,
            layers=layers,
            figure_title='Extracellular E Field',
            colorbar_title='Electric Field [V/m]',

            # Prefer an alternative colormap.
            colormap=self._phase.p.background_cm,
        )

    # ..................{ EXPORTERS ~ fluid                  }..................
    @piperunner(
        categories=('Fluid Flow', 'Intracellular',),
        requirements=phasereqs.FLUID,
    )
    def export_fluid_intra(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Animate all intracellular fluid flow field lines for the cell cluster
        over all sampled time steps.
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
        requirements=phasereqs.FLUID_EXTRA,
    )
    def export_fluid_extra(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Animate all extracellular fluid flow field lines for the cell cluster
        environment over all sampled time steps.
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
        requirements=phasereqs.ION_CALCIUM,
    )
    def export_ion_calcium_intra(
        self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Animate all intracellular calcium ion (Ca2+) concentrations over all
        sampled time steps.
        '''

        # Layer sequence containing only a single layer animating these ion
        # concentrations.
        layers = (lyrvecabc.make_layer(
            phase=self._phase,
            vector=self._phase.cache.vector.ion_calcium_intra),)

        # Animate these layers.
        AnimCellsAfterSolvingLayered(
            phase=self._phase,
            conf=conf,
            layers=layers,
            figure_title='Cytosolic Ca2+',
            colorbar_title='Concentration [nmol/L]',
        )

    # ..................{ EXPORTERS ~ junction               }..................
    # This exporter is solver- and feature-agnostic.
    @piperunner(categories=('Gap Junction', 'Relative Permeability',))
    def export_gj_permeability(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Animate all **gap junction connectivity states** (i.e., relative
        permeabilities of the gap junctions connecting all cell membranes) for
        the cell cluster over all sampled time steps.
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

    # ..................{ EXPORTERS ~ microtubules           }..................
    @piperunner(
        categories=('Microtubules', 'Coherence',),
        requirements=phasereqs.MICROTUBULE,
    )
    def export_microtubule(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Animate the coherence of all cellular microtubules for the cell cluster
        over all sampled time steps.
        '''

        # Sequence containing only a layer animating the vector field of all
        # cellular microtubules over all sampled time steps.
        layers = (LayerCellsFieldQuiverMembranes(
            field=self._phase.cache.vector_field.microtubule),)

        # Animate these layers.
        AnimCellsAfterSolvingLayered(
            phase=self._phase,
            conf=conf,
            layers=layers,
            figure_title='Microtubule arrangement in cells',
        )

    # ..................{ EXPORTERS ~ pressure               }..................
    @piperunner(
        categories=('Pressure', 'Total',),
        requirements=phasereqs.PRESSURE_TOTAL,
    )
    def export_pressure_total(self, conf: SimConfVisualCellsListItem) -> None:
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
        categories=('Pressure', 'Osmotic',),
        requirements=phasereqs.PRESSURE_OSMOTIC,
    )
    def export_pressure_osmotic(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Animate the cellular osmotic pressure over all sampled time steps.
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
    # @piperunner(
    #     categories=('Ion Pump', 'Density Factor',),
    #     requirements={phasereqs.ELECTROOSMOSIS,},
    # )
    # def export_pump_density(self, conf: SimConfVisualCellsListItem) -> None:
    #     '''
    #     Animate all cell membrane ion pump density factors for the cell cluster
    #     over all sampled time steps.
    #     '''
    #
    #     # Animate this animation.
    #     AnimMembraneTimeSeries(
    #         phase=self._phase,
    #         conf=conf,
    #         time_series=self._phase.sim.rho_pump_time,
    #         label='rhoPump',
    #         figure_title='Pump Density Factor',
    #         colorbar_title='mol fraction/m2',
    #     )

    # ..................{ EXPORTERS ~ voltage                }..................
    @piperunner(
        categories=('Voltage', 'Extracellular',),
        requirements=phasereqs.VOLTAGE_EXTRA,
    )
    def export_voltage_extra(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Animate all extracellular voltages for the cell cluster environment over
        all sampled time steps.
        '''

        # Layer sequence containing only a single layer animating these
        # voltages.
        layers = (LayerCellsVectorSmoothGrids(
            vector=self._phase.cache.vector.voltage_extra),)

        # Animate these layers.
        AnimCellsAfterSolvingLayered(
            phase=self._phase,
            conf=conf,
            layers=layers,
            figure_title='Environmental Voltage',
            colorbar_title='Voltage [mV]',
        )


    @piperunner(categories=('Voltage', 'Transmembrane',))
    def export_voltage_membrane(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Animate all transmembrane voltages (Vmem) for the cell cluster over all
        time steps.
        '''

        AnimCellsAfterSolvingLayered(
            phase=self._phase,
            conf=conf,
            layers=(self._phase.cache.vector.layer.voltage_membrane,),
            figure_title='Transmembrane Voltage',
            colorbar_title='Voltage [mV]',
        )


    @piperunner(
        categories=('Voltage', 'Polarity',),
        requirements=phasereqs.VOLTAGE_POLARITY,
    )
    def export_voltage_polarity(
        self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Animate all cellular voltage polarities for the cell cluster over all
        time steps.
        '''

        # Layer sequence containing...
        layers = (
            # A lower layer animating all transmembrane voltages.
            self._phase.cache.vector.layer.voltage_membrane,

            # A higher layer animating all transmembrane voltage polarities.
            LayerCellsFieldQuiverCells(
                field=self._phase.cache.vector_field.voltage_polarity),
        )

        # Animate these layers.
        AnimCellsAfterSolvingLayered(
            phase=self._phase,
            conf=conf,
            layers=layers,
            figure_title='Cell Vmem polarity',
            colorbar_title='Voltage [mV]',
        )
