#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level facilities for **pipelining** (i.e., iteratively displaying and/or
exporting) post-simulation animations.
'''

#FIXME: This module would be a *GREAT* candidate for testing out Python 3.5-
#based asynchronicity and parallelization. Ideally, we'd be able to segregate
#the generation of each animation to its own Python process. Verdant shimmers!

# ....................{ IMPORTS                           }....................
from betse.science.config.export.visual.confexpvisanim import (
    SimConfExportAnimCells)
from betse.science.math.vector.veccls import VectorCellsCache
from betse.science.phase.phasecls import SimPhase
from betse.science.phase.require import phasereqs
from betse.science.pipe.export.pipeexpabc import SimPipeExportABC
from betse.science.pipe.piperun import piperunner
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
from betse.util.type.descriptor.descs import classproperty_readonly
from betse.util.type.types import type_check, SequenceTypes

# ....................{ SUBCLASSES                        }....................
class SimPipeExportAnimCells(SimPipeExportABC):
    '''
    **Post-simulation animation pipeline** (i.e., class iteratively creating
    all post-simulation animations requested by the current simulation
    configuration).
    '''

    # ..................{ SUPERCLASS ~ properties           }..................
    @classproperty_readonly
    def _NOUN_SINGULAR(cls) -> str:
        return 'animation'

    # ..................{ SUPERCLASS ~ methods              }..................
    @type_check
    def iter_runners_conf(self, phase: SimPhase) -> SequenceTypes:
        return phase.p.anim.anims_after_sim

    @type_check
    def _is_enabled(self, phase: SimPhase) -> bool:
        return phase.p.anim.is_after_sim

    # ..................{ EXPORTERS ~ current               }..................
    @piperunner(
        categories=('Current Density', 'Intracellular',),
        requirements=phasereqs.ELECTRIC_CURRENT,
    )
    def export_currents_intra(
        self, phase: SimPhase, conf: SimConfExportAnimCells) -> None:
        '''
        Animate all intracellular current densities for the cell cluster for
        all time steps.
        '''

        # Intracellular current density field.
        field = phase.cache.vector_field.currents_intra

        # Vector cache of all intracellular current density field magnitudes
        # over all time steps, spatially situated at cell centres.
        field_magnitudes = VectorCellsCache(
            phase=phase,
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
            phase=phase,
            conf=conf,
            layers=layers,
            figure_title='Intracellular Current',
            colorbar_title='Current Density [uA/cm2]',

            # Prefer an alternative colormap.
            colormap=phase.p.background_cm,
        )


    @piperunner(
        categories=('Current Density', 'Extracellular',),
        requirements=phasereqs.ELECTRIC_CURRENT_EXTRA,
    )
    def export_currents_extra(
        self, phase: SimPhase, conf: SimConfExportAnimCells) -> None:
        '''
        Animate all extracellular current densities for the cell cluster
        environment over all sampled time steps.
        '''

        # Extracellular current density field over all sampled time steps.
        field = phase.cache.vector_field.currents_extra

        # Vector of all extracellular current density magnitudes over all time
        # steps, spatially situated at environmental grid space centres.
        field_magnitudes = VectorCellsCache(
            phase=phase,
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
            phase=phase,
            conf=conf,
            layers=layers,
            figure_title='Extracellular Current',
            colorbar_title='Current Density [uA/cm2]',

            # Prefer an alternative colormap.
            colormap=phase.p.background_cm,
        )

    # ..................{ EXPORTERS ~ deform                }..................
    @piperunner(
        categories=('Deformation', 'Total',),
        requirements=phasereqs.DEFORM,
    )
    def export_deform_total(
        self, phase: SimPhase, conf: SimConfExportAnimCells) -> None:
        '''
        Animate all **total cellular displacements** (i.e., summations of all
        cellular deformations due to galvanotropic and osmotic pressure body
        forces) for the cell cluster over all sampled time steps.
        '''

        # Total cellular displacements.
        field = phase.cache.vector_field.deform_total

        # Vector cache of all total cellular displacement magnitudes over all
        # time steps, spatially situated at cell centres.
        field_magnitudes = VectorCellsCache(
            phase=phase,
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
            phase=phase,
            conf=conf,
            layers=layers,
            figure_title='Deformation',
            colorbar_title='Displacement [um]',
        )

    # ..................{ EXPORTERS ~ electric              }..................
    @piperunner(
        categories=('Electric Field', 'Intracellular',),
        requirements=phasereqs.ELECTRIC_FIELD,
    )
    def export_electric_intra(
        self, phase: SimPhase, conf: SimConfExportAnimCells) -> None:
        '''
        Animate all intracellular electric field lines for the cell cluster
        over all sampled time steps.
        '''

        # Intracellular electric field.
        field = phase.cache.vector_field.electric_intra

        # Vector cache of all intracellular electric field magnitudes over all
        # time steps, spatially situated at cell centres.
        field_magnitudes = VectorCellsCache(
            phase=phase,
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
            phase=phase,
            conf=conf,
            layers=layers,
            figure_title='Intracellular E Field',
            colorbar_title='Electric Field [V/m]',

            # Prefer an alternative colormap.
            colormap=phase.p.background_cm,
        )


    @piperunner(
        categories=('Electric Field', 'Extracellular',),
        requirements=phasereqs.ELECTRIC_FIELD_EXTRA,
    )
    def export_electric_extra(
        self, phase: SimPhase, conf: SimConfExportAnimCells) -> None:
        '''
        Animate all extracellular electric field lines for the cell cluster
        environment over all sampled time steps.
        '''

        # Extracellular electric field over all sampled time steps.
        field = phase.cache.vector_field.electric_extra

        # Vector of all extracellular electric field magnitudes over all time
        # steps, spatially situated at environmental grid space centres.
        field_magnitudes = VectorCellsCache(
            phase=phase,
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
            phase=phase,
            conf=conf,
            layers=layers,
            figure_title='Extracellular E Field',
            colorbar_title='Electric Field [V/m]',

            # Prefer an alternative colormap.
            colormap=phase.p.background_cm,
        )

    # ..................{ EXPORTERS ~ fluid                 }..................
    @piperunner(
        categories=('Fluid Flow', 'Intracellular',),
        requirements=phasereqs.FLUID,
    )
    def export_fluid_intra(
        self, phase: SimPhase, conf: SimConfExportAnimCells) -> None:
        '''
        Animate all intracellular fluid flow field lines for the cell cluster
        over all sampled time steps.
        '''

        # Animate this animation.
        AnimVelocityIntracellular(
            phase=phase,
            conf=conf,
            kind='Velocity_gj',
            figure_title='Intracellular Fluid Velocity',
            colorbar_title='Fluid Velocity [nm/s]',
        )


    @piperunner(
        categories=('Fluid Flow', 'Extracellular',),
        requirements=phasereqs.FLUID_EXTRA,
    )
    def export_fluid_extra(
        self, phase: SimPhase, conf: SimConfExportAnimCells) -> None:
        '''
        Animate all extracellular fluid flow field lines for the cell cluster
        environment over all sampled time steps.
        '''

        # Animate this animation.
        AnimVelocityExtracellular(
            phase=phase,
            conf=conf,
            kind='Velocity_ecm',
            figure_title='Extracellular Fluid Velocity',
            colorbar_title='Fluid Velocity [um/s]',
        )

    # ..................{ EXPORTERS ~ ion                   }..................
    @piperunner(
        categories=('Ion Concentration', 'Calcium',),
        requirements=phasereqs.ION_CALCIUM,
    )
    def export_ion_calcium_intra(
        self, phase: SimPhase, conf: SimConfExportAnimCells) -> None:
        '''
        Animate all intracellular calcium ion (Ca2+) concentrations over all
        sampled time steps.
        '''

        # Layer sequence containing only a single layer animating these ion
        # concentrations.
        layers = (lyrvecabc.make_layer(
            phase=phase,
            vector=phase.cache.vector.ion_calcium_intra),)

        # Animate these layers.
        AnimCellsAfterSolvingLayered(
            phase=phase,
            conf=conf,
            layers=layers,
            figure_title='Cytosolic Ca2+',
            colorbar_title='Concentration [nmol/L]',
        )

    # ..................{ EXPORTERS ~ junction              }..................
    # This exporter is solver- and feature-agnostic.
    @piperunner(categories=('Gap Junction', 'Relative Permeability',))
    def export_gj_permeability(
        self, phase: SimPhase, conf: SimConfExportAnimCells) -> None:
        '''
        Animate all **gap junction connectivity states** (i.e., relative
        permeabilities of the gap junctions connecting all cell membranes) for
        the cell cluster over all sampled time steps.
        '''

        # Animate this animation.
        AnimGapJuncTimeSeries(
            phase=phase,
            conf=conf,
            time_series=phase.sim.gjopen_time,
            kind='Vmem_gj',
            figure_title='Gap Junction State over Vmem',
            colorbar_title='Voltage [mV]',
        )

    # ..................{ EXPORTERS ~ microtubules          }..................
    #FIXME: Restore after proper support for microtubules is implemented.
    # @piperunner(
    #     categories=('Microtubules', 'Coherence',),
    #     requirements=phasereqs.MICROTUBULE,
    # )
    # def export_microtubule(
    #     self, phase: SimPhase, conf: SimConfExportAnimCells) -> None:
    #     '''
    #     Animate the coherence of all cellular microtubules for the cell cluster
    #     over all sampled time steps.
    #     '''
    #
    #     # Sequence containing only a layer animating the vector field of all
    #     # cellular microtubules over all sampled time steps.
    #     layers = (LayerCellsFieldQuiverMembranes(
    #         field=phase.cache.vector_field.microtubule),)
    #
    #     # Animate these layers.
    #     AnimCellsAfterSolvingLayered(
    #         phase=phase,
    #         conf=conf,
    #         layers=layers,
    #         figure_title='Microtubule arrangement in cells',
    #     )

    # ..................{ EXPORTERS ~ pressure              }..................
    @piperunner(
        categories=('Pressure', 'Total',),
        requirements=phasereqs.PRESSURE_TOTAL,
    )
    def export_pressure_total(
        self, phase: SimPhase, conf: SimConfExportAnimCells) -> None:
        '''
        Animate all **cellular pressure totals** (i.e., summations of all
        cellular mechanical and osmotic pressures) for the cell cluster over
        all time steps.
        '''

        # Animate this animation.
        AnimFlatCellsTimeSeries(
            phase=phase,
            conf=conf,
            time_series=phase.sim.P_cells_time,
            kind='Pcell',
            figure_title='Pressure in Cells',
            colorbar_title='Pressure [Pa]',
        )


    @piperunner(
        categories=('Pressure', 'Osmotic',),
        requirements=phasereqs.PRESSURE_OSMOTIC,
    )
    def export_pressure_osmotic(
        self, phase: SimPhase, conf: SimConfExportAnimCells) -> None:
        '''
        Animate the cellular osmotic pressure over all sampled time steps.
        '''

        # Animate this animation.
        AnimFlatCellsTimeSeries(
            phase=phase,
            conf=conf,
            time_series=phase.sim.osmo_P_delta_time,
            kind='Osmotic Pcell',
            figure_title='Osmotic Pressure in Cells',
            colorbar_title='Pressure [Pa]',
        )

    # ..................{ EXPORTERS ~ pump                   }..................
    # @piperunner(
    #     categories=('Ion Pump', 'Density Factor',),
    #     requirements={phasereqs.ELECTROOSMOSIS,},
    # )
    # def export_pump_density(self, conf: SimConfExportAnimCells) -> None:
    #     '''
    #     Animate all cell membrane ion pump density factors for the cell cluster
    #     over all sampled time steps.
    #     '''
    #
    #     # Animate this animation.
    #     AnimMembraneTimeSeries(
    #         phase=phase,
    #         conf=conf,
    #         time_series=phase.sim.rho_pump_time,
    #         kind='rhoPump',
    #         figure_title='Pump Density Factor',
    #         colorbar_title='mol fraction/m2',
    #     )

    # ..................{ EXPORTERS ~ voltage               }..................
    @piperunner(
        categories=('Voltage', 'Extracellular',),
        requirements=phasereqs.VOLTAGE_EXTRA,
    )
    def export_voltage_extra(
        self, phase: SimPhase, conf: SimConfExportAnimCells) -> None:
        '''
        Animate all extracellular voltages for the cell cluster environment
        over all sampled time steps.
        '''

        # Layer sequence containing only a single layer animating these
        # voltages.
        layers = (LayerCellsVectorSmoothGrids(
            vector=phase.cache.vector.voltage_extra),)

        # Animate these layers.
        AnimCellsAfterSolvingLayered(
            phase=phase,
            conf=conf,
            layers=layers,
            figure_title='Environmental Voltage',
            colorbar_title='Voltage [mV]',
        )


    @piperunner(categories=('Voltage', 'Transmembrane',))
    def export_voltage_membrane(
        self, phase: SimPhase, conf: SimConfExportAnimCells) -> None:
        '''
        Animate all transmembrane voltages (Vmem) for the cell cluster over all
        time steps.
        '''

        AnimCellsAfterSolvingLayered(
            phase=phase,
            conf=conf,
            layers=(phase.cache.vector.layer.voltage_membrane,),
            figure_title='Transmembrane Voltage',
            colorbar_title='Voltage [mV]',
        )


    @piperunner(
        categories=('Voltage', 'Polarity',),
        requirements=phasereqs.VOLTAGE_POLARITY,
    )
    def export_voltage_polarity(
        self, phase: SimPhase, conf: SimConfExportAnimCells) -> None:
        '''
        Animate all cellular voltage polarities for the cell cluster over all
        time steps.
        '''

        # Layer sequence containing...
        layers = (
            # A lower layer animating all transmembrane voltages.
            phase.cache.vector.layer.voltage_membrane,

            # A higher layer animating all transmembrane voltage polarities.
            LayerCellsFieldQuiverCells(
                field=phase.cache.vector_field.voltage_polarity),
        )

        # Animate these layers.
        AnimCellsAfterSolvingLayered(
            phase=phase,
            conf=conf,
            layers=layers,
            figure_title='Cell Vmem polarity',
            colorbar_title='Voltage [mV]',
        )
