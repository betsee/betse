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
from betse.exceptions import BetseSimConfigException
from betse.science.simulate.simphase import SimPhaseABC
from betse.science.vector import vectormake
from betse.science.vector.field import fieldmake
from betse.science.vector.vectorcls import VectorCells
from betse.science.visual import visualutil
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
from betse.util.io.log import logs
from betse.util.type.obj import objects
from betse.util.type.types import type_check

# ....................{ CLASSES                            }....................
#FIXME: Actually use below.
class AnimCellsPipelayer(object):
    '''
    Animation factory running the currently requested **post-simulation
    animation pipeline** (i.e., list of all animations to be animated after
    simulation solving, requested by this simulation configuration).

    Attributes (Private)
    ----------
    _phase : SimPhaseABC
        Current simulation phase.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self, phase: SimPhaseABC) -> None:
        '''
        Initialize this animation pipeline.

        Parameters
        ----------
        phase : SimPhaseABC
            Current simulation phase.
        '''

        # Classify all passed parameters.
        self._phase = phase

    # ..................{ RUNNERS                            }..................
    def run(self) -> None:
        '''
        Run the currently requested post-simulation animation pipeline.

        This method iteratively animates all post-simulation animations
        requested by this simulation configuration (in the order in which .

        Raises
        ----------
        BetseSimConfigException
        '''

        # For each post-simulation animation in this pipeline...
        for anim in self._phase.p.anim.postsim_pipeline:
            # Name of the method animating this animation.
            anim_method_name = 'anim_' + anim.kind

            # Method animating this animation *OR* None if this type of
            # animation is unrecognized.
            anim_method = objects.get_method_or_none(
                obj=self, method_name=anim_method_name)

            # If this type of animation is unrecognized, raise an exception.
            if anim_method is None:
                raise BetseSimConfigException(
                    'Animation type "{}" unrecognized.'.format(
                        anim_method_name))

            # Else, this type of animation is recognized. Run this animation.
            anim_method()

    # ..................{ ANIMATORS ~ electric               }..................
    def anim_electric_intra(self) -> None:
        '''
        Animate the intracellular electric field for all time steps.
        '''

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


    def anim_electric_total(self) -> None:
        '''
        Animate the total electric field (i.e., both intra- and extracellular)
        for all time steps.
        '''

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

    # ..................{ ANIMATORS ~ voltage                }..................
    def anim_voltage_intra(self) -> None:
        '''
        Animate all intracellular voltages for all time steps.
        '''

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

    # If animating cell membrane voltage, do so.
    if p.ani_vm2d is True:
        _anim_voltage_intra(sim=sim, cells=cells, p=p)

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
        _anim_electric_intra(sim=sim, cells=cells, p=p)

        # Also animate the extracellular spaces electric field if desired.
        if p.sim_ECM is True:
            _anim_electric_total(sim=sim, cells=cells, p=p)

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

# ....................{ PRIVATE ~ electric field           }....................
def _anim_electric_intra(
    sim: 'betse.science.sim.Simulator',
    cells: 'betse.science.cells.Cells',
    p: 'betse.science.parameters.Parameters',
) -> None:
    '''
    Animate the intracellular electric field for all time steps.
    '''

    # Vector field cache of the intracellular electric field for all time steps.
    field = fieldmake.make_electric_intra(sim=sim, cells=cells, p=p)

    # Vector of all intracellular electric field magnitudes for all time steps,
    # spatially situated at cell centres.
    field_magnitudes = VectorCells(
        cells=cells, p=p,
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
        sim=sim, cells=cells, p=p, layers=layers,
        label='Efield_gj',
        figure_title='Intracellular E Field',
        colorbar_title='Electric Field [V/m]',
        is_color_autoscaled=p.autoscale_Efield_ani,
        color_min=p.Efield_ani_min_clr,
        color_max=p.Efield_ani_max_clr,

        # Prefer an alternative colormap.
        colormap=p.background_cm,
    )


def _anim_electric_total(
    sim: 'betse.science.sim.Simulator',
    cells: 'betse.science.cells.Cells',
    p: 'betse.science.parameters.Parameters',
) -> None:
    '''
    Animate the total electric field (i.e., both intra- and extracellular) for
    all time steps.
    '''

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
# ....................{ PRIVATE ~ voltage                  }....................
def _anim_voltage_intra(
    sim: 'betse.science.sim.Simulator',
    cells: 'betse.science.cells.Cells',
    p: 'betse.science.parameters.Parameters',
) -> None:
    '''
    Animate all intracellular voltages for all time steps.
    '''

    # Vector of all cell membrane voltages for all time steps.
    vector = vectormake.make_voltages_intra(sim=sim, cells=cells, p=p)

    # Sequence of layers, consisting of only one layer animating these voltages
    # as a Gouraud-shaded surface.
    layers = (layervectorsurface.make(p=p, vector=vector),)

    # Animate these layers.
    AnimCellsAfterSolvingLayered(
        sim=sim, cells=cells, p=p, layers=layers,
        label='Vmem',
        figure_title='Transmembrane Voltage',
        colorbar_title='Voltage [mV]',
        is_color_autoscaled=p.autoscale_Vmem_ani,
        color_min=p.Vmem_ani_min_clr,
        color_max=p.Vmem_ani_max_clr,
    )

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
        return visualutil.upscale_cell_data(sim.vm_time)
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
        return visualutil.upscale_cell_data(sim.vcell_time)
