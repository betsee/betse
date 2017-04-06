#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level facilities for **pipelining** (i.e., iteratively displaying and/or
exporting) post-simulation plots.
'''

# ....................{ IMPORTS                            }....................
from betse.science.simulate.pipe import piperunreq
from betse.science.visual.plot.pipe.plotpipecell import PlotCellsPipelinerCell
from betse.science.visual.plot.pipe.plotpipecells import PlotCellsPipelinerCells
from betse.science.simulate.simphase import SimPhaseABC
from betse.util.type.types import type_check

# ....................{ OBSOLETE                           }....................
#FIXME: Replace *ALL* functionality defined below with pipeliner run() calls.
@type_check
def pipeline(phase: SimPhaseABC) -> None:
    '''
    Display and/or save all currently enabled plots for the passed simulation
    phase.

    Parameters
    -----------
    phase: SimPhaseABC
        Current simulation phase.
    '''

    #FIXME: Replace *ALL* logic below with the following single call:
    #    pipeliner.run()
    #FIXME: When doing so, note that *ALL* uses of hardcoded plot-specific
    #parameter options (e.g., "self._phase.p.I_ani_min_clr") will need to be
    #refactored to use the general-purpose settings for the current plot.
    #FIXME: Likewise, refactor tests to exercise the new dynamic pipeline schema
    #rather than the obsolete hardcoded schema.

    # If post-simulation plots are disabled, noop.
    if not phase.p.plot.is_after_sim:
       return

    # ..................{ EXPORTERS ~ cell                   }..................
    #FIXME: Consider shifting all single-cell plots into a separate plot
    #pipeline associated with a different YAML key for the following reasons:
    #
    #* Single-cell plots require a different YAML configuration from
    #  multi-cell plots. Specifically:
    #  * Single-cell plots require a "cell index" key.
    #  * Multi-cell plots require a "colorbar" key.
    #* Separating this conjoined pipeline into two distinct pipelines permits
    #  each associated "SimConfList" to be assigned a unique
    #  "SimConfVisualListable" subclass.
    #* Attempting to construct a GUI around the current conjoined pipeline will
    #  be awkward at best, due to the different configuration needs of the two
    #  types of plots.
    #* This conjoined pipeline is already much too large. Attempting to document
    #  this conjoined pipeline in the YAML file is cumbersome and error-prone.

    if phase.p.plot_single_cell_graphs:
        # Post-simulation single-cell plot pipeline producing all such plots.
        pipeliner_cell = PlotCellsPipelinerCell(phase=phase)

        # Plot all cell transmembrane voltages.
        pipeliner_cell.export_voltage_membrane()
        pipeliner_cell.export_voltage_membrane_fft()

        # Plot all cell transmembrane current densities.
        pipeliner_cell.export_currents_membrane()

        # Plot all Na-K-ATPase pump rates.
        pipeliner_cell.export_pump_nakatpase()

        # If calcium is enabled, plot all cell calcium concentrations.
        if phase.p.ions_dict['Ca'] == 1:
            pipeliner_cell.export_ion_calcium()

        # If M anions are enabled, plot all cell M anion concentrations.
        if phase.p.ions_dict['M'] == 1:
            pipeliner_cell.export_ion_m_anion()

        # If potassium is enabled, plot all cell potassium concentrations.
        if phase.p.ions_dict['K'] == 1:
            pipeliner_cell.export_ion_potassium()

        # If sodium is enabled, plot all cell sodium concentrations.
        if phase.p.ions_dict['Na'] == 1:
            pipeliner_cell.export_ion_sodium()

        # If deformations are enabled, plot all cell deformations.
        if phase.p.deformation:
            pipeliner_cell.export_deform()

        # If osmotic pressure is enabled, plot all cell osmotic pressures.
        if phase.p.deform_osmo:
            pipeliner_cell.export_pressure_osmotic()

        # If any pressure is enabled, plot all cell pressure totals.
        if piperunreq.PRESSURE_TOTAL.is_satisfied(phase):
            pipeliner_cell.export_pressure_total()

    # ..................{ EXPORTERS ~ cells                  }..................
    # Post-simulation cell cluster plot pipeline producing all such plots.
    pipeliner_cells = PlotCellsPipelinerCells(phase=phase)

    pipeliner_cells.export_junction_state()
    pipeliner_cells.export_microtubule()
    pipeliner_cells.export_pump_nakatpase()

    # If plotting voltages, do so.
    if phase.p.plot_vm2d:
        pipeliner_cells.export_voltage_membrane()
        pipeliner_cells.export_voltage_membrane_average()

        if phase.p.GHK_calc:
            pipeliner_cells.export_voltage_membrane_ghk()

        if phase.p.sim_ECM:
            pipeliner_cells.export_voltage_extra()

    if phase.p.plot_I2d:
        pipeliner_cells.export_currents_intra()

        if phase.p.sim_ECM:
            pipeliner_cells.export_currents_extra()

    if phase.p.plot_Efield:
        pipeliner_cells.export_electric_intra()

        if phase.p.sim_ECM:
            pipeliner_cells.export_electric_extra()

    if phase.p.plot_ca2d and phase.p.ions_dict['Ca'] == 1:
        pipeliner_cells.export_ion_calcium_intra()

        if phase.p.sim_ECM:
            pipeliner_cells.export_ion_calcium_extra()

    if phase.p.plot_Vel and phase.p.fluid_flow:
        pipeliner_cells.export_fluid_intra()

        if phase.p.sim_ECM:
            pipeliner_cells.export_fluid_extra()

    if phase.p.sim_eosmosis:
        pipeliner_cells.export_channel_density()
        pipeliner_cells.export_pump_density()

    if phase.p.plot_pH2d and phase.p.ions_dict['H'] == 1:
        pipeliner_cells.export_ion_hydrogen_intra()

    if phase.p.plot_P and piperunreq.PRESSURE_TOTAL.is_satisfied(phase):
        pipeliner_cells.export_pressure_total()

    if phase.p.plot_Deformation and phase.p.deformation:
        pipeliner_cells.export_deform()
