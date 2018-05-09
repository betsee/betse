#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level facilities for **pipelining** (i.e., iteratively displaying and/or
exporting) all post-simulation exports -- including plots, animations, and
spreadsheets.
'''

# ....................{ IMPORTS                            }....................
from betse.science.phase.phasecls import SimPhase
from betse.science.phase.phaseenum import SimPhaseKind
from betse.science.export.csv.csvpipe import SimPipeExportCSVs
from betse.science.visual.anim.animpipe import AnimCellsPipe
from betse.science.visual.plot.pipe.plotpipecell import PlotCellPipe
from betse.science.visual.plot.pipe.plotpipecells import PlotCellsPipe
from betse.util.io.log import logs
from betse.util.type.types import type_check

# ....................{ PIPELINES                          }....................
@type_check
def pipeline(phase: SimPhase) -> None:
    '''
    Display and/or save all exports (e.g., comma-separated value (CSV) files,
    plots, animations) enabled in a pipelined fashion for the passed simulation
    phase by the current simulation configuration.

    Parameters
    ----------------------------
    phase: SimPhase
        Current simulation phase.
    '''

    #FIXME: This is atrocious, but unfortunately still required by "plotutil"
    #utility functions transitively called below. After refactoring all such
    #plots to leverage the "PlotCellsABC" and "LayerCellsABC" APIs instead,
    #remove this shameful kludgery.

    # String constant corresponding to the current simulation phase. Valid
    # values include:
    # * `init`, for plotting initialization phase results.
    # * `sim`, for plotting simulation phase results.
    phase.p.plot_type = None
    if phase.kind is SimPhaseKind.SEED:
        phase.p.plot_type = 'seed'
    elif phase.kind is SimPhaseKind.INIT:
        phase.p.plot_type = 'init'
    elif phase.kind is SimPhaseKind.SIM:
        phase.p.plot_type = 'sim'

    # Display and/or save all exports enabled by this configuration in
    # increasing order of expected duration (i.e., from fastest to slowest),
    # trivially improving aesthetic responsiveness for end users:
    #
    # * CSV files, typically exported faster than visuals.
    # * Single-cell plots, typically exported faster than cell cluster plots.
    # * Cell cluster plots, typically exported faster than animations.
    # * Animations.
    SimPipeExportCSVs(phase).run()
    PlotCellPipe(phase).run()
    PlotCellsPipe(phase).run()
    AnimCellsPipe(phase).run()

    # Log the directory to which all results were exported.
    logs.log_info('Results exported to: %s', phase.export_dirname)
