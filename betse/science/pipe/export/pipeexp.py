#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level facilities for **pipelining** (i.e., iteratively displaying and/or
exporting) all post-simulation exports -- including plots, animations, and
spreadsheets.
'''

# ....................{ IMPORTS                           }....................
from betse.science.phase.phasecls import SimPhase
from betse.util.io.log import logs
from betse.util.type.types import type_check

#FIXME: Refactor away all explicit usage of individual pipeline subclasses.
#Instead, we want this submodule to be able to generically query some master
#singleton object or some such for the set of all currently enabled pipeline
#runners -- which is all this submodule generally cares about. See also similar
#commentary in the "betse_test.fixture.simconf.simconfwrapper" submodule.
from betse.science.pipe.export.pipeexpcsv import SimPipeExportCSVs
from betse.science.pipe.export.pipeexpanim import SimPipeExportAnimCells
from betse.science.pipe.export.plot.pipeexpplotcell import (
    SimPipeExportPlotCell)
from betse.science.pipe.export.plot.pipeexpplotcells import (
    SimPipeExportPlotCells)

# ....................{ PIPELINES                         }....................
@type_check
def pipeline(phase: SimPhase) -> None:
    '''
    Display and/or save all exports (e.g., comma-separated value (CSV) files,
    plots, animations) enabled in a pipelined fashion for the passed simulation
    phase by the current simulation configuration.

    Parameters
    ----------
    phase: SimPhase
        Current simulation phase.
    '''

    # Sequence of all export pipeline subclasses in increasing order of
    # expected duration (i.e., from fastest to slowest), trivially improving
    # aesthetic responsiveness for end users:
    #
    # * CSV files, typically exported faster than visuals.
    # * Single-cell plots, typically exported faster than cell cluster plots.
    # * Cell cluster plots, typically exported faster than animations.
    # * Animations.
    EXPORT_PIPES_TYPE = (
        SimPipeExportCSVs,
        SimPipeExportPlotCell,
        SimPipeExportPlotCells,
        SimPipeExportAnimCells,
    )

    #FIXME: Improve this extremely coarse-grained measure of export progress.
    #Rather than merely hard-coding this to a small magic number, this range of
    #progress should instead be dynamically computed as the total number of
    #CSV files, plots, and animations to be exported by the current collection
    #of export pipelines.

    # Notify the caller of the range of work performed by this subcommand.
    # Specifically, notify the caller of the total number of times that this
    # method calls the SimCallbacksBC.progressed() callback or a callback
    # calling that callback (e.g., SimCallbacksBC.progressed_next()).
    phase.callbacks.progress_ranged(progress_max=len(EXPORT_PIPES_TYPE))

    # For each export pipeline subclass...
    for export_pipe_type in EXPORT_PIPES_TYPE:
        # Export pipeline of this subclass.
        export_pipe = export_pipe_type()

        # Display and/or save all exports enabled by this phase.
        export_pipe.run(phase)

        # Notify of the caller of the completion of these exports.
        phase.callbacks.progressed_next()

    # Log the directory to which all results were exported.
    logs.log_info('Simulation results exported to:')
    logs.log_info('\t%s', phase.export_dirname)
