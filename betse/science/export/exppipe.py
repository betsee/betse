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
from betse.science.phase.phaseenum import SimPhaseKind
from betse.science.export.csv.csvpipe import SimPipeExportCSVs
from betse.science.visual.anim.animpipe import AnimCellsPipe
from betse.science.visual.plot.pipe.plotpipecell import PlotCellPipe
from betse.science.visual.plot.pipe.plotpipecells import PlotCellsPipe
from betse.util.io.log import logs
from betse.util.type.types import type_check

# ....................{ PIPELINES                         }....................
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

    #FIXME: Improve this extremely coarse-grained measure of export progress.
    #Rather than merely hard-coding this to a small magic number, this range of
    #progress should instead be dynamically computed as the total number of
    #CSV files, plots, and animations to be exported by the current collection
    #of export pipelines.

    # Cuumulative number of times that each call of this subcommand calls
    # the SimCallbacksBC.progressed() callback or a callback calling that
    # callback (e.g., SimCallbacksBC.progressed_next()).
    #
    # This magic number *must* be manually synchronized with the implementation
    # of both this method and methods transitively called by this subcommand.
    # Failure to do so *will* result in fatal exceptions. Sadly, there exists
    # no reasonable means of either automating or enforcing this constraint.
    PROGRESS_TOTAL = (
        # Number of progress callbacks performed directly in this method.
        4
    )

    # Notify the caller of the range of work performed by this subcommand.
    phase.callbacks.progress_ranged(progress_max=PROGRESS_TOTAL)

    # Display and/or save all exports enabled by this configuration in
    # increasing order of expected duration (i.e., from fastest to slowest),
    # trivially improving aesthetic responsiveness for end users:
    #
    # * CSV files, typically exported faster than visuals.
    # * Single-cell plots, typically exported faster than cell cluster plots.
    # * Cell cluster plots, typically exported faster than animations.
    # * Animations.
    SimPipeExportCSVs(phase).run()
    phase.callbacks.progressed_next()

    PlotCellPipe(phase).run()
    phase.callbacks.progressed_next()

    PlotCellsPipe(phase).run()
    phase.callbacks.progressed_next()

    AnimCellsPipe(phase).run()
    phase.callbacks.progressed_last()

    # Log the directory to which all results were exported.
    logs.log_info('Results exported to: %s', phase.export_dirname)
