#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level facilities for **pipelining** (i.e., iteratively displaying and/or
exporting) all post-simulation exports -- including plots, animations, and
spreadsheets.
'''

# ....................{ IMPORTS                            }....................
from betse.science.export.csv import csvpipe
from betse.science.simulate.simphase import SimPhaseABC, SimPhaseKind
from betse.science.visual.anim import animpipe
from betse.science.visual.plot import plotpipe
from betse.util.io.log import logs
from betse.util.type.types import type_check

# ....................{ PIPELINES                          }....................
@type_check
def pipeline(phase: SimPhaseABC) -> None:
    '''
    Display and/or save all currently enabled simulation exports (e.g., plots,
    animations, spreadsheets) for the passed simulation phase.

    Parameters
    ----------------------------
    phase: SimPhaseABC
        Current simulation phase.
    '''

    #FIXME: This is atrocious, but unfortunately still required by "plotutil"
    #utility functions called below. After refactoring all such plots to
    #leverage the "PlotCellsABC" and "LayerCellsABC" APIs instead, remove this
    #shameful kludgery.

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

    # Display and/or save all animations.
    animpipe.pipeline(phase)

    # Display and/or save all plots.
    plotpipe.pipeline(phase)

    # Display and/or save all spreadsheets.
    csvpipe.pipeline(phase)

    # Log the directory to which all results were exported.
    logs.log_info('Results exported to: %s', phase.save_dirname)
