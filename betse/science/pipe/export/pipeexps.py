#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level **export pipeline container** (i.e., container transparently
aggregating runners defined by all available export pipelines) facilities.
'''

# ....................{ IMPORTS                           }....................
from betse.science.phase.phasecls import SimPhase
from betse.util.io.log import logs
from betse.util.type.types import type_check, IterableTypes

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

# ....................{ CONSTANTS                         }....................
_PIPES_EXPORT_TYPE = (
    SimPipeExportCSVs,
    SimPipeExportPlotCell,
    SimPipeExportPlotCells,
    SimPipeExportAnimCells,
)
'''
Sequence of all export pipeline subclasses in increasing order of expected
duration (i.e., from fastest to slowest), trivially improving aesthetic
responsiveness for end users:

* CSV files, typically exported faster than visuals.
* Single-cell plots, typically exported faster than cell cluster plots.
* Cell cluster plots, typically exported faster than animations.
* Animations.
'''

# ....................{ CLASSES                           }....................
class SimPipesExport(object):
    '''
    High-level **export pipeline container** (i.e., container transparently
    aggregating runners defined by all available export pipelines).

    Attributes
    ----------
    _PIPES_EXPORT: IterableTypes
        Iterable of all available export pipelines.
    '''

    # ..................{ INITIALIZERS                      }..................
    @type_check
    def __init__(self) -> None:
        '''
        Initialize this export pipeline container.
        '''

        # Iterable of all available export pipelines. To avoid dynamically (and
        # inefficiently) regenerating this iterable on each access, this
        # iterable is defined as a tuple rather than generator comprehension.
        self._PIPES_EXPORT = tuple(
            # Export pipeline of this subclass.
            pipe_export_type()
            # For each export pipeline subclass...
            for pipe_export_type in _PIPES_EXPORT_TYPE
        )

    # ..................{ PROPERTIES                        }..................
    # Read-only properties, preventing callers from resetting these attributes.

    @property
    def PIPES_EXPORT(self) -> IterableTypes:
        '''
        Iterable of all available export pipelines.
        '''

        return self._PIPES_EXPORT

    # ..................{ EXPORTERS                         }..................
    @type_check
    def export(self, phase: SimPhase) -> None:
        '''
        Export (e.g., display, save) all available exports enabled by the
        passed simulation phase.

        Parameters
        ----------
        phase: SimPhase
            Current simulation phase.
        '''

        #FIXME: Improve this extremely coarse-grained measure of export progress.
        #Rather than merely hard-coding this to a small magic number, this range of
        #progress should instead be dynamically computed as the total number of
        #CSV files, plots, and animations to be exported by the current collection
        #of export pipelines.

        # Notify the caller of the range of work performed by this subcommand.
        # namely, notify the caller of the total number of times that this
        # method calls the SimCallbacksBC.progressed() callback or a callback
        # calling that callback (e.g., SimCallbacksBC.progressed_next()).
        phase.callbacks.progress_ranged(progress_max=len(self._PIPES_EXPORT))

        # For each export pipeline...
        for pipe_export in self._PIPES_EXPORT:
            # Display and/or save all exports enabled by this phase.
            pipe_export.run(phase)

            # Notify of the caller of the completion of these exports.
            phase.callbacks.progressed_next()

        # Log the directory to which all results were exported.
        logs.log_info('Simulation results exported to:')
        logs.log_info('\t%s', phase.export_dirname)
