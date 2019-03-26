#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level **export pipeline container** (i.e., container transparently
aggregating runners defined by all available export pipelines) facilities.
'''

# ....................{ IMPORTS                           }....................
from betse.exceptions import BetseSimPipeRunnerUnsatisfiedException
from betse.lib.matplotlib import mplfigure
from betse.science.phase.phasecls import SimPhase
from betse.science.pipe.export.pipeexpcsv import SimPipeExportCSVs
from betse.science.pipe.export.pipeexpanim import SimPipeExportAnimCells
from betse.science.pipe.export.plot.pipeexpplotcell import (
    SimPipeExportPlotCell)
from betse.science.pipe.export.plot.pipeexpplotcells import (
    SimPipeExportPlotCells)
from betse.util.io.log import logs
from betse.util.type.types import type_check, IterableTypes

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
    High-level **simulation export pipeline container** (i.e., container
    transparently aggregating each pipeline runner defined by each available
    simulation export pipeline).

    Attributes
    ----------
    _PIPES_EXPORT: IterableTypes
        Iterable of all available simulation export pipelines.
    '''

    # ..................{ INITIALIZERS                      }..................
    @type_check
    def __init__(self) -> None:
        '''
        Initialize this simulation export pipeline container.
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
        Iterable of all available simulation export pipelines.
        '''

        return self._PIPES_EXPORT

    # ..................{ EXPORTERS                         }..................
    @type_check
    def export(self, phase: SimPhase) -> None:
        '''
        Export (e.g., interactively display, non-interactively save) all
        exports enabled for all export pipelines enabled by the passed
        simulation phase.

        Specifically, for each such pipeline:

        * If that pipeline is enabled (i.e., the :meth:`SimPipeABC.is_enabled`
          method for that pipeline returns ``True`` when passed that phase):

          * For each **pipeline runner subconfiguration** (i.e.,
            :class:`SimConfExportABC` instance in the sequence of these
            instances listed by the :meth:`iter_runners_enabled` method for
            that pipeline):

            * Call the method defined by that pipeline implementing this
              runner, passed both this phase and configuration.
            * If that method reports this runner's requirements to be
              unsatisfied (e.g., due to the current simulation configuration
              disabling extracellular spaces), this runner is ignored with a
              non-fatal warning.

        * Else, log an informative message and ignore that pipeline.

        Parameters
        ----------
        phase: SimPhase
            Current simulation phase.
        '''

        # List of 2-tuples "(runner_method, runner_conf)" yielding the method
        # and configuration of each enabled pipeline runner for each export
        # pipeline for this phase, aggregating the 2-tuples yielded by each
        # SimPipeABC.iter_runners_enabled() generator.
        runners_enabled = []

        # For each available export pipeline...
        for pipe_export in self._PIPES_EXPORT:
            # Initialize this pipeline for this phase.
            pipe_export.init(phase)

            # Append all pipeline runners enabled for this pipeline and phase.
            runners_enabled.extend(pipe_export.iter_runners_enabled(phase))

        # Notify the caller of the range of work performed by this subcommand.
        # namely, notify the caller of the total number of times that this
        # method calls the SimCallbacksBC.progressed() callback or a callback
        # calling that callback (e.g., SimCallbacksBC.progressed_next()).
        phase.callbacks.progress_ranged(progress_max=len(runners_enabled))

        # For the method and configuration of each enabled pipeline runner...
        for runner_method, runner_conf in runners_enabled:
            # Metadata associated with this runner.
            runner_metadata = runner_method.metadata

            # Attempt to...
            try:
                # Run this runner with this phase and configuration.
                runner_method(phase, runner_conf)

                #FIXME: Refactor this low-level kludge from the BETSE codebase
                #into a high-level implementation in the BETSEE codebase. See
                #the prominent "FIXME" comment in the "pipeabc" submodule for
                #preliminary work required to begin doing so. For now, this
                #tragically suffices.

                # Notify the caller of the successful completion of this
                # runner. Since the prior call failed to raise an exception,
                # this runner necessarily succeeded.
                phase.callbacks.progressed_next(
                    status='Exported {} "{}".'.format(
                        runner_metadata.noun_singular_lowercase,
                        runner_metadata.kind))
            # If this runner's requirements are unsatisfied (e.g., due to the
            # current simulation configuration disabling fluid flow), notify
            # the caller of this non-fatal condition and continue.
            except BetseSimPipeRunnerUnsatisfiedException as exception:
                phase.callbacks.progressed_next(
                    status='Excluding {} "{}", as {}.'.format(
                        runner_metadata.noun_singular_lowercase,
                        runner_metadata.kind,
                        exception.reason))
            # Else if this runner raises any other exception, permit this
            # exception to propagate up the callstack without intervention.

        # Unconditionally close all currently open matplotlib figures
        # regardless of whether any of the above runners invoked matplotlib.
        #
        # In theory, each of the runner_method() methods called above should
        # explicitly do so already; in practice, the matplotlib API is
        # sufficiently non-deterministic *AND* resource-consumptive to justify
        # such precautions.
        mplfigure.close_figures_all()

        # Log the directory to which all results were exported.
        logs.log_info('Simulation results exported to:')
        logs.log_info('\t%s', phase.export_dirname)
