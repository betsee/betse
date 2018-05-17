#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level **simulation phase** (e.g., cell cluster seeding, initialization, and
simulation) functionality.
'''

# ....................{ IMPORTS                            }....................
from betse.exceptions import BetseSimPhaseException
from betse.science.phase import phasecallbacks
from betse.science.phase.phasecallbacks import (
    SimCallbacksABCOrNoneTypes, SimCallbacksNoop)
from betse.science.phase.phaseenum import SimPhaseKind
from betse.util.type.types import type_check, NoneType

# ....................{ CLASSES                            }....................
class SimPhase(object):
    '''
    High-level simulation phase, encapsulating all lower-level objects required
    to perform a single phase (e.g., seed, initialization, simulation) of a
    given cell cluster, configuration, and simulation.

    This object principally behaves as a simple container whose:

    * Direct parent is the root-level :class:`betse.science.simrunner.SimRunner`
      object owning all objects pertaining to the current cell cluster,
      configuration, and simulation.
    * Direct children are:
      * The current cell cluster.
      * The current simulation.
      * The current simulation configuration.

    Attributes
    ----------
    kind : SimPhaseKind
        Type of this simulation phase.

    Attributes (High-level)
    ----------
    cells : betse.science.cells.Cells
        Cell cluster for this phase.
    dyna : betse.science.tissue.tishandler.TissueHandler
        Tissue handler for this phase.
    p : betse.science.parameters.Parameters
        Simulation configuration for this phase.
    sim : betse.science.sim.Simulator
        Simulation for this phase.
    cache : betse.science.phase.cache.cacheabc.SimPhaseCaches
        Simulation cache for this phase.

    Attributes (Caller)
    ----------
    callbacks : SimCallbacksABC
        Caller-defined object whose methods are periodically called during this
        phase (e.g., to notify this caller of phase progress).

    Attributes (Path)
    ----------
    export_dirname : StrOrNoneTypes
        Absolute path of the top-level directory containing all exported results
        (e.g., plots, animations, CSVs) for this simulation phase if this phase
        is either an initialization or simulation *or* ``None`` otherwise (i.e.,
        if this phase is a seed).
    '''

    # ..................{ INITIALIZORS                       }..................
    @type_check
    def __init__(
        self,

        # Mandatory parameters.
        kind: SimPhaseKind,
        p:     'betse.science.parameters.Parameters',

        # Optional parameters.
        cells: ('betse.science.cells.Cells', NoneType) = None,
        sim:   ('betse.science.sim.Simulator', NoneType) = None,
        callbacks: SimCallbacksABCOrNoneTypes = None,
    ) -> None:
        '''
        Initialize this simulation phase.

        Parameters
        ----------
        kind : SimPhaseKind
            Current simulation phase type.
        p : betse.science.parameters.Parameters
            Current simulation configuration.
        cells : (betse.science.cells.Cells, NoneType)
            Current cell cluster. Defaults to ``None``, in which case this
            defaults to the empty cell cluster for this configuration.
        sim : (betse.science.sim.Simulation, NoneType)
            Current simulation. Defaults to ``None``, in which case this
            defaults to an uninitialized simulation for this configuration.
        callbacks : SimCallbacksABCOrNoneTypes
            Caller-defined object whose methods are periodically called during
            this phase (e.g., to notify this caller of phase progress). Defaults
            to ``None``, in which case this defaults to a placeholder object
            whose methods all silently reduce to noops.
        '''

        # Avoid circular import dependencies.
        from betse.science.cells import Cells
        from betse.science.math.cache.cacheabc import SimPhaseCaches
        from betse.science.sim import Simulator
        from betse.science.tissue.tishandler import TissueHandler

        # Default all unpassed parameters to sane defaults.
        if callbacks is None:
            callbacks = phasecallbacks.make_default()
        if cells is None:
            cells = Cells(p=p)
        if sim is None:
            sim = Simulator(p=p)

        # Classify all passed parameters.
        self.callbacks = callbacks
        self.cells = cells
        self.kind = kind
        self.p = p
        self.sim = sim

        #FIXME: To ensure the "dyna" object is fully initialized, refactor the
        #"TissueHandler" class as follows:
        #
        #* The TissueHandler.tissueProfiles() method should be refactored into a
        #  private TissueHandler._map_tissues_to_cells() method having
        #  the following signature:
        #      @type_check
        #      def _map_tissues_to_cells(self, phase: SimPhase) -> None:
        #* The TissueHandler.__init__() method should be refactered to:
        #  * As the very first logic in that method, validate that the passed
        #    "cells" and "sim" objects have been *FULLY* initialized.
        #  * As the very last logic in that method, call the newly refactored
        #    self._map_tissues_to_cells() method.
        #* The TissueHandler._map_tissues_to_cells() method should *NOT* be
        #  called by any other class. In practice, we probably will need to
        #  break encapsulation and call the newly private
        #  TissueHandler._map_tissues_to_cells() method elsewhere. When we do
        #  so, add a FIXME comment suggesting this to be bad.
        #
        #Alternately, if the _map_tissues_to_cells() method truly *DOES* need
        #to be called elsewhere, simply make it public. *sigh*

        # Classify all remaining high-level objects for this phase.
        self.cache = SimPhaseCaches(phase=self)
        self.dyna = TissueHandler(p=p)

        #FIXME: Isolate exports produced by the "seed" phase to their own
        #directory; for simplicity, these exports currently reuse the same
        #directory as that of the "init" phase.

        # Absolute path of the top-level exports directory for this phase.
        if kind is SimPhaseKind.SEED:
            self.export_dirname = p.init_export_dirname
        elif kind is SimPhaseKind.INIT:
            self.export_dirname = p.init_export_dirname
        elif kind is SimPhaseKind.SIM:
            self.export_dirname = p.sim_export_dirname
        else:
            raise BetseSimPhaseException(
                'Simulation phase "{}" unrecognized.'.format(kind.name))

    # ..................{ EXCEPTIONS                         }..................
    @type_check
    def die_unless_kind(self, kind: SimPhaseKind) -> None:
        '''
        Raise an exception unless the kind of this simulation phase is exactly
        the passed kind of such phases (e.g., seed, initialization, simulation).
        '''

        if self.kind is not kind:
            raise BetseSimPhaseException(
                'Simulation phase "{}" not "{}".'.format(self.kind, kind))
