#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level simulation phase classes.
'''

# ....................{ IMPORTS                           }....................
from betse.exceptions import BetseSimPhaseException
from betse.science.phase import phasecallbacks
from betse.science.phase.phasecallbacks import SimCallbacksABCOrNoneTypes
from betse.science.phase.phaseenum import SimPhaseKind
from betse.util.type.types import type_check, NoneType

# ....................{ CLASSES                           }....................
#FIXME: Generalize to support dynamically changing cell structure as follows:
#
#* Define a new "cells_time" list indexed by time step for the current phase,
#  each item of which is a "Cells" instance. Most such items are simply
#  references to "cells"; all other items if any will be references to a
#  post-cutting event "Cells" instance. Edge cases include:
#  * For the seed phase, "cells_time" should be "None".
#  * For the init phase, "cells_time" should be a list of the expected length
#    (i.e., the number of initialization time steps), each of whose items is
#    an unconditional reference to "cells".
#  * For the sim phase, "cells_time" should be a list of the expected length
#    (i.e., the number of simulation time steps), only the first of whose items
#    is guaranteed to be a reference to "cells".
#FIXME: While effectively mandatory, the above generalization opens up the
#proverbial can of worms. In particular: pickling. We currently only pickle
#"cells", "p", and "sim". That's an obvious problem. Fortunately, the
#answer is simple (albeit arguably non-ideal): just define a new
#"Simulator.cells_time" list in lieu of a "SimPhase.cells_time" list, but
#otherwise defined exactly as above.
class SimPhase(object):
    '''
    High-level simulation phase, encapsulating all lower-level objects required
    to perform a single phase (e.g., seed, initialization, simulation) of a
    given cell cluster, configuration, and simulation.

    This object principally behaves as a simple container whose:

    * Direct parent is the root-level
      :class:`betse.science.simrunner.SimRunner` object owning all objects
      pertaining to the current cell cluster, configuration, and simulation.
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
        Cell cluster for this phase. If this phase is currently being
        simulated, this object refers to the cluster at the time step being
        simulated; else (i.e., if this phase has already been simulated), this
        object refers to the cluster at the last simulation time step for this
        phase. Since the cluster may change with time while simulating (e.g.,
        due to surgical interventions like cutting events), this object applies
        *only* to this time step. Callers should *not* assume this object to
        uniformly apply to any other time steps of either this or other phases;
        callers should index the :attr:`sim.cells_time` list by time step for
        the cluster at the time step immediately following that time step.
    dyna : betse.science.tissue.tishandler.TissueHandler
        Tissue handler for this phase.
    p : betse.science.parameters.Parameters
        Simulation configuration for this phase.
    sim : betse.science.sim.Simulator
        Simulation for this phase.
    cache : betse.science.phase.cache.cacheabc.SimPhaseCaches
        Simulation cache for this phase.

    Attributes (Low-level: Caller)
    ----------
    callbacks : SimCallbacksABC
        Caller-defined object whose methods are periodically called during this
        phase (e.g., to notify this caller of phase progress).

    Attributes (Low-level: Path)
    ----------
    export_dirname : StrOrNoneTypes
        Absolute path of the top-level directory containing all exported
        results (e.g., plots, animations, CSVs) for this simulation phase if
        this phase is either an initialization or simulation *or* ``None``
        otherwise (i.e., if this phase is a seed).
    '''

    # ..................{ INITIALIZORS                      }..................
    @type_check
    def __init__(
        self,

        # Mandatory parameters.
        kind: SimPhaseKind,
        p: 'betse.science.parameters.Parameters',

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
            this phase (e.g., to notify this caller of phase progress).
            Defaults to ``None``, in which case this defaults to a placeholder
            object whose methods all silently reduce to noops.
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

        #FIXME: Eliminate this crude hack by refactoring all references to
        #"sim.dyna" throughout the codebase to "phase.dyna" instead. Naturally,
        #this will require refactoring all methods referencing "sim.dyna" to
        #accept a "phase: SimPhase" parameter. After doing so, remove this line.
        #Praise be to the multifoliate rose!
        #FIXME: O.K.; thanks to intrepid obsessiveness, only one submodule
        #referencing "sim.dyna" exists in the codebase: the
        #"betse.science.chemistry.networks" submodule. (You know what to do.)
        self.sim.dyna = self.dyna

        #FIXME: Eliminate this crude hack by refactoring all references to
        #"p._run_sim" throughout the codebase to test "phase.kind" instead.
        #Although there only exists one remaining reference, refactoring that
        #reference will be a living nightmare. Why? Because the reference
        #resides in a modulator function -- which are currently called in an
        #extremely dynamic, inscrutable, non-greppable manner. Refactoring this
        #single reference to test "phase.kind" instead will require:
        #
        #* Defining a new high-level
        #  betse.science.math.modulate.make_modulator() function resembling:
        #      def make_modulator(modulator_name: str) -> CallableTypes
        #* Refactoring *ALL* modulator functions to accept a single "phase:
        #  SimPhase" parameter rather than the two "cells, p" parameters
        #  currently accepted by modulator functions.
        #* Refactoring *ALL* calls to modulator functions to instead:
        #  * First call the make_modulator() function. This will ensure that
        #    this *NEVER* happens again, by ensuring that we're able to
        #    reliably manage modulator functions from a single point of
        #    failure.
        #  * Next call the callable returned by that function, passing that
        #    callable the current "phase" object. Of course, this will in turn
        #    necessitate refactoring each callable performing such a call to
        #    transitively accept the current "phase" object. This is Hell.
        #* Refactoring the f_sweep() modulator function to test "phase.kind"
        #  rather than "p._run_sim".
        #
        #To grep for existing calls to modulator functions, grep for both the
        #strings "\bmodulator function\b" and "\bgetattr\(". Our YAML structure
        #appears to reliably refer to the desired modulator function for a
        #given operation via "'modulator function'". *sigh*
        self.p._run_sim = kind is SimPhaseKind.SIM

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

    # ..................{ EXCEPTIONS                        }..................
    @type_check
    def die_unless_kind(self, kind: SimPhaseKind) -> None:
        '''
        Raise an exception unless the kind of this simulation phase is exactly
        the passed kind of such phases (e.g., initialization, simulation).
        '''

        if self.kind is not kind:
            raise BetseSimPhaseException(
                'Simulation phase "{}" not "{}".'.format(self.kind, kind))
