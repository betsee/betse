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
from betse.science.phase.phaseenum import SimPhaseKind
from betse.util.type.types import type_check

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
        Current simulation phase type.

    Attributes (Objects)
    ----------
    cells : betse.science.cells.Cells
        Current cell cluster.
    dyna : betse.science.tissue.tishandler.TissueHandler
        Current tissue handler.
    p : betse.science.parameters.Parameters
        Current simulation configuration.
    sim : betse.science.sim.Simulator
        Current simulation.
    cache : betse.science.phase.cache.cacheabc.SimPhaseCaches
        Current simulation cache.

    Attributes (Path)
    ----------
    save_dirname : StrOrNoneTypes
        Absolute path of the top-level directory containing all exported results
        (e.g., plots, animations, CSVs) for this simulation phase if this phase
        is either an initialization or simulation *or* ``None`` otherwise (i.e.,
        if this phase is a seed).
    '''

    # ..................{ INITIALIZORS                       }..................
    @type_check
    def __init__(
        self,
        kind: SimPhaseKind,

        # Avoid circular import dependencies.
        cells: 'betse.science.cells.Cells',
        p:     'betse.science.parameters.Parameters',
        sim:   'betse.science.sim.Simulator',
    ) -> None:
        '''
        Initialize this simulation phase.

        Parameters
        ----------
        kind : SimPhaseKind
            Current simulation phase type.
        cells : betse.science.cells.Cells
            Current cell cluster.
        p : betse.science.parameters.Parameters
            Current simulation configuration.
        sim : betse.science.sim.Simulation
            Current simulation.
        '''

        # Avoid circular import dependencies.
        from betse.science.math.cache.cacheabc import SimPhaseCaches
        from betse.science.tissue.tishandler import TissueHandler

        # Classify all passed parameters.
        self.kind = kind
        self.cells = cells
        self.p = p
        self.sim = sim

        #FIXME: To ensure the "dyna" object is fully initialized, refactor the
        #"TissueHandler" class as follows:
        #
        #* The TissueHandler.tissueProfiles() method should be refactored into a
        #  private TissueHandler._map_tissue_profiles_to_cells() method having
        #  the following signature:
        #      @type_check
        #      def _map_tissue_profiles_to_cells(self, phase: SimPhase) -> None:
        #* The TissueHandler.__init__() method should be refactered to:
        #  * As the very first logic in that method, validate that the passed
        #    "cells" and "sim" objects have been *FULLY* initialized.
        #  * As the very last logic in that method, call the newly refactored
        #    self._map_tissue_profiles_to_cells() method.
        #* The old TissueHandler.tissueProfiles() method should *NOT* be called
        #  any other class. In practice, we probably will need to break
        #  encapsulation and call the newly private
        #  TissueHandler._map_tissue_profiles_to_cells() method elsewhere. When
        #  we do so, add a FIXME comment suggesting this to be bad.
        #
        #Alternately, if the _map_tissue_profiles_to_cells() method truly *DOES*
        #need to be called elsewhere, simply make it public. *sigh*

        # Classify all remaining high-level objects for this phase.
        self.cache = SimPhaseCaches(self)
        self.dyna = TissueHandler(p)

        #FIXME: Rename the "save_dirname" variable to "export_dirname".
        #FIXME: Isolate exports produced by the "seed" phase to their own
        #directory; for simplicity, such exports currently reuse that of the
        #"init" phase.

        # Absolute path of the top-level exports directory for this phase.
        if kind is SimPhaseKind.SEED:
            self.save_dirname = p.init_export_dirname
        elif kind is SimPhaseKind.INIT:
            self.save_dirname = p.init_export_dirname
        elif kind is SimPhaseKind.SIM:
            self.save_dirname = p.sim_export_dirname
        else:
            raise BetseSimPhaseException(
                'Simulation phase {} unrecognized.'.format(kind.name))

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
