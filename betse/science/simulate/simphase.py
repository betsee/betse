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
from betse.util.type.enums import EnumOrdered
from betse.util.type.types import type_check

# ....................{ ENUMS                              }....................
SimPhaseKind = EnumOrdered('SimPhaseKind', ('SEED', 'INIT', 'SIM',))
'''
Ordered enumeration of all possible simulation phases.

Each member of this enumeration is arbitrarily comparable to each other member.
Each member's value is less than that of another member's value if and only if
the former simulation phase is performed _before_ the latter. Specifically,
this enumeration is a total ordering such that:

    >>> SimPhaseKind.SEED < SimPhaseKind.INIT < SimPhaseKind.SIM
    True

Attributes
----------
SEED : enum
    Seed simulation phase, as implemented by the
    :meth:`betse.science.simrunner.SimRunner.seed` method. This phase creates
    the cell cluster and caches this cluster to an output file.
INIT : enum
    Initialization simulation phase, as implemented by the
    :meth:`betse.science.simrunner.SimRunner.init` method. This phase
    initializes the previously created cell cluster from a cached input file
    and caches this initialization to an output file.
SIM : enum
    Proper simulation phase, as implemented by the
    :meth:`betse.science.simrunner.SimRunner.sim` method. This phase simulates
    the previously initialized cell cluster from a cached input file and caches
    this simulation to an output file.
'''

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
    p : betse.science.parameters.Parameters
        Current simulation configuration.
    sim : betse.science.sim.Simulator
        Current simulation.
    cache : betse.science.simulate.cache.cacheabc.SimPhaseCaches
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
        sim:   'betse.science.sim.Simulator',
        cells: 'betse.science.cells.Cells',
        p:     'betse.science.parameters.Parameters',
    ) -> None:
        '''
        Initialize this simulation phase instance.

        Parameters
        ----------
        kind : SimPhaseKind
            Current simulation phase type.
        sim : betse.science.sim.Simulation
            Current simulation.
        cells : betse.science.cells.Cells
            Current cell cluster.
        p : betse.science.parameters.Parameters
            Current simulation configuration.
        '''

        # Avoid circular import dependencies.
        from betse.science.math.cache.cacheabc import SimPhaseCaches

        # Classify all passed parameters.
        self.kind = kind
        self.sim = sim
        self.cells = cells
        self.p = p

        # Classify the cache for this phase.
        self.cache = SimPhaseCaches(self)

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
