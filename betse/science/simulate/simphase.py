#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
High-level **simulation phase** (e.g., cell cluster seeding, initialization, and
simulation) functionality.
'''

# ....................{ IMPORTS                            }....................
from abc import ABCMeta, abstractmethod
from betse.util.py import references
from betse.util.type.enums import EnumOrdered
from betse.util.type.types import type_check

# ....................{ ENUMS                              }....................
SimPhaseType = EnumOrdered('SimPhaseType', ('SEED', 'INIT', 'SIM',))
'''
Ordered enumeration of all possible simulation phases.

Each member of this enumeration is arbitrarily comparable to each other member.
Each member's value is less than that of another member's value if and only if
the former simulation phase is performed _before_ the latter. Specifically,
this enumeration is a total ordering such that:

    >>> SimPhaseType.SEED < SimPhaseType.INIT < SimPhaseType.SIM
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

# ....................{ SUPERCLASSES                       }....................
class SimPhaseABC(object, metaclass=ABCMeta):
    '''
    Abstract base class of all simulation phase subclasses, each encapsulating
    all lower-level objects required for running a single simulation phase in
    a particular manner (e.g., with either strong or weak references).

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
    cells : betse.science.cells.Cells
        Current cell cluster.
    p : betse.science.parameters.Parameters
        Current simulation configuration.
    sim : betse.science.sim.Simulator
        Current simulation.
    '''

    # ..................{ INITIALIZORS                       }..................
    @abstractmethod
    def __init__(
        self,
        sim:   'betse.science.sim.Simulator',
        cells: 'betse.science.cells.Cells',
        p:     'betse.science.parameters.Parameters',
    ) -> None:
        '''
        Initialize this simulation phase instance.

        Parameters
        ----------
        cells : Cells
            Current cell cluster.
        p : Parameters
            Current simulation configuration.
        sim : Simulator
            Current simulation.
        '''

        pass

# ....................{ SUBCLASSES                         }....................
class SimPhaseStrong(SimPhaseABC):
    '''
    Top-level simulation phase, encapsulating all lower-level objects required
    for running a single simulation phase with **strong references** (i.e.,
    references "keeping" their referents alive).

    Caveats
    ----------
    This class is intended to be instantiated *only* by the root-level
    :class:`betse.science.simrunner.SimRunner` class. All other callers should
    instantiate the :class:`SimRunnerWeak` class instead, preventing harmful
    circular object references from being formed. Circular object references
    inhibit garbage collection and hence encourage resource bottlenecks in
    memory allocation for long-running applications (e.g., the BETSE GUI).
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(
        self,
        sim:   'betse.science.sim.Simulator',
        cells: 'betse.science.cells.Cells',
        p:     'betse.science.parameters.Parameters',
    ) -> None:

        # Strongly classify all passed parameters.
        self.cells = cells
        self.p = p
        self.sim = sim


class SimPhaseWeak(SimPhaseABC):
    '''
    Top-level simulation phase, encapsulating all lower-level objects required
    for running a single simulation phase with **weak references** (i.e.,
    references permitting their referents to die at any time).

    Caveats
    ----------
    This class is intended to be instantiated by *all* callers other than the
    root-level :class:`betse.science.simrunner.SimRunner` class. This class
    prevents harmful circular object references from being formed. Circular
    object references inhibit garbage collection and hence encourage resource
    bottlenecks in memory allocation for long-running applications (e.g., the
    BETSE GUI).

    Usage of this class is guaranteed to *always* be safe. The root-level
    :class:`betse.science.simrunner.SimRunner` class owns *all* other
    simulation-centric classes, including directly owning all objects weakly
    referenced by this class and indirectly owning all instances of this class.
    Since the :class:`betse.science.simrunner.SimRunner` class necessarily lives
    longer than this class and thus all weak references internally maintained by
    this class, no complications arise. Such references *always* yield their
    expected referents rather than non-deterministically yielding ``None``
    (e.g., when the referenced objects are unexpectedly garbage-collected).
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(
        self,
        sim:   'betse.science.sim.Simulator',
        cells: 'betse.science.cells.Cells',
        p:     'betse.science.parameters.Parameters',
    ) -> None:

        # Weakly classify all passed parameters.
        self.cells = references.proxy_weak(cells)
        self.p     = references.proxy_weak(p)
        self.sim   = references.proxy_weak(sim)
