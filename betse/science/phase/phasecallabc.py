#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
**Simulation phase callback** (i.e., caller-defined object whose methods are
periodically called while simulating one or more phases) class hierarchy.
'''

# ....................{ IMPORTS                            }....................
from abc import ABCMeta, abstractmethod
from betse.util.io.log import logs
from betse.util.type.types import NoneType  #type_check

# ....................{ SUPERCLASSES                       }....................
class SimCallbacksABC(metaclass=ABCMeta):
    '''
    Abstract base class of all **simulation phase callbacks** (i.e.,
    caller-defined object whose methods are periodically called while simulating
    one or more phases) subclasses.
    '''

    # ..................{ CALLBACKS ~ progress               }..................
    # Subclasses are required to implement the following abstract callbacks.

    @abstractmethod
    def progress_ranged(self, progress_min: int, progress_max: int) -> None:
        '''
        Callback passed the range of progress values subsequently passed to the
        :meth:`progressed` callback by the current simulation subcommand
        (e.g., :meth:`SimRunner.seed`, :meth:`SimRunner.init`).

        Each simulation subcommand calls this callback once *before* simulating
        that simulation phase. This callback is guaranteed to be called *before*
        the :meth:`progressed` callback is called by that subcommand, enabling
        callers to initialize external high-level objects (e.g., progress bars)
        requiring this range.

        Parameters
        ----------
        progress_min : int
            Minimum value subsequently passed to the :meth:`progressed` callback
            by this simulation subcommand. While this value is typically 0,
            callers should *not* necessarily assume this to be the case.
        progress_max : int
            Maximum value subsequently passed to the :meth:`progressed` callback
            by this simulation subcommand.
        '''

        pass


    @abstractmethod
    def progressed(self, progress: int) -> None:
        '''
        Callback passed the range of progress values subsequently passed to the
        :meth:`progressed` callback by the current simulation subcommand
        (e.g., :meth:`SimRunner.seed`, :meth:`SimRunner.init`

        Each simulation subcommand calls this callback repeatedly while
        simulating that simulation phase, enabling callers to incrementally
        update external high-level objects (e.g., progress bars).

        Parameters
        ----------
        progress : int
            Current state of progress for this simulation subcommand. This
            integer is guaranteed to be in the range ``[progress_min,
            progress_max]``, where ``progress_min`` and ``progress_max`` are the
            pair of integers previously passed to the :attr:`progress_ranged`
            callback for this simulation subcommand.
        '''

        pass

# ....................{ TYPES                              }....................
SimCallbacksABCOrNoneTypes = (SimCallbacksABC, NoneType)
'''
Tuple of both the simulation phase callbacks type *and* that of the ``None``
singleton.
'''

# ....................{ SUBCLASSES                         }....................
class SimCallbacksNoop(SimCallbacksABC):
    '''
    **Noop simulation phase callbacks** (i.e., simulation phase callbacks whose
    methods all silently reduce to noops).

    This callbacks subclass is typically instantiated as a convenience fallback
    in the event that no other callbacks subclass is passed by an external
    caller to a simulation subcommand.
    '''

    # ..................{ CALLBACKS                          }..................
    def progress_ranged(self, progress_min: int, progress_max: int) -> None:
        pass

    def progressed(self, progress: int) -> None:
        pass
