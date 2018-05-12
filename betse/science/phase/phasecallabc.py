#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
**Simulation phase callback** (i.e., caller-defined object whose methods are
periodically called while simulating one or more phases) class hierarchy.
'''

# ....................{ IMPORTS                            }....................
# from abc import ABCMeta, abstractmethod
# from betse.util.io.log import logs
from betse.util.type.call.callbacks import CallbacksABC
from betse.util.type.types import NoneType  #type_check

# ....................{ SUPERCLASSES                       }....................
# This subclass is currently an empty placeholder but will be subsequently
# extended with subclass-specific behaviour.
class SimCallbacksABC(CallbacksABC):
    '''
    Abstract base class of all **simulation phase callbacks** (i.e.,
    caller-defined object whose methods are periodically called while
    simulating one or more simulation phases) subclasses.
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
    methods all silently reduce to noops for efficiency).

    This callbacks subclass is typically instantiated as a convenience fallback
    in the event that no other callbacks subclass is passed by an external
    caller to a simulation subcommand.

    Design
    ----------
    For efficiency, all callbacks reimplemented by this subclass intentionally:

    * Do *not* type-check any parameters passed by the source callable.
    * Do *not* call their superclass implementations.
    '''

    # ..................{ CALLBACKS                          }..................
    def progress_ranged(self, progress_min: int, progress_max: int) -> None:
        pass

    def progressed(self, progress: int) -> None:
        pass

    def progressed_next(self) -> None:
        pass
