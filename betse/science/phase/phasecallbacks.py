#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
**Simulation phase callback** (i.e., caller-defined object whose methods are
periodically called while simulating one or more phases) class hierarchy.
'''

# ....................{ IMPORTS                           }....................
from betse.util.io.log import logs
from betse.util.py import pys
from betse.util.type.call.callbacks import CallbacksABC
from betse.util.type.types import NoneType  #type_check

# ....................{ SUPERCLASSES                      }....................
# This subclass is currently an empty placeholder but will be subsequently
# extended with subclass-specific behaviour.
class SimCallbacksABC(CallbacksABC):
    '''
    Abstract base class of all **simulation phase callbacks** (i.e.,
    caller-defined object whose methods are periodically called while
    simulating one or more simulation phases) subclasses.
    '''

    pass

# ....................{ TYPES                             }....................
SimCallbacksABCOrNoneTypes = (SimCallbacksABC, NoneType)
'''
Tuple of both the simulation phase callbacks type *and* that of the ``None``
singleton.
'''

# ....................{ SUBCLASSES                        }....................
class SimCallbacksNoop(SimCallbacksABC):
    '''
    **Noop simulation phase callbacks** (i.e., simulation phase callbacks whose
    methods all silently reduce to noops for efficiency).

    This callbacks subclass is typically instantiated as a convenience fallback
    in the event that no other callbacks subclass is passed by an external
    caller to a simulation subcommand. For efficiency, all callbacks
    reimplemented by this subclass intentionally:

    * Do *not* type-check any parameters passed by the source callable.
    * Do *not* call their superclass implementations.

    See Also
    ----------
    :func:`make_default`
        Factory function internally instantiating and returning a new instance
        of this subclass, which should typically be called in lieu of manually
        instantiating instances of this subclass. This function handles edge
        cases intentionally excluded by this subclass, including the need to
        instantiate an instance of a different subclass when running tests.
    '''

    # ..................{ CALLBACKS                         }..................
    def progress_ranged(
        self, progress_max: int, progress_min: int = 0) -> None:
        pass

    def progressed(self, progress: int) -> None:
        pass

    def progressed_next(self) -> None:
        pass

# ....................{ MAKERS                            }....................
def make_default() -> SimCallbacksABC:
    '''
    Create and return a new simulation phase callbacks object suitable for use
    as an efficient fallback in the event that a caller fails to supply a
    caller-specific callbacks object.

    Specifically, this factory function:

    * If tests are currently being run, this function creates and returns an
      instance of the :class:`SimCallbacksABC` subclass. While less efficient
      than the comparable :class:`SimCallbacksNoop` subclass, doing so ensures
      that contractual guarantees maintained by the :class:`CallbacksABC`
      superclass are properly exercised.
    * Else, this function creates and returns an instance of the
      :class:`SimCallbacksNoop` subclass. While more efficient than the
      comparable :class:`SimCallbacksABC` subclass, doing so avoids contractual
      guarantees maintained by the :class:`CallbacksABC` superclass.
    '''

    # If tests are currently being run...
    if pys.is_testing():
        # Log this behaviour.
        logs.log_debug('Enabling test-specific callbacks...')

        # Return a new simulation phase callbacks object preserving (and hence
        # testing) superclass API guarantees.
        return SimCallbacksABC()
    # Else, tests are *NOT* currently being run. In this case...
    else:
        # Log this behaviour.
        logs.log_debug('Enabling optimized noop callbacks...')

        # Return a new simulation phase callbacks object maximizing efficiency.
        return SimCallbacksNoop()
