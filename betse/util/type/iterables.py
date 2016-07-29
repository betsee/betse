#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level **non-string iterable** (i.e., non-string object implementing the
abstract base class `collections.abc.Iterable`) facilities.

See Also
----------
betse.util.type.types.is_iterable
    Further details on what constitutes iterables and non-string iterables.
'''

# ....................{ IMPORTS                            }....................
from betse.util.type import types
from betse.util.type.types import type_check, IterableType
from collections import deque

# ....................{ GETTERS                            }....................
@type_check
def consume(iterable: IterableType, iterations: int) -> object:
    '''
    Consume the passed number of iterations from the passed iterable by
    advancing this iterable forward by this number of iterations.

    For efficiency, this iterable is consumed at C speeds via standard objects
    implemented in low-level C rather than high-level Python.

    Parameters
    ----------
    iterable : Iterable
        Iterable to be consumed.
    iterations : int
        Number of iterations to advance this iterable. This number should be
        strictly positive (i.e., `iterations >= 1`).

    Returns
    ----------
    object
        Object yielded by the last iterable iteration (i.e., the last `next()`
        method called on this iterable) if any _or_ `None` if this iterable was
        already exhausted (i.e., empty) when passed.

    See Also
    ----------
    http://docs.python.org/3/library/itertools.html
        Official documentation strongly inspiring this function.
    '''
    assert types.is_int_positive(iterations), (
        types.assert_not_int_positive(iterations))

    # Advance to the empty slice starting (and hence immediately ending) at the
    # passed iterable position and return the value of this iterable at this
    # position if any.
    return next(slice(iterable, iterations, iterations), None)


@type_check
def exhaust(iterable: IterableType) -> object:
    '''
    Exhaust the passed iterable by advancing this iterable directly past its
    last iteration.

    For efficiency, this iterable is exhausted at C speeds via standard objects
    implemented in low-level C rather than high-level Python.

    Caveats
    ----------
    **This function should _only_ be called for finite iterables.** If the
    passed iterable:

    * Explicitly halts with a `StopIteration` exception and hence is finite,
      this function exhausts this iterable as expected.
    * Does _not_ explicitly halt with a `StopIteration` exception and hence is
      infinite, this function reduces to an **infinite loop.**

    Parameters
    ----------
    iterable : Iterable
        Iterable to be exhausted.

    Returns
    ----------
    object
        Object yielded by the last iterable iteration (i.e., the last `next()`
        method called on this iterable) if any _or_ `None` if this iterable was
        already exhausted (i.e., empty) when passed.

    See Also
    ----------
    http://docs.python.org/3/library/itertools.html
        Official documentation strongly inspiring this function.
    '''

    # For efficiency, feed this iterable into a zero-length deque retaining only
    # the last iterated value if any.
    iterable_deque = deque(iterable, maxlen=1)

    # If this iterable was *NOT* already exhausted, return its last value.
    # Else, this iterable was already exhausted. Return nothing.
    if iterable_deque:
        return iterable_deque[0]
    else:
        return None
