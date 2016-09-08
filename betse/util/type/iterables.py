#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level **non-string iterable** (i.e., non-string object implementing the
abstract base class `collections.abc.Iterable`) facilities.

See Also
----------
:func:`betse.util.type.types.is_iterable`
    Further details on what constitutes iterables and non-string iterables.
'''

# ....................{ IMPORTS                            }....................
import itertools
from betse.exceptions import BetseIterableException
from betse.util.type import types
from betse.util.type.types import (
    type_check, GeneratorType, IterableTypes, SizedType,)
from collections import deque

# ....................{ CLASSES                            }....................
class Sentinel(object):
    '''
    Class encapsulating sentinel objects of arbitrary (albeit unique) value.

    Instances of this class are intended to be used as placeholder objects in
    iterables, typically to identify erroneous or edge-case algorithm input.
    '''

    def __repr__(self) -> str:
        '''
        Human- and machine-readable representation of this sentinel.

        This method has been overridden purely to improve the debuggability of
        algorithms requiring instances of this class.
        '''

        return 'Sentinel()'

# ....................{ CONSTANTS                          }....................
SENTINEL = Sentinel()
'''
Sentinel object of arbitrary value.

This object is internally leveraged by various utility functions (e.g.,
:func:`zip_isometric`) to identify erroneous and edge-case input (e.g.,
iterables of insufficient length).
'''

# ....................{ GETTERS                            }....................
@type_check
def consume(iterable: IterableTypes, iterations: int) -> object:
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
def exhaust(iterable: IterableTypes) -> object:
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
    if iterable_deque:
        return iterable_deque[0]
    # Else, this iterable was already exhausted. Return nothing.
    else:
        return None

# ....................{ SORTERS                            }....................
@type_check
def sort_lexicographic_ascending(iterable: IterableTypes) -> IterableTypes:
    '''
    Sort the passed non-string iterable into a new non-string iterable in
    **ascending lexicographic order** (i.e., traditional order of dead-tree
    dictionaries and encyclopedias).

    Parameters
    ----------
    iterable : Iterable
        Unsorted iterable to be returned sorted. For generality, this iterable
        is _not_ modified by this function.

    Returns
    ----------
    IterableTypes
        Iterable sorted from the passed iterable. For efficiency, this iterable
        is only a shallow rather than deep copy of the passed iterable.
    '''

    return sorted(iterable)   # Well, that was easy.

# ....................{ ZIPPERS                            }....................
#FIXME: Unit test us up.
@type_check
def zip_isometric(*iterables: IterableTypes) -> GeneratorType:
    '''
    Generator zipping the passed iterables of the same length.

    Specifically, this generator iteratively yields an `n`-tuple, where:

    * `n` is the length of each passed iterable.
    * The `i`-th element of this tuple is in the `i`-th passed iterable.

    Parameters
    ----------
    iterables : Iterable
        Tuple of iterables of the same length to be zipped.

    Returns
    ----------
    GeneratorType
        Generator zipping these iterables.

    Raises
    ----------
    :exc:`BetseIterableException`
        If any passed iterable differs in length from any other passed iterable.

    See Also
    ----------
    https://stackoverflow.com/a/32954700/2809027
        Stackoverflow answer strongly inspiring this implementation.
    '''

    # For efficiency, declare this frequently accessed variable to be global,
    # preventing Python from continually attempting to access this variable as a
    # local in the loop below. (This has been timed.)
    global SENTINEL

    # Iteratively zip and yield each n-tuple from the passed n iterables. To
    # efficiently detect iterables of insufficient length, the C-based
    # zip_longest() function is called to fill all iterables of insufficient
    # length with sentinel objects to the expected length.
    #
    # After zipping but before yielding each n-tuple, this n-tuple is then
    # manually searched for sentinel objects. Since these objects may reside at
    # any index of this n-tuple, the entire n-tuple *MUST* be searched in an
    # O(n) manner. While unfortunate, this approach remains substantially more
    # efficient than all alternatives -- largely due to the efficacy of the
    # zip_longest() function. Surprisingly, this "for" loop-based approach has
    # been timed to be faster than the following generator expression:
    #
    #     return (
    #         ntuple if SENTINEL not in ntuple else (
    #             _zip_isometric_error(iterables, ntuple))
    #         for ntuple in zip_longest(
    #             *iterables, fillvalue=SENTINEL)
    #     )
    for ntuple in itertools.zip_longest(
        *iterables, fillvalue=SENTINEL):
        # If this n-tuple contains a sentinel, at least one passed iterable is
        # of insufficient length. Raise a human-readable exception indicating
        # the index and contents of this iterable.
        if SENTINEL in ntuple:
            raise _zip_isometric_error(iterables, ntuple)

        # Else, this n-tuple is valid. Yield it up!
        yield ntuple


def _zip_isometric_error(iterables: tuple, ntuple: tuple) -> None:
    '''
    Raise an exception indicating that the iterable of the passed tuple of
    iterables identified by the passed `n`-tuple is of smaller length than other
    iterables in this tuple of iterables.

    This private function is _only_ intended to be called by the
    :func:`zip_isometric` function.
    '''

    # Index of the erroneously short iterable in this tuple of iterables,
    # identical to the index of the first sentinel in the passed zipped tuple.
    for iterable_short_index, item in enumerate(ntuple):
        if item is SENTINEL:
            break

    # This erroneously short iterable.
    iterable_short = iterables[iterable_short_index]
    # print('\n!!!!faulty n-tuple: {}'.format(ntuple))
    # print('!!!!faulty iterable index: {}'.format(iterable_short_index))

    # If the first iterable is of predefined length (e.g., is *NOT* a generator
    # of dynamic length), end this exception message with this length.
    if isinstance(iterables[0], SizedType):
        exception_suffix = 'length {} of prior iterables'.format(
            len(iterables[0]))
    # Else, end this exception message as is.
    else:
        exception_suffix = 'length of prior iterables'

    # If this erroneously short iterable is of predefined length (e.g., is *NOT*
    # a generator of dynamic length), begin this exception message with this
    # length and end this message with this iterable's contents.
    if isinstance(iterable_short, SizedType):
        exception_prefix = "Length {} of iterable {}".format(
            len(iterable_short), iterable_short_index)
        exception_suffix += ': {!r}'.format(iterable_short)
    # Else, begin this exception message as is.
    else:
        exception_prefix = "Length of iterable {}".format(iterable_short_index)

    # Raise this exception.
    raise BetseIterableException('{} differs from {}'.format(
        exception_prefix, exception_suffix))
