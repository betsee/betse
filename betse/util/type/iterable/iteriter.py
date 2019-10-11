#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **non-string iterable iterators** (i.e., utility functions converting
non-string objects implementing the abstract :class:`collections.abc.Iterable`
base class into similar objects for iteration purposes) facilities.
'''

# ....................{ IMPORTS                           }....................
import itertools
# from betse.util.io.log import logs
from betse.util.type.types import type_check, IterableTypes

# ....................{ ITERATORS ~ combinations          }....................
@type_check
def iter_pairs(iterable: IterableTypes) -> IterableTypes:
    '''
    Iterable of all possible pairs of items in the passed iterable.

    Parameters
    ----------
    iterable : IterableTypes
        Iterable to be inspected.

    Returns
    ----------
    IterableTypes
        Iterable of all possible pairs of items in the passed iterable.

    See Also
    ----------
    :func:`iter_combinations`
        Further details.
    '''

    # Before the one-liner, there was only shared misery.
    return iter_combinations(iterable=iterable, k=2)


@type_check
def iter_combinations(iterable: IterableTypes, k: int) -> IterableTypes:
    '''
    Iterable of all possible **k-combinations** (i.e., subiterables of ``k``
    distinct items) of the passed iterable.

    The iterable created and returned by this function is guaranteed to contain
    exactly ``(n k)`` tuples of ``k`` items each, where:

    * ``n`` is the length of the passed iterable.
    * ``k`` is the passed parameter ``k``.
    * ``(n k)`` is the binomial coefficient -- a non-negative integer either:

      * If ``k > n``, 0.
      * If ``k <= n``, ``n!/[k!(n-k)!]``.

    Moreover, the order of both the k-combinations of this iterable *and* the
    items of each such k-combination is guaranteed to be the exact same as that
    of the passed iterable. Ergo, passing a sorted iterable returns a similarly
    sorted iterable.

    Motivation
    ----------
    This function is a thin wrapper about the standard C-based
    :func:`itertools.combinations` function, whose official documentation and
    examples decidedly leave something to be desired.

    Parameters
    ----------
    iterable : IterableTypes
        Iterable to be inspected.
    k : int
        Length of each k-combination in the iterable to be returned. If this
        length is negative, an exception is raised.

    Returns
    ----------
    IterableTypes
        Iterable of all possible k-combinations of the passed iterable.

    Raises
    ----------
    ValueError
        If ``k < 0``.

    Examples
    ----------
        >>> from betse.util.type.iterable import iteriter
        >>> letter_to_woodburn_harris = (
        ...     'I love to dream, ',
        ...     'but I never try to dream ',
        ...     'and think at the same time.')
        >>> a_lover_of_dreams = iteriter.iter_combinations(
        ...     iterable=letter_to_woodburn_harris, k=2)
        >>> tuple(a_lover_of_dreams)
        (('I love to dream, ', 'but I never try to dream '),
         ('I love to dream, ', 'and think at the same time.'),
         ('but I never try to dream ', 'and think at the same time.'))
    '''

    # And then there was one... liners.
    return itertools.combinations(iterable, k)
