#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **set** (i.e., :class:`set`-like types or instances) functionality.
'''

# ....................{ IMPORTS                           }....................
from betse.util.type.types import type_check, IterableTypes, SetType

# ....................{ MAKERS                            }....................
@type_check
def make_union(*iterables: IterableTypes) -> SetType:
    '''
    **Union** (i.e., set of all objects contained in all passed iterables) of
    all passed iterables.

    This function generalizes the :meth:`set.union` method, which implicitly
    requires the first iterable to be a set, to arbitrary iterables, none of
    which (including the first iterable) are required to be a set.

    Parameters
    ----------
    iterables : Tuple[IterableTypes]
        Tuple of all iterables to be united.

    Raises
    ----------
    BetseSequenceException
        If less than two iterables are passed.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.iterable import generators, sequences

    # If passed no iterables, raise an exception.
    sequences.die_if_empty(sequence=iterables)
    # Else, at least one iterable is passed.

    # If passed only one generator, implicitly expand this generator into the
    # sequence of all items yielded by this generator.
    if len(iterables) == 1 and generators.is_generator(iterables[0]):
        iterables = tuple(iterables[0])

    # First iterable.
    iterable_first = iterables[0]

    # First iterable coerced into a set if nor already a set *OR* preserved as
    # is otherwise.
    #
    # Note that testing whether this iterable is a "SetType" does *NOT*
    # suffice, as that abstract base class sadly fails to require that
    # subclasses adhere to the public API of the "set" builtin (e.g., by
    # defining a union() method).
    iterable_first_set = (
        iterable_first if isinstance(iterable_first, set) else
        set(iterable_first))

    # If passed only one iterable, simply return this set.
    if len(iterables) == 1:
        return iterable_first_set
    # Else, two or more iterables are passed.

    # Generator comprehension yielding each subsequent iterable.
    iterables_rest = (iterable for iterable in iterables[1:])

    # Union of all passed iterables.
    iterables_united = iterable_first_set.union(*iterables_rest)

    # Return this union.
    return iterables_united
