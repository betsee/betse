#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level **type testers** (i.e., functions testing the types of passed
objects).
'''
# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To avoid non-halting recursive imports when imported at the top-level
# of other modules in the `betse.util` package, this module may import *ONLY*
# from stock Python packages. (By definition, this excludes all BETSE packages.)
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# ....................{ IMPORTS                            }....................
from collections.abc import Iterable, Sequence

# ....................{ TESTERS                            }....................
def is_string(obj: object) -> bool:
    '''
    `True` if the passed object is a **string** (i.e., instance of either the
    `str` class _or_ a subclass thereof).
    '''
    return isinstance(obj, str)

def is_string_nonempty(obj: object) -> bool:
    '''
    `True` if the passed object is a **nonempty string* (i.e., string comprising
    one or more characters and hence _not_ the empty string).
    '''
    return is_string(obj) and len(obj)

# ....................{ TESTERS ~ iterable                 }....................
def is_iterable(obj: object) -> bool:
    '''
    `True` if the passed object is an **iterable**.

    Iterables are objects capable of returning their members one at a time.
    Equivalently, iterables implement the abstract base class
    `collections.Iterable` and hence define the `__iter__()` method.
    '''
    return isinstance(obj, Iterable)

def is_iterable_nonstring(obj: object) -> bool:
    '''
    `True` if the passed object is a **non-string iterable** (i.e., implements
    the abstract base class `collections.Iterable` _and_ is not a string).
    '''
    return is_iterable(obj) and not is_string(obj)

# ....................{ TESTERS ~ sequence                 }....................
def is_sequence(obj: object) -> bool:
    '''
    `True` if the passed object is a **sequence**.

    Sequences are iterables supporting efficient element access via integer
    indices. Equivalently, sequences implement the abstract base class
    `collections.Sequence` and hence define the `__getitem__()` and `__len__()`
    methods (among numerous others).

    While all sequences are iterables, not all iterables are sequences.
    Generally speaking, sequences correspond to the proper subset of iterables
    whose elements are ordered. `dict` and `OrderedDict` are the canonical
    examples. `dict` implements `collections.Iterable` but _not_
    `collections.Sequence`, due to _not_ supporting integer index-based lookup;
    `OrderedDict` implements both, due to supporting such lookup.
    '''
    return isinstance(obj, Sequence)

def is_sequence_nonstring(obj: object) -> bool:
    '''
    `True` if the passed object is a **non-string sequence** (i.e., implements
    the abstract base class `collections.Sequence` _and_ is not a string).
    '''
    return is_sequence(obj) and not is_string(obj)

def is_sequence_nonstring_nonempty(obj: object) -> bool:
    '''
    `True` if the passed object is a **nonempty non-string sequence** (i.e.,
    implements the abstract base class `collections.Sequence`, is not a string,
    and contains at least one element).
    '''
    return is_sequence_nonstring(obj) and len(obj)

# ....................{ ASSERTERS                          }....................
def assert_is_not_string(obj: object) -> str:
    '''
    String asserting the passed object to _not_ be a string.
    '''
    return '"{}" not a string.'.format(obj)

def assert_is_not_string_nonempty(obj: object, label: str) -> str:
    '''
    String asserting the passed object categorized by the passed human-readable
    label to _not_ be a nonempty string.
    '''
    return assert_is_not_string(obj) if not is_string(obj) else\
        '{} empty.'.format(label.capitalize())

# ....................{ ASSERTERS                          }....................
def assert_is_not_sequence_nonstring(obj: object) -> str:
    '''
    String asserting the passed object to _not_ be a non-string sequence.
    '''
    return '"{}" not a non-string sequence (e.g., list, tuple).'.format(obj)

def assert_is_not_sequence_nonstring_nonempty(obj: object, label: str) -> str:
    '''
    String asserting the passed object categorized by the passed human-readable
    label to _not_ be a nonempty non-string sequence.
    '''
    return assert_is_not_sequence_nonstring(obj) if not is_sequence_nonstring(
        obj) else '{} empty.'.format(label.capitalize())

# --------------------( WASTELANDS                         )--------------------
    # if not is_string(obj):
    #     return assert_is_not_string(obj)
    # else:
    #     return '{} empty.'.format(label.capitalize())

