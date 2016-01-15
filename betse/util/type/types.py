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
# WARNING: To raise human-readable exceptions on missing mandatory dependencies
# *AND* avoid non-halting recursive imports when imported at the top-level
# of other modules in the `betse.util` package, this module may import *ONLY*
# from stock Python packages. (By definition, this excludes both BETSE and
# third-party packages.)
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

import re
from collections.abc import Iterable, Mapping, Sequence
from enum import Enum

# ....................{ FORMATTER                          }....................
# This string-centric function is defined in this module rather than in the
# arguably more appropriate "strs" module to drastically simplify assertion
# handlers defined below.
def trim(obj: object) -> str:
    '''
    Convert the passed object to a terse human-readable string suitable for use
    in end user messages.
    '''

    # Maximum length of the string to be returned, defined to be the customary
    # line length of 80 characters less output indentation of 4 characters.
    MAX_LEN = 76

    # Uncompiled regular expression grouping zero or more trailing delimiters.
    # For usability, the expression *ALWAYS* captures exactly one group.
    LAST_DELIMETERS_REGEX = r'.*?([\])}>\'"]*)$'

    # Uncompiled regular expression grouping zero or more leading characters
    # preceding a newline *AND* zero or more trailing delimiters.
    PRENEWLINE_CHARS_LAST_DELIMITERS_REGEX = r'^(.*?)\n{}'.format(
        LAST_DELIMETERS_REGEX)

    # String describing the passed object. For debuggability, the verbose
    # (albeit less human-readable) output of repr() is preferred to the terse
    # (albeit more human-readable) output of str().
    obj_synopsis = repr(obj)

    # If this string contains one or more newlines, replace the substring in
    # this string from the first newline to the string end (excluding any
    # optional trailing delimiters) with an ellipses.
    if '\n' in obj_synopsis:
        obj_synopsis = re.sub(
            PRENEWLINE_CHARS_LAST_DELIMITERS_REGEX,
            r'\1...\2',
            obj_synopsis,
            flags=re.DOTALL
        )

    # If this string still exceeds the maximum, replace this string again.
    if len(obj_synopsis) > MAX_LEN:
        # Optional trailing delimiters if any.
        last_delimiters = re.match(
            LAST_DELIMETERS_REGEX, flags=re.DOTALL).groups()[0]
        last_delimiters_len = len(last_delimiters)

        # If these delimiters exist, remove them from this string.
        if last_delimiters_len:
            obj_synopsis = obj_synopsis[-last_delimiters_len:]

        # Truncate this string to the first leading characters less the length
        # of these delimiters, which will be reappended below.
        obj_synopsis = obj_synopsis[:MAX_LEN - last_delimiters_len]

        # Reappend these delimiters.
        obj_synopsis += last_delimiters

    # Get this synopsis.
    return obj_synopsis

# ....................{ TESTERS                            }....................
def is_bool(obj: object) -> bool:
    '''
    `True` if the passed object is **boolean** (i.e., either `True` or `False`).
    '''
    return isinstance(obj, bool)


def is_char(obj: object) -> bool:
    '''
    `True` if the passed object is a **character** (i.e., a string of length 1).
    '''
    return is_str(obj) and len(obj) == 1


def is_mapping(obj: object) -> bool:
    '''
    `True` if the passed object is a **mapping** (i.e., indexable by strings).

    The canonical examples are `dict` and `OrderedDict` instances.
    '''
    return isinstance(obj, Mapping)

# ....................{ TESTERS ~ iterable                 }....................
def is_iterable(obj: object) -> bool:
    '''
    `True` if the passed object is an **iterable**.

    Iterables are objects capable of returning their members one at a time.
    Equivalently, iterables implement the abstract base class
    `collections.Iterable` and hence define the `__iter__()` method.
    '''
    return isinstance(obj, Iterable)


def is_iterable_nonstr(obj: object) -> bool:
    '''
    `True` if the passed object is a **non-string iterable** (i.e., implements
    the abstract base class `collections.Iterable` _and_ is not a string).
    '''
    return is_iterable(obj) and not is_str(obj)

# ....................{ TESTERS ~ lib                      }....................
def is_matplotlib_colormap(obj: object) -> bool:
    '''
    `True` if the passed object is a Matplotlib colormap.
    '''
    # Avoid importing third-party packages at the top level.
    from matplotlib.colors import Colormap
    return isinstance(obj, Colormap)

# ....................{ TESTERS ~ numeric                  }....................
def is_int(obj: object) -> bool:
    '''
    `True` if the passed object is an integer.
    '''
    return isinstance(obj, int)


def is_int_ge(obj: object, ge: int) -> bool:
    '''
    `True` if the passed object is an integer greater than or equal to the
    second passed integer.
    '''
    assert is_int(ge), assert_not_int(ge)
    return is_int(obj) and obj >= ge


def is_numeric(obj: object) -> bool:
    '''
    `True` if the passed object is **numeric** (i.e., instance of either the
    `int` or `float` classes).
    '''
    return isinstance(obj, (int, float))

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


def is_sequence_nonstr(obj: object) -> bool:
    '''
    `True` if the passed object is a **non-string sequence** (i.e., implements
    the abstract base class `collections.Sequence` _and_ is not a string).

    For generality, this functions returns `True` for both pure-Python
    non-string sequences _and_ non-Pythonic Fortran-based `numpy` arrays and
    matrices (which fail to subclass the `collections.abc.Sequence` API despite
    implementing all methods defined by that subclass).
    '''

    # Avoid importing third-party packages at the top level.
    import numpy

    # Let's do this.
    return (
        # Is this a pure-Python non-string sequence?
        (is_sequence(obj) and not is_str(obj)) or
        # Is this a non-Pythonic Fortran-based numpy array or matrix, all of
        # which subclass the "ndarray" superclass?
        isinstance(obj, numpy.ndarray)
    )


def is_sequence_nonstr_nonempty(obj: object) -> bool:
    '''
    `True` if the passed object is a **nonempty non-string sequence** (i.e.,
    implements the abstract base class `collections.Sequence`, is not a string,
    and contains at least one element).
    '''
    return is_sequence_nonstr(obj) and len(obj)

# ....................{ TESTERS ~ science                  }....................
def is_cells(obj: object) -> bool:
    '''
    `True` if the passed object is an instance of the BETSE-specific `Cells`
    class.
    '''
    # Avoid circular import dependencies.
    from betse.science.cells import Cells
    return isinstance(obj, Cells)


def is_parameters(obj: object) -> bool:
    '''
    `True` if the passed object is an instance of the BETSE-specific
    `Parameters` class.
    '''
    # Avoid circular import dependencies.
    from betse.science.parameters import Parameters
    return isinstance(obj, Parameters)


def is_simulator(obj: object) -> bool:
    '''
    `True` if the passed object is an instance of the BETSE-specific `Simulator`
    class.
    '''
    # Avoid circular import dependencies.
    from betse.science.sim import Simulator
    return isinstance(obj, Simulator)


def is_tissue_picker(obj: object) -> bool:
    '''
    `True` if the passed object is an instance of the BETSE-specific
    `TissuePicker` class.
    '''
    # Avoid circular import dependencies.
    from betse.science.tissue.picker import TissuePicker
    return isinstance(obj, TissuePicker)

# ....................{ TESTERS ~ str                      }....................
def is_str(obj: object) -> bool:
    '''
    `True` if the passed object is a **string** (i.e., instance of the `str`
    class).
    '''
    return isinstance(obj, str)


def is_str_nonempty(obj: object) -> bool:
    '''
    `True` if the passed object is a **nonempty string* (i.e., string comprising
    one or more characters and hence _not_ the empty string).
    '''
    return is_str(obj) and len(obj)

# ....................{ ASSERTERS                          }....................
def assert_not_bool(obj: object) -> str:
    '''
    String asserting the passed object to _not_ be boolean.
    '''
    return '"{}" not boolean (i.e., neither "True" nor "False").'.format(
        trim(obj))


def assert_not_char(obj: object) -> str:
    '''
    String asserting the passed object to _not_ be character.
    '''
    return '"{}" not a character (i.e., string of length 1).'.format(trim(obj))


def assert_not_mapping(obj: object) -> str:
    '''
    String asserting the passed object to _not_ be a mapping.
    '''
    return '"{}" not a mapping (e.g., "dict", "OrderedDict").'.format(trim(obj))

# ....................{ ASSERTERS ~ contains               }....................
def assert_not_in_enum(obj: object, enum: Enum) -> str:
    '''
    String asserting the passed object to _not_ be in the passed enumeration.
    '''
    return '"{}" not in the "{}" enumeration.'.format(trim(obj), trim(enum))

# ....................{ ASSERTERS ~ iterable               }....................
def assert_not_iterable_nonstr(obj: object) -> str:
    '''
    String asserting the passed object to _not_ be a non-string iterable.
    '''
    return '"{}" not a non-string iterable (e.g., dict, list).'.format(
        trim(obj))


def assert_not_iterable_nonstr_nonempty(obj: object, label: str) -> str:
    '''
    String asserting the passed object categorized by the passed human-readable
    label to _not_ be a nonempty non-string iterable.
    '''
    return assert_not_iterable_nonstr(obj) if not is_iterable_nonstr(
        obj) else '{} empty.'.format(label.capitalize())

# ....................{ ASSERTERS ~ lib                    }....................
def assert_not_matplotlib_colormap(obj: object) -> bool:
    '''
    String asserting the passed object to _not_ be a Matplotlib colormap.
    '''
    return '"{}" not a Matplotlib colormap.'.format(trim(obj))

# ....................{ ASSERTERS ~ numeric                }....................
def assert_not_int(obj: object) -> str:
    '''
    String asserting the passed object to _not_ be an integer.
    '''
    return '"{}" not an integer.'.format(trim(obj))


def assert_not_int_ge(obj: object, ge: int) -> str:
    '''
    String asserting the passed object to _not_ be an integer greater than or
    equal to the second passed integer.
    '''
    return '"{}" not an integer or not >= {}.'.format(trim(obj), ge)


def assert_not_numeric(obj: object) -> str:
    '''
    String asserting the passed object to _not_ be numeric.
    '''
    return '"{}" not numeric (i.e., neither an integer nor float).'.format(
        trim(obj))

# ....................{ ASSERTERS ~ sequence               }....................
def assert_not_sequence_nonstr(obj: object) -> str:
    '''
    String asserting the passed object to _not_ be a non-string sequence.
    '''
    return '"{}" not a non-string sequence (e.g., list).'.format(trim(obj))


def assert_not_sequence_nonstr_nonempty(obj: object, label: str) -> str:
    '''
    String asserting the passed object categorized by the passed human-readable
    label to _not_ be a nonempty non-string sequence.
    '''
    return assert_not_sequence_nonstr(obj) if not is_sequence_nonstr(
        obj) else '{} empty.'.format(label.capitalize())

# ....................{ ASSERTERS ~ science                }....................
def assert_not_cells(obj: object) -> str:
    '''
    String asserting the passed object to _not_ be an instance of the BETSE-
    specific `Cells` class.
    '''
    return '"{}" not a "Cells" instance.'.format(trim(obj))


def assert_not_parameters(obj: object) -> str:
    '''
    String asserting the passed object to _not_ be an instance of the BETSE-
    specific `Parameters` class.
    '''
    return '"{}" not a "Parameters" instance.'.format(trim(obj))


def assert_not_simulator(obj: object) -> str:
    '''
    String asserting the passed object to _not_ be an instance of the BETSE-
    specific `Simulator` class.
    '''
    return '"{}" not a "Simulator" instance.'.format(trim(obj))


def assert_not_tissue_picker(obj: object) -> str:
    '''
    String asserting the passed object to _not_ be an instance of the BETSE-
    specific `TissuePicker` class.
    '''
    return '"{}" not a "TissuePicker" instance.'.format(trim(obj))

# ....................{ ASSERTERS ~ str                    }....................
def assert_not_str(obj: object) -> str:
    '''
    String asserting the passed object to _not_ be a string.
    '''
    return '"{}" not a string.'.format(trim(obj))


def assert_not_str_nonempty(obj: object, label: str) -> str:
    '''
    String asserting the passed object categorized by the passed human-readable
    label to _not_ be a nonempty string.
    '''
    return (assert_not_str(obj) if not is_str(obj) else
        '{} empty.'.format(label.capitalize()))
