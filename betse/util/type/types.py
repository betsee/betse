#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
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
from collections.abc import Container, Iterable, Mapping, Sequence
from enum import Enum, EnumMeta

# ....................{ FORMATTER                          }....................
# This string-centric function is defined in this module rather than in the
# arguably more appropriate "strs" module to drastically simplify assertion
# handlers defined below.
def trim(obj: object) -> str:
    '''
    Convert the passed object to a terse human-readable string suitable for use
    in end-user messages.
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
    `True` only if the passed object is **boolean** (i.e., either `True` or
    `False`).
    '''
    return isinstance(obj, bool)


def is_char(obj: object) -> bool:
    '''
    `True` only if the passed object is a **character** (i.e., string of length
    1).
    '''
    return is_str(obj) and len(obj) == 1


# ....................{ TESTERS ~ callable                 }....................
def is_callable(obj: object) -> bool:
    '''
    `True` only if the passed object is **callable** (e.g., function, method,
    class defining the special `__call__()` method).
    '''
    return callable(obj)

# ....................{ TESTERS ~ collection               }....................
def is_mapping(obj: object) -> bool:
    '''
    `True` only if the passed object is a **mapping** (i.e., indexable by
    strings).

    Customary mappings include `dict` and `OrderedDict` instances.
    '''
    return isinstance(obj, Mapping)

# ....................{ TESTERS ~ collection               }....................
def is_container(obj: object) -> bool:
    '''
    `True` only if the passed object is a **container** (i.e., implements the
    `__contains__` special method returning `True` only if that container
    contains the passed element).

    Most collections of interest (e.g., `dict`, `list`, `set`) are containers.
    '''
    return isinstance(obj, Container)


def is_container_nonstr(obj: object) -> bool:
    '''
    `True` only if the passed object is a **non-string container** (i.e.,
    implements the `__contains__` special method _and_ is not a string).
    '''
    return is_container(obj) and not is_str(obj)

# ....................{ TESTERS ~ collection : iterable    }....................
def is_iterable(obj: object) -> bool:
    '''
    `True` only if the passed object is an **iterable**.

    Iterables are objects capable of returning their members one at a time.
    Equivalently, iterables implement the abstract base class
    `collections.Iterable` and hence define the `__iter__()` method.
    '''
    return isinstance(obj, Iterable)


def is_iterable_nonstr(obj: object) -> bool:
    '''
    `True` only if the passed object is a **non-string iterable** (i.e.,
    implements the abstract base class `collections.Iterable` _and_ is not a
    string).
    '''
    return is_iterable(obj) and not is_str(obj)

# ....................{ TESTERS ~ sequence                 }....................
def is_sequence(obj: object) -> bool:
    '''
    `True` only if the passed object is a **sequence**.

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
    `True` only if the passed object is a **non-string sequence** (i.e.,
    implements the abstract base class `collections.Sequence` _and_ is not a
    string).

    For generality, this functions returns `True` for both pure-Python
    non-string sequences _and_ non-Pythonic Fortran-based `numpy` arrays and
    matrices (which fail to subclass the `collections.abc.Sequence` API despite
    implementing all methods defined by that subclass).
    '''

    # Let's do this.
    return (
        # Is this a pure-Python non-string sequence?
        (is_sequence(obj) and not is_str(obj)) or
        # Is this a non-Pythonic Fortran-based numpy array or matrix?
        is_numpy_array(obj)
    )


def is_sequence_nonstr_nonempty(obj: object) -> bool:
    '''
    `True` only if the passed object is a **nonempty non-string sequence**
    (i.e., implements the abstract base class `collections.Sequence`, is not a
    string, _and_ contains at least one element).
    '''
    return is_sequence_nonstr(obj) and len(obj)

# ....................{ TESTERS ~ enum                     }....................
def is_enum(obj: object) -> bool:
    '''
    `True` only if the passed object is an enumeration.
    '''
    return isinstance(obj, EnumMeta)


def is_in_enum(obj: object, enum: Enum) -> bool:
    '''
    `True` only if the passed object is in the passed enumeration.

    While trivial, this tester is provided for orthogonality with the
    `assert_not_in_enum()` function.
    '''
    assert is_enum(enum), assert_not_enum(enum)
    return obj in enum

# ....................{ TESTERS ~ error                    }....................
def is_exception(obj: object) -> bool:
    '''
    `True` only if the passed object is an **exception**.
    '''
    return isinstance(obj, Exception)

# ....................{ TESTERS ~ lib                      }....................
def is_numpy_array(obj: object) -> bool:
    '''
    `True` only if the passed object is a Numpy array or matrix.

    This function returns true if the passed object is an instance of the
    Numpy-specific `ndarray` superclass.
    '''
    # Avoid importing third-party packages at the top level.
    from numpy import ndarray
    return isinstance(obj, ndarray)

# ....................{ TESTERS ~ lib : matplotlib         }....................
def is_matplotlib_collection(obj: object) -> bool:
    '''
    `True` only if the passed object is a Matplotlib collection.
    '''
    # Avoid importing third-party packages at the top level.
    from matplotlib.collections import Collection
    return isinstance(obj, Collection)


def is_matplotlib_colormap(obj: object) -> bool:
    '''
    `True` only if the passed object is a Matplotlib colormap.
    '''
    # Avoid importing third-party packages at the top level.
    from matplotlib.colors import Colormap
    return isinstance(obj, Colormap)


def is_matplotlib_mappable(obj: object) -> bool:
    '''
    `True` only if the passed object is a Matplotlib mappable.
    '''
    # Avoid importing third-party packages at the top level.
    from matplotlib.cm import ScalarMappable
    return isinstance(obj, ScalarMappable)


def is_matplotlib_polycollection(obj: object) -> bool:
    '''
    `True` only if the passed object is a Matplotlib **polygon collection**
    (i.e., an instance of the `PolyCollection` class).
    '''
    # Avoid importing third-party packages at the top level.
    from matplotlib.collections import PolyCollection
    return isinstance(obj, PolyCollection)


def is_matplotlib_streamplot(obj: object) -> bool:
    '''
    `True` only if the passed object is a Matplotlib **streamplot** (i.e., an
    object returned by the `matplotlib.plt.streamplot()` function).
    '''
    # Avoid importing third-party packages at the top level.
    from matplotlib.streamplot import StreamplotSet
    return isinstance(obj, StreamplotSet)


def is_matplotlib_trimesh(obj: object) -> bool:
    '''
    `True` only if the passed object is a Matplotlib **triangle mesh** (i.e., an
    instance of the `TriMesh` class).
    '''
    # Avoid importing third-party packages at the top level.
    from matplotlib.collections import TriMesh
    return isinstance(obj, TriMesh)

# ....................{ TESTERS ~ lib : py                 }....................
def is_py_path_local(obj: object) -> bool:
    '''
    `True` only if the passed object is a `py.path.local` instance.
    '''
    # Avoid importing third-party packages at the top level.
    from py._path.local import LocalPath
    return isinstance(obj, LocalPath)

# ....................{ TESTERS ~ numeric                  }....................
def is_int(obj: object) -> bool:
    '''
    `True` only if the passed object is an integer.
    '''
    return isinstance(obj, int)


def is_int_ge(obj: object, ge: int) -> bool:
    '''
    `True` only if the passed object is an integer greater than or equal to the
    second passed integer.
    '''
    assert is_int(ge), assert_not_int(ge)
    return is_int(obj) and obj >= ge


def is_numeric(obj: object) -> bool:
    '''
    `True` only if the passed object is **numeric** (i.e., instance of either
    the `int` or `float` types).
    '''
    return isinstance(obj, (int, float))

# ....................{ TESTERS ~ science                  }....................
def is_cells(obj: object) -> bool:
    '''
    `True` only if the passed object is an instance of the BETSE-specific
    `Cells` class.
    '''
    # Avoid circular import dependencies.
    from betse.science.cells import Cells
    return isinstance(obj, Cells)


def is_parameters(obj: object) -> bool:
    '''
    `True` only if the passed object is an instance of the BETSE-specific
    `Parameters` class.
    '''
    # Avoid circular import dependencies.
    from betse.science.parameters import Parameters
    return isinstance(obj, Parameters)


def is_simulator(obj: object) -> bool:
    '''
    `True` only if the passed object is an instance of the BETSE-specific
    `Simulator` class.
    '''
    # Avoid circular import dependencies.
    from betse.science.sim import Simulator
    return isinstance(obj, Simulator)


def is_tissue_picker(obj: object) -> bool:
    '''
    `True` only if the passed object is an instance of the BETSE-specific
    `TissuePicker` class.
    '''
    # Avoid circular import dependencies.
    from betse.science.tissue.picker import TissuePicker
    return isinstance(obj, TissuePicker)

# ....................{ TESTERS ~ str                      }....................
def is_str(obj: object) -> bool:
    '''
    `True` only if the passed object is a **string** (i.e., instance of the
    `str` class).
    '''

    return isinstance(obj, str)


def is_str_or_none(obj: object) -> bool:
    '''
    `True` only if the passed object is either a string _or_ `None`.
    '''

    return isinstance(obj, str) or obj is None


def is_str_nonempty(obj: object) -> bool:
    '''
    `True` only if the passed object is a **nonempty string* (i.e., string
    comprising one or more characters and hence _not_ the empty string).
    '''

    return is_str(obj) and len(obj)


def is_str_nonempty_or_none(obj: object) -> bool:
    '''
    `True` only if the passed object is either a nonempty string _or_ `None`.
    '''

    return is_str_nonempty(obj) or obj is None

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

# ....................{ TESTERS ~ callable                 }....................
def assert_not_callable(obj: object) -> bool:
    '''
    String asserting the passed object to _not_ be callable.
    '''
    return '"{}" not callable.'.format(trim(obj))

# ....................{ ASSERTERS ~ collection             }....................
def assert_not_mapping(obj: object) -> str:
    '''
    String asserting the passed object to _not_ be a mapping.
    '''
    return '"{}" not a mapping (e.g., "dict", "OrderedDict").'.format(trim(obj))

# ....................{ ASSERTERS ~ collection : container }....................
def assert_not_container_nonstr(obj: object) -> str:
    '''
    String asserting the passed object to _not_ be a non-string container.
    '''
    return '"{}" not a non-string container (e.g., "dict", "set").'.format(
        trim(obj))

# ....................{ ASSERTERS ~ collection : iterable  }....................
def assert_not_iterable_nonstr(obj: object) -> str:
    '''
    String asserting the passed object to _not_ be a non-string iterable.
    '''
    return '"{}" not a non-string iterable (e.g., "dict", "list").'.format(
        trim(obj))


def assert_not_iterable_nonstr_nonempty(obj: object, label: str) -> str:
    '''
    String asserting the passed object categorized by the passed human-readable
    label to _not_ be a nonempty non-string iterable.
    '''
    return assert_not_iterable_nonstr(obj) if not is_iterable_nonstr(
        obj) else '{} empty.'.format(label.capitalize())

# ....................{ ASSERTERS ~ collection : sequence  }....................
def assert_not_sequence_nonstr(obj: object) -> str:
    '''
    String asserting the passed object to _not_ be a non-string sequence.
    '''
    return '"{}" not a non-string sequence (e.g., "list").'.format(trim(obj))


def assert_not_sequence_nonstr_nonempty(obj: object, label: str) -> str:
    '''
    String asserting the passed object categorized by the passed human-readable
    label to _not_ be a nonempty non-string sequence.
    '''
    return assert_not_sequence_nonstr(obj) if not is_sequence_nonstr(
        obj) else '{} empty.'.format(label.capitalize())

# ....................{ ASSERTERS ~ enum                   }....................
def assert_not_enum(obj: object) -> str:
    '''
    String asserting the passed object to _not_ be an enumeration.
    '''
    return '"{}" not an enumeration "{}".'.format(trim(obj))


def assert_not_in_enum(obj: object, enum: Enum) -> str:
    '''
    String asserting the passed object to _not_ be in the passed enumeration.
    '''
    return '"{}" not in the enumeration "{}".'.format(trim(obj), trim(enum))

# ....................{ ASSERTERS ~ error                  }....................
def assert_not_exception(obj: object) -> str:
    '''
    String asserting the passed object to _not_ be an exception.
    '''
    return '"{}" not an exception.'.format(trim(obj))

# ....................{ ASSERTERS ~ lib                    }....................
def assert_not_numpy_array(obj: object) -> bool:
    '''
    String asserting the passed object to _not_ be a Numpy array or matrix.
    '''
    return '"{}" not a array or matrix.'.format(trim(obj))

# ....................{ ASSERTERS ~ lib : matplotlib       }....................
def assert_not_matplotlib_collection(obj: object) -> bool:
    '''
    String asserting the passed object to _not_ be a Matplotlib collection.
    '''
    return '"{}" not a Matplotlib collection.'.format(trim(obj))


def assert_not_matplotlib_colormap(obj: object) -> bool:
    '''
    String asserting the passed object to _not_ be a Matplotlib colormap.
    '''
    return '"{}" not a Matplotlib colormap.'.format(trim(obj))


def assert_not_matplotlib_mappable(obj: object) -> bool:
    '''
    String asserting the passed object to _not_ be a Matplotlib mappable.
    '''
    return '"{}" not a Matplotlib mappable.'.format(trim(obj))


def assert_not_matplotlib_polycollection(obj: object) -> bool:
    '''
    String asserting the passed object to _not_ be a Matplotlib polygon
    collection.
    '''
    return '"{}" not a Matplotlib polygon collection.'.format(trim(obj))


def assert_not_matplotlib_streamplot(obj: object) -> bool:
    '''
    String asserting the passed object to _not_ be a Matplotlib streamplot.
    '''
    return '"{}" not a Matplotlib streamplot.'.format(trim(obj))


def assert_not_matplotlib_trimesh(obj: object) -> bool:
    '''
    String asserting the passed object to _not_ be a Matplotlib triangle mesh.
    '''
    return '"{}" not a Matplotlib triangle mesh.'.format(trim(obj))

# ....................{ ASSERTERS ~ lib : py               }....................
def assert_not_py_path_local(obj: object) -> bool:
    '''
    String asserting the passed object to _not_ be a `py.path.local` instance.
    '''
    return '"{}" not a "py.path.local" instance.'.format(trim(obj))

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


def assert_not_str_or_none(obj: object) -> str:
    '''
    String asserting the passed object to be neither a string _nor_ `None`.
    '''

    return '"{}" not a string or "None".'.format(trim(obj))


def assert_not_str_nonempty(obj: object, label: str) -> str:
    '''
    String asserting the passed object categorized by the passed human-readable
    label to _not_ be a nonempty string.
    '''

    if not is_str(obj):
        return assert_not_str(obj)
    else:
        return '{} empty.'.format(label.capitalize())


def assert_not_str_nonempty_or_none(obj: object, label: str) -> str:
    '''
    String asserting the passed object categorized by the passed human-readable
    label to be neither a nonempty string _nor_ `None`.
    '''

    if not is_str_or_none(obj):
        return assert_not_str_or_none(obj)
    else:
        return assert_not_str_nonempty(obj, label)
