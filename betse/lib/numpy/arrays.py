#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level Numpy array and matrix facilities.
'''

# ....................{ IMPORTS                            }....................
import numpy as np
from betse.util.type import strs
from betse.util.type.types import type_check, SequenceTypes
from collections import OrderedDict
from numpy import ndarray

# ....................{ CONVERTERS                         }....................
@type_check
def from_sequence(sequence: SequenceTypes) -> ndarray:
    '''
    Convert the passed sequence into a Numpy array.

    If this sequence is already a Numpy array, this function returns this array
    unmodified; else, this function converts this sequence into a Numpy array
    and returns that array.

    Caveats
    ----------
    **This high-level function should always be called in lieue of the low-level
    :func:`np.asarray` function,** which is fundamentally unsafe and should
    _never_ be called directly. Unlike this function, the :func:`np.asarray`
    function unsafely converts _any_ arbitrary non-sequence into a Numpy array.
    This function safely wraps the unsafe :func:`np.asarray` function with type
    and sanity checking, preventing that function's overly permissive design
    imperatives from corrupting the fragile purity of this codebase: e.g.,

        >>> import numpy as np
        >>> np.asarray(None)
        array(None, dtype=object)
        >>> np.asarray(
        ...     'We are the Bug. '
        ...     'Your computational and technological distinctiveness '
        ...     'will be added to our own. Resistance is futile.')
        array('We are the Bug. Your computational and technological distinctiveness will be added to our own. Resistance is futile.',
        dtype='<U116')

    Parameters
    ----------
    sequence : SequenceTypes
        Sequence to be converted into a Numpy array.

    Returns
    ----------
    ndarray
        Numpy array converted from this sequence.
    '''

    # Thanks to the @type_check decorator invoked above, the passed parameter is
    # guaranteed to be a sequence safely passable as is to this function.
    return np.asarray(sequence)

# ....................{ SAVERS                             }....................
@type_check
def save_csv(filename: str, column_name_to_values: OrderedDict) -> None:
    '''
    Serialize each key-value pair of the passed ordered dictionary into a new
    column in comma-separated value (CSV) format to the plaintext file with the
    passed filename.

    Parameters
    ----------
    filename : str
        Absolute or relative path of the plaintext file to be written. If this
        file already exists, this file is silently overwritten.
    column_name_to_values: OrderedDict
        Ordered dictionary of all columns to be serialized such that:
        * Each key of this dictionary is a **column name** (i.e., terse string
          describing the type of data contained in this column).
        * Each value of this dictionary is **column data** (i.e.,
          one-dimensional sequence of all row data comprising this column).
          Ideally, each such value is a one-dimensional Numpy array of floats.
    '''

    # Comma-separated string listing all column names.
    columns_name = strs.join_on(column_name_to_values.keys(), delimiter=',')

    # Two-dimensional Numpy array of all row data converted from this column
    # data, whose:
    #
    # * First dimension indexes each sampled time step such that each element is
    #   a one-dimensional Numpy array of length the number of columns (i.e., the
    #   number of key-value pairs in the passed dictionary).
    # * Second dimension indexes each column data point for this time step
    columns_values = np.column_stack(column_name_to_values.values())

    # Serialize these sequences to this file in CSV format.
    np.savetxt(
        fname=filename,
        X=columns_values,
        header=columns_name,
        delimiter=',',

        # Prevent Numpy from prefixing the above header by "# ". Most popular
        # software importing CSV files implicitly supports a comma-delimited
        # first line listing all column names.
        comments='',
    )
