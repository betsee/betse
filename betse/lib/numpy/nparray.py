#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level Numpy array and matrix facilities.
'''

#FIXME: Donate the write_csv() function back to Numpy as a new np.savecsv()
#function paralleling the existing np.savetxt() function.

# ....................{ IMPORTS                            }....................
import numpy as np
from betse.exceptions import BetseSequenceException, BetseStrException
from betse.util.path import dirs
from betse.util.type import sequences, types
from betse.util.type.text import strs
from betse.util.type.types import type_check, ClassType, IterableTypes
from collections import OrderedDict
from numpy import ndarray

# ....................{ CONVERTERS                         }....................
@type_check
def from_iterable(iterable: IterableTypes) -> ndarray:
    '''
    Convert the passed iterable into a Numpy array.

    If this iterable is:

    * A non-sequence (e.g., generator), this iterable is first converted into a
      sequence and finally into a Numpy array.
    * A non-Numpy sequence (e.g., :class:`list`), this sequence is directly
      converted into a Numpy array.
    * A Numpy array, this array is returned unmodified.

    Caveats
    ----------
    **This high-level function should always be called in lieue of the low-level
    :func:`np.asarray` function,** which is fundamentally unsafe and should
    *never* be called directly. Unlike this function, the :func:`np.asarray`
    function unsafely converts *any* arbitrary non-sequence into a Numpy array.
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
    iterable : IterableTypes
        Iterable to be converted into a Numpy array.

    Returns
    ----------
    ndarray
        Numpy array converted from this iterable.
    '''

    # Sequence converted from this iterable. If this iterable is already a
    # sequence, reuse this sequence as is; else, convert this iterable into the
    # most space- and time-efficient pure-Python sequence available: a tuple.
    sequence = (
        iterable if sequences.is_sequence(iterable) else tuple(iterable))

    # Numpy array converted from this sequence.
    return np.asarray(sequence)


@type_check
def to_iterable(array: ndarray, cls: ClassType) -> IterableTypes:
    '''
    Convert the passed Numpy array into an iterable of the passed type.

    If the type of this iterable is that of:

    * A non-Numpy iterable (e.g., :class:`list`), this array is first converted
      into a :class:`list` and then into an iterable of this type. To do so,
      this type's ``__init__`` method is expected to accept this :class:`list`
      as a single positional argument.
    * A Numpy array, this array is returned unmodified.

    Parameters
    ----------
    array : ndarray
        Numpy array to be converted.
    cls : ClassType
        Type of the iterable to convert this array into.

    Returns
    ----------
    IterableTypes
        Iterable converted from this Numpy array.
    '''

    # If the iterable to be returned is a Numpy array, return this array as is.
    if cls is ndarray:
        return array

    # List converted from this array.
    array_list = array.tolist()

    # If the iterable to be returned is a list, return this list as is.
    if cls is list:
        return array_list

    # Else, return an iterable converted from this list.
    return cls(array_list)

# ....................{ WRITERS                            }....................
@type_check
def write_csv(filename: str, column_name_to_values: OrderedDict) -> None:
    '''
    Serialize each key-value pair of the passed ordered dictionary into a new
    column in comma-separated value (CSV) format to the plaintext file with the
    passed filename.

    Caveats
    ----------
    To ensure that all columns of this file have the same number of rows, all
    values of this dictionary *must* be one-dimensional sequences of:

    * The same length. If this is not the case, an exception is raised.
    * Any type satisfying the :class:`SequenceTypes` API, including:
      * Numpy arrays.
      * Lists.
      * Tuples.

    Typically, each such value is a one-dimensional Numpy array of floats.

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
          one-dimensional sequence of all arbitrary data comprising this
          column).

    Raises
    ----------
    BetseSequenceException
        If one or more values of this dictionary are either:
        * *Not* sequences.
        * Sequences whose length differs from that of any preceding value
          sequences of this dictionary.
    BetseStrException
        If this column name contains one or more characters reserved for use by
        the CSV non-standard, including:
        * Double quotes, reserved for use as the CSV quoting character.
        * Newlines, reserved for use as the CSV row delimiting character.
    '''

    # Validate the contents of this dictionary. While the np.column_stack()
    # function called below also does, the exceptions raised by the latter are
    # both ambiguous and non-human-readable and hence effectively useless.
    #
    # Length of all prior columns or None if no columns have yet to be iterated.
    columns_prior_len = None

    # List of all column names sanitized such that each name containing one or
    # more comma characters is double-quoted.
    column_names = []

    # For each passed column...
    for column_name, column_values in column_name_to_values.items():
        # If this column is *NOT* a sequence, raise a human-readable exception.
        if not types.is_sequence_nonstr(column_values):
            raise BetseSequenceException(
                'Column "{}" type {!r} not a sequence.'.format(
                    column_name, type(column_values)))

        # Length of this column.
        column_len = len(column_values)

        # If this is the first column to be iterated, require all subsequent
        # columns be of the same length.
        if columns_prior_len is None:
            columns_prior_len = column_len
        # Else if this column's length differs from that of all prior columns,
        # raise a human-readable exception.
        elif column_len != columns_prior_len:
            raise BetseSequenceException(
                'Column "{}" length {} differs from '
                'length {} of prior columns.'.format(
                    column_name, column_len, columns_prior_len))

        # If this column name contains one or more reserved characters, raise an
        # exception. This includes:
        #
        # * Double quotes, reserved for use as the CSV quoting character.
        # * Newlines, reserved for use as the CSV row delimiting character.
        if '"' in column_name:
            raise BetseStrException(
                'Column name {} contains '
                "one or more reserved '\"' characters.".format(column_name))
        if '\n' in column_name:
            raise BetseStrException(
                'Column name {} contains '
                'one or more newline characters.'.format(column_name))

        # If this column name contains one or more commas (reserved for use as
        # the CSV delimiter), double-quote this name. Since the prior logic
        # guarantees this name to *NOT* contain double quotes, no further logic
        # is required
        if ',' in column_name:
            column_name = '"{}"'.format(column_name)

        # Append this sanitized column name to this list of such names.
        column_names.append(column_name)

    # Comma-separated string listing all column names.
    columns_name = strs.join_on(column_names, delimiter=',')

    # Two-dimensional Numpy array of all row data converted from this column
    # data, whose:
    #
    # * First dimension indexes each sampled time step such that each element is
    #   a one-dimensional Numpy array of length the number of columns (i.e., the
    #   number of key-value pairs in the passed dictionary).
    # * Second dimension indexes each column data point for this time step
    columns_values = np.column_stack(column_name_to_values.values())

    # Create the directory containing this file if needed.
    dirs.make_parent_unless_dir(filename)

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
