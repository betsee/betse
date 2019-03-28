#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **Numpy type conversion** (i.e., functions converting Numpy arrays to
and from various types) facilities.
'''

# ....................{ IMPORTS                           }....................
import numpy as np
# from betse.util.io.log import logs
from betse.util.type.iterable import sequences
from betse.util.type.types import (
    type_check, ClassType, IterableTypes, NumpyArrayType,)

# ....................{ GLOBALS                           }....................
DTYPE_UNSIGNED_TO_SIGNED = {
    np.uint0:  np.int8,
    np.uint8:  np.int16,
    np.uint16: np.int32,
    np.uint32: np.int64,
}
'''
Dictionary mapping from all unsigned Numpy data types to the next largest
signed Numpy data types, guaranteed to preserve the contents of arrays whose
data types are these unsigned types in the most space- and time-efficient
manner.

See Also
----------
:func:`to_signed`
    Further discussion.
'''

# ....................{ TESTERS                           }....................
def is_array(obj: object) -> bool:
    '''
    ``True`` only if the passed object is a Numpy array.
    '''

    return isinstance(obj, NumpyArrayType)

# ....................{ GETTERS                           }....................
#FIXME: Shift elsewhere, as this function does no conversions and hence has no
#commonality with the remainder of this submodule -- say, into a new
#"betse.lib.numpy.npindex" submodule.
@type_check
def get_subarray_indices(
    superarray: NumpyArrayType, subarray: NumpyArrayType) -> NumpyArrayType:
    '''
    One-dimensional Numpy array of the indices of the first passed Numpy array
    of all items also found in the second passed Numpy array.

    Parameters
    ----------
    superarray : NumpyArrayType
        Numpy array to return a subset of the indices of.
    subarray : NumpyArrayType
        Numpy array to find the items of in the first passed Numpy array.

    Parameters
    ----------
    NumpyArrayType
        One-dimensional Numpy array of the indices of the superarray of all
        items also found in the subarray.

    See Also
    ----------
    https://stackoverflow.com/a/8251757/2809027
        StackOverflow answer strongly inspiring this implementation.
    '''

    #FIXME: Document this non-trivial behaviour, please.
    superarray_sorted_indices = np.argsort(superarray)
    superarray_sorted = superarray[superarray_sorted_indices]
    superarray_sorted_subarray_indices = np.searchsorted(
        superarray_sorted, subarray)

    superarray_sorted_subarray_extant_indices = np.take(
        superarray_sorted_indices,
        superarray_sorted_subarray_indices,
        mode="clip")
    superarray_subarray_mask = (
        superarray[superarray_sorted_subarray_extant_indices] != subarray)
    superarray_subarray_mask_indices = np.ma.array(
        superarray_sorted_subarray_extant_indices,
        mask=superarray_subarray_mask)

    # Return a non-masked rather than masked array for generality.
    return superarray_subarray_mask_indices.data

# ....................{ CONVERTERS ~ iterable             }....................
@type_check
def from_iterable(iterable: IterableTypes) -> NumpyArrayType:
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
    **This high-level function should always be called in lieue of the
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
        ...     'will be added to our own. Resistance is futile.'
        ... )
        array('We are the Bug. Your computational and technological distinctiveness will be added to our own. Resistance is futile.',
        dtype='<U116')

    Parameters
    ----------
    iterable : IterableTypes
        Iterable to be converted into a Numpy array.

    Returns
    ----------
    NumpyArrayType
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
def to_iterable(array: NumpyArrayType, cls: ClassType) -> IterableTypes:
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
    array : NumpyArrayType
        Numpy array to be converted.
    cls : ClassType
        Type of the iterable to convert this array into.

    Returns
    ----------
    IterableTypes
        Iterable converted from this Numpy array.
    '''

    # If the iterable to be returned is a Numpy array, return this array as is.
    if cls is NumpyArrayType:
        return array

    # List converted from this array.
    array_list = array.tolist()

    # If the iterable to be returned is a list, return this list as is.
    if cls is list:
        return array_list

    # Else, return an iterable converted from this list.
    return cls(array_list)

# ....................{ CONVERTERS ~ signed               }....................
@type_check
def to_signed(array: NumpyArrayType) -> NumpyArrayType:
    '''
    Convert the passed possibly unsigned Numpy array into a signed Numpy array
    exactly containing the same elements.

    If this array's data type is:

    * Already signed (e.g., :attr:`np.int32`, :attr:`np.float64`), this array
      is returned unmodified.
    * Unsigned (e.g., :attr:`np.uint8`), a new signed array preserving the
      contents of this array is returned.

    In the latter case, this array's contents are preserved across this data
    type conversion by doubling the precision of this array. Since signed and
    unsigned data types of the same bit width cover different integer ranges,
    creating a signed array of same bit width does *not* suffice to preserve
    this array's contents in the general case.

    For example, consider the following data types:

    * :attr:`np.uint16`, covering the integer range ``[0, 65535]``.
    * :attr:`np.int16`, covering the integer range ``[−32768, 32767]``.
    * :attr:`np.int32`, covering the integer range
      ``[−2147483648, 2147483647]``.

    If the data type of the passed array is :attr:`np.uint16`, this implies
    that only a data type of :attr:`np.int32` or wider (e.g., :attr:`np.int64`)
    suffices to preserve the contents of this array. Coercing these contents
    into an array of data type :attr:`np.int16` would destroy (e.g, truncate,
    wrap) *all* array elements in the integer range ``[32768, 65535]``.

    Parameters
    ----------
    array: NumpyArrayType
        Possibly unsigned Numpy array to be converted into a signed array.

    Returns
    ----------
    NumpyArrayType
        Signed Numpy array converted from this possibly unsigned array.
    '''

    # Next largest signed data type guaranteed to preserve array contents if
    # this array is unsigned *OR* "None" otherwise.
    dtype_signed = DTYPE_UNSIGNED_TO_SIGNED.get(array.dtype, None)

    # Return either:
    #
    # * If this array is already signed, this array unmodified.
    # * If this array is unsigned, this array cast to this signed data type.
    return array if dtype_signed is None else dtype_signed(array)
