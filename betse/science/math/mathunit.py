#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **units** (i.e., standard for quantitative measurement of observable
quantities) facilities.
'''

# ....................{ IMPORTS                           }....................
from betse.lib.numpy import nparray
from betse.util.type import types
# from betse.util.type.numeric import ints
from betse.util.type.types import (
    type_check, NumericSimpleTypes, NumericOrIterableTypes,)

# ....................{ CONSTANTS ~ inverse               }....................
INVERSE_CENTI = 1e2
'''
Inverse of the centi- unit prefix (i.e., ``10**−2``), commonly used as a
multiplicative factor for upscaling quantities from centi-prefixed units to
unprefixed units.
'''


INVERSE_MILLI = 1e3
'''
Inverse of the milli- unit prefix (i.e., ``10**−3``), commonly used as a
multiplicative factor for upscaling quantities from milli-prefixed units to
unprefixed units.
'''


INVERSE_MICRO = 1e6
'''
Inverse of the micro- unit prefix (i.e., ``10**−6``), commonly used as a
multiplicative factor for upscaling quantities from micro-prefixed units to
unprefixed units.

Specifically, this is equivalent to the multiplicative factor for upscaling
cell cluster coordinates from micro-prefixed units to unprefixed units (e.g.,
to improve the readability of these coordinates in user-friendly exports).
'''

# ....................{ UPSCALERS ~ cell : data           }....................
@type_check
def upscale_units_centi(
    data: NumericOrIterableTypes) -> NumericOrIterableTypes:
    '''
    Upscale the contents of the passed number or iterable of numbers whose
    units are assumed to be centi-prefixed (i.e., factors of ``10**-2``).

    This function does *not* modify the passed object. If this object is:

    * A scalar number (e.g., :class:`int`, :class:`float`), this function
      returns a new scalar number of the same type multiplied by a positive
      constant.
    * A non-Numpy iterable (e.g., :class:`list`), this function returns a new
      Numpy array first converted from this iterable and then multiplied by a
      positive constant.
    * A Numpy array, this function returns a new Numpy array first copied from
      the passed array and then multiplied by a positive constant.

    Parameters
    ----------
    data : NumericOrIterableTypes
        Either a number of iterable of numbers to upscale.

    Returns
    ----------
    object
        Upscaled object as described above.
    '''

    return _upscale_data_in_units(data=data, factor=INVERSE_CENTI)


@type_check
def upscale_units_milli(
    data: NumericOrIterableTypes) -> NumericOrIterableTypes:
    '''
    Upscale the contents of the passed number or iterable of numbers whose
    units are assumed to be milli-prefixed (i.e., factors of ``10**-3``).

    See Also
    ----------
    :func:`upscale_units_centi`
        Further details.
    '''

    return _upscale_data_in_units(data=data, factor=INVERSE_MILLI)


@type_check
def upscale_units_micro(
    data: NumericOrIterableTypes) -> NumericOrIterableTypes:
    '''
    Upscale the contents of the passed number or iterable of numbers whose
    units are assumed to be micro-prefixed (i.e., factors of ``10**-6``).

    See Also
    ----------
    :func:`upscale_units_centi`
        Further details.
    '''

    return _upscale_data_in_units(data=data, factor=INVERSE_MICRO)

# ....................{ UPSCALERS ~ cell : coordinates    }....................
@type_check
def upscale_coordinates(
    coordinates: NumericOrIterableTypes) -> NumericOrIterableTypes:
    '''
    Upscale the contents of the passed number or sequence of numbers whose
    units are assumed to be denominated in micrometers (i.e., ``10**-6``),
    typically for converting spatial data into X and Y coordinates suitable for
    plot and animation display.

    Parameters
    ----------
    coordinates : NumericOrIterableTypes
        Number or sequence to be upscaled.

    Returns
    ----------
    NumericOrIterableTypes
        Either:

        * If this data is numeric, this number upscaled by this multiplier.
        * If this data is sequential, this sequence converted into a Numpy
          array whose elements are upscaled by this multiplier.

    See Also
    ----------
    :func:`upscale_units_milli`
        Further details.
    '''

    return upscale_units_micro(coordinates)


@type_check
def upscale_coordinates_tuple(
    *cells_coordinates: NumericOrIterableTypes) -> tuple:
    '''
    Upscale the contents of each passed number or sequence of numbers whose
    units are assumed to be denominated in micrometers (i.e., ``10**-6``),
    returning an n-tuple of each upscaled number or sequence of numbers.

    Parameters
    ----------
    cells_coordinates : tuple[NumericOrIterableTypes]
        Tuple of all numbers and sequences to be upscaled.

    Returns
    ----------
    tuple
        N-tuple for ``N`` the number of passed positional arguments such that
        each element is the corresponding positional argument upscaled by the
        :func:`upscale_coordinates` function.
    '''

    return tuple(
        upscale_coordinates(cell_coordinates)
        for cell_coordinates in cells_coordinates
    )

# ....................{ UPSCALERS ~ private               }....................
@type_check
def _upscale_data_in_units(
    data: NumericOrIterableTypes, factor: NumericSimpleTypes) -> (
    NumericOrIterableTypes):
    '''
    Upscale the contents of the passed object by the passed multiplier,
    typically under the assumption that these contents are denominated in the
    reciprocal units of this multiplier (i.e., ``10**6`` for micrometers).

    Parameters
    ----------
    data : NumericOrIterableTypes
        Number or sequence to be upscaled.
    factor : NumericSimpleTypes
        Reciprocal of the units this object is denominated in.

    Returns
    ----------
    NumericOrIterableTypes
        Either:

        * If this data is numeric, this number upscaled by this multiplier.
        * If this data is sequential, this sequence converted into a Numpy
          array whose elements are upscaled by this multiplier.

    See Also
    ----------
    :func:`upscale_units_milli`
        Further details.
    '''

    # If the passed object is numeric, return this number upscaled.
    if types.is_numeric(data):
        return factor * data
    # Else, this object is a sequence. Return this sequence converted into a
    # Numpy array and then upscaled.
    else:
        return factor * nparray.from_iterable(data)
