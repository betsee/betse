#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Utility functions of general-purpose relevance to all plots and animations.
'''

# ....................{ IMPORTS                            }....................
from betse.lib.numpy import arrays
from betse.util.type import types
from betse.util.type.types import (
    type_check, NumericOrSequenceTypes, NumericTypes)

# ....................{ UPSCALERS                          }....................
#FIXME: Use everywhere in the "betse.science.visual.anim.animpipe" submodule.
@type_check
def upscale_cell_data(
    cell_data: NumericOrSequenceTypes) -> NumericOrSequenceTypes:
    '''
    Upscale the contents of the passed object whose units are assumed to be
    denominated in milli- (i.e., `10**-3`).

    This function does _not_ modify the passed object. If this object is:

    * A scalar number (e.g., of type :class:`int` or :class:`float`), this
      function returns a new scalar number of the same type multiplied by a
      positive constant.
    * A NumPy array, this function returns a new NumPy array equal to the
      passed array multiplied by a positive constant.
    * A Python sequence (e.g., :class:`list`), this function returns a new NumPy
      array equal to the passed sequence converted into a new NumPy array and
      then multiplied by a positive constant.

    Parameters
    ----------
    cell_data : NumericOrSequenceTypes
        Number or sequence to be upscaled.

    Returns
    ----------
    object
        Upscaled object as described above.
    '''

    return _upscale_data_in_units(cell_data, 1000)


@type_check
def upscale_cell_coordinates(
    cell_coordinates: NumericOrSequenceTypes) -> NumericOrSequenceTypes:
    '''
    Upscale the contents of the passed object whose units are assumed to be
    denominated in micrometers (i.e., `10**-6`), typically to convert spatial
    data into X and Y coordinates for plot and animation artists.

    Parameters
    ----------
    cell_coordinates : NumericOrSequenceTypes
        Number or sequence to be upscaled.

    Returns
    ----------
    object
        Upscaled object as described above.

    See Also
    ----------
    :func:`upscale_cell_data`
        Further details.
    '''

    return _upscale_data_in_units(cell_coordinates, 1000000)


@type_check
def _upscale_data_in_units(
    data: NumericOrSequenceTypes, multiplier: NumericTypes) -> (
    NumericOrSequenceTypes):
    '''
    Upscale the contents of the passed object by the passed multiplier under
    the assumption that the units of these contents are the reciprocal of this
    multiplier (i.e., `10**6` for micrometers).

    Parameters
    ----------
    data : NumericOrSequenceTypes
        Number or sequence to be upscaled.
    multiplier : NumericTypes
        Reciprocal of the units this object is denominated in.

    Returns
    ----------
    object
        Upscaled object as described above.

    See Also
    ----------
    :func:`upscale_cell_data`
        Further details.
    '''

    # If the passed object is numeric, return this number upscaled.
    if types.is_numeric(data):
        return data * multiplier
    # Else, this object is a sequence. Return this sequence converted into a
    # Numpy array and then upscaled.
    else:
        return arrays.from_sequence(data) * multiplier
