#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level Numpy array and matrix facilities.
'''

# ....................{ IMPORTS                            }....................
import numpy as np
from betse.util.type.types import type_check, SequenceTypes
from numpy import ndarray

# ....................{ GETTERS                            }....................
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
