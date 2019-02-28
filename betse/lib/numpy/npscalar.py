#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **Numpy-specific scalar** (i.e., primitive non-standard data types
specific to the Numpy API) facilities.
'''

# ....................{ IMPORTS                           }....................
# import numpy as np
# from betse.util.io.log import logs
from betse.util.type.types import type_check, NumpyScalarType, ScalarTypes

# ....................{ CONVERTERS                        }....................
@type_check
def to_python(scalar: NumpyScalarType) -> ScalarTypes:
    '''
    Standard Numpy-agnostic scalar (e.g., :class:`bool`) losslessly coerced
    from the passed non-standard Numpy-specific scalar (e.g.,
    :class:`np.bool_`).

    This function converts the passed Numpy scalar to the corresponding Python
    scalar such that the latter is guaranteed to be a **lossless copy** (i.e.,
    copy with *no* appreciable loss in numerical precision) of the former.

    Parameters
    ----------
    scalar : NumpyScalarType
        Numpy-specific scalar to be converted into the corresponding
        Numpy-agnostic scalar.

    Returns
    ----------
    ScalarTypes
        Numpy-agnostic scalar converted from this Numpy-specific scalar.

    See Also
    ----------
    https://stackoverflow.com/a/11389998/2809027
        StackOverflow answer strongly inspiring this implementation.

    Examples
    ----------
        >>> from betse.lib.numpy import npscalar
        >>> import numpy as np
        >>> python_bool = True
        >>> numpy_bool = np.bool_(python_bool)
        >>> numpy_bool
        True
        >>> numpy_bool is True
        False
        >>> numpy_bool == True
        False
        >>> python_numpy_bool = npscalar.to_python(numpy_bool)
        >>> python_numpy_bool
        True
        >>> python_numpy_bool is True
        True
        >>> python_numpy_bool == True
        True
    '''

    # This is insanity. This is Numpy. For further details, see the official
    # documentation for the numpy.ndarray.item() method at:
    #     https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.item.html
    return scalar.item()
