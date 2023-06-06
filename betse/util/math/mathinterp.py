#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2023 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Project-wide **interpolators** (i.e., low-level callables interpolating two- and
three-dimensional meshes defined from one set of points onto another set of
points).
'''

# ....................{ IMPORTS                            }....................
from beartype import beartype
from beartype.typing import Callable
from functools import partial
from numpy import ndarray
from scipy.interpolate import (
    bisplev,
    bisplrep,
)

# ....................{ INTERPOLATORS                      }....................
@beartype
def interp2d_linear(x: ndarray, y: ndarray, z: ndarray) -> (
    Callable[[ndarray, ndarray], ndarray]):
    '''
    Create and return **two-dimensional interpolator** (i.e., low-level callable
    interpolating a source mesh defined on the passed set of X and Y coordinates
    and Z values (e.g., colours) onto a target mesh defined a new of X and Y
    coordinates passed to that callable, which then returns the new Z values
    produced by that interpolation).

    This function effectively implements the obsolete
    :class:`scipy.interpolate.interp2d` class (passed the parameters
    ``kind='linear'`` and ``fill_value=0.0``) in terms of non-obsolete
    functionality defined by the :mod:`scipy.interpolate` subpackage, improving
    backward compatibility with existing code expecting that class to exist.

    Parameters
    ----------
    x: ndarray
        Source X coordinates to interpolate from.
    y: ndarray
        Source Y coordinates to interpolate from.
    z: ndarray
        Source Z values to interpolate from.

    Returns
    ----------
    Callable[[ndarray, ndarray], ndarray]
        Two-dimensional interpolator as described above. The signature of this
        callable should be assumed to be:

        .. code-block:: python

           def interpolator(x_new: ndarray, y_new: ndarray) -> ndarray: ...
    '''

    #FIXME: Comment this out to test our new interp2d() replacement, which is
    #currently broken for unknown reasons. For example, this test fails:
    #    ./pytest -k test_cli_sim_compat
    from scipy.interpolate import interp2d
    return interp2d(x, y, z, kind='linear', fill_value=0.0)

    # Linear B-spline representation of this two-dimensional surface.
    #
    # Note that the "kx=1" and "ky=1" parameters effectively produce a linear
    # interpolation.
    tck = bisplrep(x, y, z, kx=1, ky=1)

    # Return the desired interpolator for this linear B-spline representation.
    #
    # Note that usage of partial() here is required to avoid infinite recursion
    # from the standard "pickle" and third-party "dill" modules. In particular,
    # attempting to replace this partial() approach with a simpler lambda
    # function or closure induces infinite recursion for unknown reasons.
    return partial(_interpolator2d_linear, _tck=tck)

# ....................{ PRIVATE ~ interpolators            }....................
#FIXME: This function is intentionally *NOT* type-checked via @beartype, as both
#the parameters and return values appear to unexpectedly expect NumPy scalars as
#well as arrays. For now, let's just permissively run with duck typing. *sigh*
def _interpolator2d_linear(x_new, y_new, _tck) -> ndarray:
    '''
    Z values (e.g., colours) produced by interpolating from the passed source
    linear B-spline representation describing the source mesh defined on the
    set of X and Y coordinates and Z values (e.g., colours) previously passed to
    the :func:`.interp2d_linear` function onto a target mesh defined on the
    passed X and Y coordinates.

    Parameters
    ----------
    x: ndarray
        Target X coordinates to interpolate onto.
    y: ndarray
        Target Y coordinates to interpolate onto.
    _tck : ndarray
        List ``[tx, ty, c, kx, ky]`` containing the knots ``(tx, ty)`` and
        coefficients ``(c)`` of the bivariate B-spline representation of the
        surface along with the degree of the spline, typically obtained by a
        prior call to the :func:`.bisplrep` function.

    Returns
    ----------
    ndarray
        Z values (e.g., colours) produced by this interpolation.
    '''

    # Perform this interpolation. Note that the transpose (i.e., ".T") mimics
    # the behaviour of the original scipy.interpolate.interp2d() function.
    return bisplev(x_new, y_new, _tck).T
