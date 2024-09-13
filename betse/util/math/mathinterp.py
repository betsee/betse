#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2025 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Project-wide **mathematical interpolators** (i.e., low-level callables
interpolating two- and three-dimensional meshes defined from one set of points
onto another set of points).
'''

# ....................{ IMPORTS                            }....................
from beartype.typing import (
    Callable,
)
from betse.exceptions import BetseMathInterpolationException
from betse.util.type.typehints import NDArrayNdim1, NDArrayNdim2
# from functools import partial
from numpy import ndarray
from scipy.interpolate import (
    RectBivariateSpline,
    # bisplev,
    # bisplrep,
)

# ....................{ INTERPOLATORS                      }....................
def interp2d_linear(
    x: NDArrayNdim1,
    y: NDArrayNdim1,
    z: NDArrayNdim2,
) -> (
    Callable[[ndarray, ndarray], ndarray]):
    '''
    Create and return a new **two-dimensional interpolator** (i.e., low-level
    callable interpolating a source mesh defined on the passed set of X and Y
    coordinates comprising a so-called "regular grid" and Z values (e.g.,
    colours) onto a target mesh defined on a new set of X and Y coordinates
    passed to that callable, which then returns the new Z values produced by
    that interpolation).

    This interpolator implements the "regular grid" code path of the obsolete
    :class:`scipy.interpolate.interp2d` class (passed the parameters
    ``kind='linear'`` and ``fill_value=0.0``) in terms of non-obsolete
    functionality defined by the :mod:`scipy.interpolate` subpackage, improving
    backward compatibility with existing code expecting that class to exist.

    Caveats
    -------
    This interpolator requires that the passed arrays satisfy these constraints:

    .. code-block:: python

       # The X and Y coordinates are *ALL* 1-dimensional: i.e.,
       x.ndim == y.ndim == 1

       # The Z values are 2-dimensional: i.e.,
       z.ndim == 2

       # The total number of X and Y coordinates differ: i.e.,
       x.size != y.size

       # The total number of Z values is the product of:
       # * The total number of X coordinates and...
       # * The total number of Y coordinates.
       z.size == x.size * y.size

    If any one of these constraints is *not* the case, an exception is raised.

    Parameters
    ----------
    x: NDArrayNdim1
        Source 1-dimensional X coordinates to interpolate from.
    y: NDArrayNdim1
        Source 1-dimensional Y coordinates to interpolate from.
    z: NDArrayNdim2
        Source 2-dimensional Z values to interpolate from.

    Returns
    -------
    Callable[[ndarray, ndarray], ndarray]
        Two-dimensional interpolator as described above. The signature of this
        callable should be assumed to be:

        .. code-block:: python

           def interpolator(x_new: ndarray, y_new: ndarray) -> ndarray: ...

    Raises
    ------
    BetseMathInterpolationException
        If any of the above constraints are violated.
    '''

    # If the total number of Z values is *NOT* the product of the total number
    # of X and Y coordinates, raise an exception.
    if z.size != x.size * y.size:
        raise BetseMathInterpolationException(
            f'Interpolation Z shape not product of X and Y coordinates shape '
            f'(i.e., "z.size = {z.size}" != '
            f'"x.size = {x.size}" * "y.size = {y.size}").'
        )
    # Else, the total number of Z values is the product of the total number
    # of X and Y coordinates.

    #FIXME: This doesn't appear to hold. Square shapes are fine. *shrug*
    # # If the number of X and Y coordinates are the same, raise an exception.
    # if x.size == y.size:
    #     raise BetseMathInterpolationException(
    #         f'Interpolation X and Y coordinates share same size (i.e., '
    #         f'"x.size = {x.size}" and "y.size = {y.size}" equal).'
    #     )
    # # Else, the number of X and Y coordinates differ.

    # #FIXME: Comment this out to test our new interp2d() replacement, which is
    # #currently broken for unknown reasons. For example, this test fails:
    # #    ./pytest -k test_cli_sim_compat
    # from scipy.interpolate import interp2d
    # return interp2d(x, y, z, kind='linear', fill_value=0.0)

    #FIXME; Non-ideal. Modern SciPy documentation strongly advises usage of:
    #    For legacy code, nearly bug-for-bug compatible replacements are
    #    `RectBivariateSpline` on regular grids, and `bisplrep`/`bisplev` for
    #    scattered 2D data.
    #    
    #    In new code, for regular grids use `RegularGridInterpolator` instead.
    #    For scattered data, prefer `LinearNDInterpolator` or
    #    `CloughTocher2DInterpolator`.
    #
    #    For more details see
    #    https://scipy.github.io/devdocs/tutorial/interpolate/interp_transition_guide.html
    #
    #We're fairly certain we should be leveraging "RectBivariateSpline" rather
    #than "bisplrep" here. The above SciPy documentation as well as Ally's
    #perfect memory suggests a similar finding. Consider refactoring, please.
    #
    #In short, there appear to exist two code paths. The first code path is
    #given by:
    #    # If the total number of X and Y coordinates differ *AND*...
    #    if x.size != y.size and z.size == len(x) * len(y):
    #       then do "RectBivariateSpline" stuff.
    #    elif len(x) == len(y) == len(z):
    #       then do "bisplrep" stuff.
    #    else:
    #       raise Exception(!)

    # Linear bivariate spline approximation over this rectangular mesh.
    #
    # Note that:
    # * The "kx=1" and "ky=1" parameters produce a linear interpolation.
    # * The "z.T" transpose reproduces the behaviour of the obsolete
    #   scipy.interpolate.interp2d() function that this higher-level function
    #   descended from and preserves backward API compatibility with.
    bivariate_spline = RectBivariateSpline(x, y, z.T, kx=1, ky=1)

    # Linear bivariate spline interpolator over this rectangular mesh.
    #
    # Note that:
    # * The ".T" transpose reproduces the behaviour of the obsolete
    #   scipy.interpolate.interp2d() function that this higher-level function
    #   descended from and preserves backward API compatibility with.
    bivariate_spline_interpolator = lambda x_new, y_new: (
        bivariate_spline(x_new, y_new).T)

    # Return this interpolator.
    return bivariate_spline_interpolator

    #FIXME: Currently unused, but preserved for posterity. *shrug*
    # # Linear B-spline representation of this two-dimensional surface.
    # #
    # # Note that the "kx=1" and "ky=1" parameters effectively produce a linear
    # # interpolation.
    # tck = bisplrep(x, y, z, kx=1, ky=1)
    #
    # # Return the desired interpolator for this linear B-spline representation.
    # #
    # # Note that usage of partial() here is required to avoid infinite recursion
    # # from the standard "pickle" and third-party "dill" modules. In particular,
    # # attempting to replace this partial() approach with a simpler lambda
    # # function or closure induces infinite recursion for unknown reasons.
    # return partial(_interpolator2d_linear, _tck=tck)

# ....................{ PRIVATE ~ interpolators            }....................
#FIXME: Currently unused, but preserved for posterity. *shrug*
# #FIXME: This function is intentionally *NOT* type-checked via @beartype, as both
# #the parameters and return values appear to unexpectedly expect NumPy scalars as
# #well as arrays. For now, let's just permissively run with duck typing. *sigh*
# def _interpolator2d_linear_scattered(x_new, y_new, _tck) -> ndarray:
#     '''
#     Z values (e.g., colours) produced by interpolating from the passed source
#     linear B-spline representation describing the source mesh defined on the
#     set of X and Y coordinates and Z values (e.g., colours) previously passed to
#     the :func:`.interp2d_linear` function onto a target mesh defined on the
#     passed X and Y coordinates.
#
#     Parameters
#     ----------
#     x: ndarray
#         Target X coordinates to interpolate onto.
#     y: ndarray
#         Target Y coordinates to interpolate onto.
#     _tck : ndarray
#         List ``[tx, ty, c, kx, ky]`` containing the knots ``(tx, ty)`` and
#         coefficients ``(c)`` of the bivariate B-spline representation of the
#         surface along with the degree of the spline, typically obtained by a
#         prior call to the :func:`.bisplrep` function.
#
#     Returns
#     -------
#     ndarray
#         Z values (e.g., colours) produced by this interpolation.
#     '''
#
#     # Perform this interpolation. Note that the transpose (i.e., ".T") mimics
#     # the behaviour of the original scipy.interpolate.interp2d() function.
#     return bisplev(x_new, y_new, _tck).T
