#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Yet Another Markup Language (YAML) **representer** (i.e., callable serializing
all objects of the same type into well-formatted YAML) functionality.
'''

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To permit YAML implementations to be conditionally imported at
# application startup, no implementations (e.g., the top-level "yaml" package
# corresponding to PyYAML) are importable here.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# from betse.util.type.types import type_check  #, MappingOrSequenceTypes

# ....................{ GLOBALS ~ private                  }....................
_TAG_NUMPY_NDARRAY = '!ndarray'
'''
YAML-formatted tag identifying each YAML sequence converted from a Numpy array.
'''

# ....................{ ADDERS                             }....................
def add_representers(representer: object) -> None:
    '''
    Add all **custom representers** (i.e., callables serializing all objects of
    the same type into well-formatted YAML) required by this application to the
    passed :class:`yaml.representer.Representer`-like object.

    Motivation
    ----------
    This function ensures that dumping arbitrary objects to YAML with the passed
    dumper will stringify all:

    * Numpy arrays as YAML sequences.
    * Numpy booleans as YAML booleans.
    * Numpy complex numbers as YAML complex numbers.
    * Numpy floats as YAML floats.
    * Numpy integers as YAML integers.

    Since most YAML implementations lack explicit support for Numpy types,
    attempting to saving Numpy objects to YAML results in non-human-readable
    output defeating the whole purpose of YAML. This output resembles:

        >>> print(yaml.dump({'a':numpy.array([3.1415,])}))
        /usr/lib64/python3.4/site-packages/yaml/representer.py:135:
        FutureWarning: comparison to `None` will result in an elementwise
        object comparison in the future.
          if data in [None, ()]:
        a: !!python/object/apply:numpy.core.multiarray._reconstruct
          args:
          - !!python/name:numpy.ndarray ''
          - !!python/tuple [0]
          - !!binary |
            Yg==
          state: !!python/tuple
          - 1
          - !!python/tuple [1]
          - !!python/object/apply:numpy.dtype
            args: [f8, 0, 1]
            state: !!python/tuple [3, <, null, null, null, -1, -1, 0]
          - false
          - !!binary |
            bxKDwMohCUA=

    Parameters
    ----------
    representer: object
        :class:`yaml.representer.Representer`-like object converting arbitrary
        Python objects to YAML-formatted strings.  This object *must* define an
        ``add_representer`` callable accepting the type to be represented and a
        callable converting objects of that type to YAML strings. This object
        may otherwise be of any type (e.g., module).
    '''

    # Defer third-party imports.
    import numpy

    # Represent all core Numpy types.
    representer.add_representer(numpy.ndarray,  _represent_numpy_ndarray)
    representer.add_representer(numpy.bool_,    _represent_numpy_bool)
    representer.add_representer(numpy.complex_, _represent_numpy_complex)
    representer.add_representer(numpy.float_,   _represent_numpy_float)
    representer.add_representer(numpy.int_,     _represent_numpy_int)

# ....................{ REPRESENTERS                       }....................
def _represent_numpy_ndarray(dumper, ndarray: 'numpy.ndarray') -> str:
    '''
    Convert the passed Numpy array into a YAML-formatted string.

    Specifically, this function returns a stringified YAML mapping whose:

    * Keys are the names of parameters accepted by the `numpy.array()` function.
    * Values are the values of these parameters required to losslessly restore
      the contents of this array on YAML deserialization of this string.

    Parameters
    ----------
    dumper: yaml.Dumper
        Object converting arbitrary Python objects to YAML-formatted strings.
    ndarray: numpy.ndarray
        Numpy array to be converted into a YAML-formatted string.

    Returns
    ----------
    str
        YAML-formatted string representing this Numpy array.
    '''

    # Dictionary mapping from numpy.array() parameter names to values.
    array_params = {
        'object': ndarray.tolist(),
        'dtype': ndarray.dtype.name,
    }

    # Convert this dictionary into a stringified YAML mapping.
    return dumper.represent_mapping(_TAG_NUMPY_NDARRAY, array_params)


def _represent_numpy_bool(dumper, npbool: 'numpy.bool_') -> str:
    '''
    Convert the passed Numpy boolean into a YAML-formatted string.

    Specifically, this function returns a stringified YAML boolean equivalent
    to the passed boolean.

    Parameters
    ----------
    dumper: yaml.Dumper
        Object converting arbitrary Python objects to YAML-formatted strings.
    npbool: numpy.ndarray
        Numpy boolean to be converted into a YAML-formatted string.

    Returns
    ----------
    str
        YAML-formatted string representing this Numpy bool.
    '''

    return dumper.represent_bool(bool(npbool))


def _represent_numpy_complex(dumper, npcomplex: 'numpy.complex_') -> str:
    '''
    Convert the passed Numpy complex number into a YAML-formatted string.

    Specifically, this function returns a stringified YAML complex number
    equivalent to the passed complex number.

    Parameters
    ----------
    dumper: yaml.Dumper
        Object converting arbitrary Python objects to YAML-formatted strings.
    npcomplex: numpy.ndarray
        Numpy complex number to be converted into a YAML-formatted string.

    Returns
    ----------
    str
        YAML-formatted string representing this Numpy complex number.
    '''

    return dumper.represent_complex(complex(npcomplex))


def _represent_numpy_float(dumper, npfloat: 'numpy.float_') -> str:
    '''
    Convert the passed Numpy float into a YAML-formatted string.

    Specifically, this function returns a stringified YAML float equivalent to
    the passed float.

    Parameters
    ----------
    dumper: yaml.Dumper
        Object converting arbitrary Python objects to YAML-formatted strings.
    npfloat: numpy.ndarray
        Numpy float to be converted into a YAML-formatted string.

    Returns
    ----------
    str
        YAML-formatted string representing this Numpy float.
    '''

    return dumper.represent_float(float(npfloat))


def _represent_numpy_int(dumper, npint: 'numpy.int_') -> str:
    '''
    Convert the passed Numpy integer into a YAML-formatted string.

    Specifically, this function returns a stringified YAML integer equivalent to
    the passed integer.

    Parameters
    ----------
    dumper: yaml.Dumper
        Object converting arbitrary Python objects to YAML-formatted strings.
    npint: numpy.ndarray
        Numpy integer to be converted into a YAML-formatted string.

    Returns
    ----------
    str
        YAML-formatted string representing this Numpy integer.
    '''

    return dumper.represent_int(int(npint))
