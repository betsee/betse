#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level **environment** (i.e., set of all external shell variables exported to
the the active Python interpreter) facilities.
'''

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: For safety, all attempts to get, set, or unset environment variables
# should act on the low-level global "os.environ" dictionary directly rather
# than calling the high-level getenv(), setenv(), or unsetenv() functions. To
# quote official Python documentation:
#
#     Calling putenv() directly does not change os.environ, so it's better to
#     modify os.environ.
#
# See http://docs.python.org/library/os.html#os.environ.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# ....................{ IMPORTS                            }....................
import os
from betse.util.type.types import type_check, MappingType

# ....................{ TESTERS                            }....................
@type_check
def is_var(*names: str) -> bool:
    '''
    `True` only if all environment variables with the passed names are defined.

    Parameters
    ----------
    names : tuple[str]
        Tuple of the names of all environment variables to test for.

    Returns
    ----------
    bool
        `True` only if all such environment variables are defined.
    '''

    # Avoid circular import dependencies.
    from betse.util.type import mappings

    # Return True only if the names of these environment variables are all keys
    # of the corresponding dictionary global.
    return mappings.is_keys(os.environ, *names)

# ....................{ GETTERS                            }....................
def get_env() -> MappingType:
    '''
    Dictionary mapping the name of each environment variable to the string value
    of this variable.

    For safety, this dictionary is guaranteed to be a deep copy of the current
    environment and hence may be safely:

    * Modified elsewhere _without_ modifying the true environment.
    * Isolate this environment to subprocesses, preventing concurrent changes in
      the environment of this process from affecting these subprocesses.
    '''

    return os.environ.copy()

# ....................{ GETTERS ~ var                      }....................
@type_check
def get_var(name: str) -> str:
    '''
    Value of the environment variable with the passed name.

    If this variable is _not_ defined, an exception is raised.

    Parameters
    ----------
    name : str
        Name of this variable.

    Returns
    ----------
    str
        Value of this variable.

    Raises
    ----------
    KeyError
        If this variable is undefined.
    '''

    return os.environ[name]


@type_check
def get_var_or_default(name: str, default: str) -> str:
    '''
    Value of the environment variable with the passed name, defaulting to the
    passed value if this variable is undefined.

    Parameters
    ----------
    name : str
        Name of this variable.
    default : str
        Default value of this variable.

    Returns
    ----------
    str
        Value of this variable if defined _or_ the passed default otherwise.
    '''

    return os.environ.get(name, default)

# ....................{ SETTERS                            }....................
@type_check
def set_var(name: str, value: str) -> None:
    '''
    Set the environment variable with the passed name to the passed value.

    Parameters
    ----------
    name : str
        Name of this variable.
    value : str
        New value of this variable.
    '''

    os.environ[name] = value

# ....................{ REMOVERS                           }....................
@type_check
def unset_var(name: str) -> None:
    '''
    Unset (i.e., remove) the environment variable with the passed name.

    If this variable is _not_ currently set, an exception is raised.

    Parameters
    ----------
    name : str
        Name of this variable.

    Raises
    ----------
    KeyError
        If this variable is unset.
    '''

    # If this variable is unset, raise an exception.
    os.environ[name]

    # Unset this variable.
    unset_var_if_set(name)


@type_check
def unset_var_if_set(name: str) -> None:
    '''
    Unset (i.e., remove) the environment variable with the passed name if set
    _or_ noop otherwise.

    Parameters
    ----------
    name : str
        Name of this variable.
    '''

    # Reduce this variable to the empty string *BEFORE* unsetting this variable,
    # thus:
    #
    # * Preventing the subsequent attempt to unset this variable from raising
    #   exceptions if undefined.
    # * Handling edge-case platforms whose kernels fail to support the
    #   os.unsetenv() operation internally invoked by deleting keys from the
    #   "os.environ" dictionary (e.g., AIX). Under such platforms, deleting
    #   environment variables in Python fails to delete these variables from the
    #   environment outside of Python. Although reducing this variable to the
    #   empty string does *NOT* delete this variable, no alternatives exist.
    os.environ[name] = ""

    # Unset this variable.
    del os.environ[name]
