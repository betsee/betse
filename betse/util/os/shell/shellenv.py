#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **environment** (i.e., set of all external shell variables exported
to the the active Python interpreter) facilities.
'''

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: For safety, all attempts to get, set, or unset environment variables
# should act on the low-level global "os.environ" dictionary directly rather
# than calling the high-level getenv(), setenv(), or unsetenv() functions. To
# quote official Python documentation:
#
#     Calling putenv() directly does not change os.environ, so it's better to
#     modify os.environ.
#
# See http://docs.python.org/library/os.html#os.environ.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# ....................{ IMPORTS                           }....................
from os import environ
from betse.exceptions import BetseShellEnvException
from betse.util.type.types import type_check, MappingType, StrOrNoneTypes

# ....................{ EXCEPTIONS                        }....................
@type_check
def die_unless_var(*names: str) -> None:
    '''
    Raise an exception unless the shell environment for the current process
    defines all environment variables with the passed names.

    Parameters
    ----------
    names : tuple[str]
        Tuple of the names of all environment variables to test for.

    Raises
    ----------
    BetseShellEnvException
        If any such variable is undefined.
    '''

    # If any such variable is undefined...
    if not is_var(*names):
        # For the name of each such variable...
        for name in names:
            # If this variable is undefined, raise a human-readable exception.
            # While directly accessing this variable via "environ[name]" would
            # implicitly raise a "KeyError" exception and hence arguably
            # suffice, that exception's message is *NOT* human-readable.
            if name not in environ:
                raise BetseShellEnvException(
                    'Environment variable "{}" undefined.'.format(name))

# ....................{ TESTERS                           }....................
@type_check
def is_var(*names: str) -> bool:
    '''
    ``True`` only if the shell environment for the current process defines all
    environment variables with the passed names.

    Parameters
    ----------
    names : tuple[str]
        Tuple of the names of all environment variables to test for.

    Returns
    ----------
    bool
        ``True`` only if all such variables are defined.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.iterable.mapping import maptest

    # Return true only if the names of these environment variables are all keys
    # of the corresponding dictionary global.
    return maptest.has_keys(mapping=environ, keys=names)

# ....................{ GETTERS                           }....................
def get_env() -> MappingType:
    '''
    Dictionary mapping the name of each environment variable to the string
    value of this variable.

    For safety, this dictionary is guaranteed to be a deep copy of the current
    environment and hence may be safely:

    * Modified elsewhere *without* modifying the true environment.
    * Isolate this environment to subprocesses, preventing concurrent changes
      in the environment of this process from affecting these subprocesses.
    '''

    return environ.copy()

# ....................{ GETTERS ~ var                     }....................
@type_check
def get_var(name: str) -> str:
    '''
    String value of the environment variable with the passed name if defined
    *or* raise an exception otherwise (i.e., if this variable is undefined).

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
    BetseShellEnvException
        If this variable is undefined.
    '''

    # If this variable is undefined, raise an exception.
    die_unless_var(name)

    # Else, return this variable's value.
    return environ[name]


@type_check
def get_var_or_default(name: str, default: str) -> str:
    '''
    String value of the environment variable with the passed name if defined
    *or* the passed default string value otherwise.

    Parameters
    ----------
    name : str
        Name of this variable.
    default : str
        Default string value of this variable.

    Returns
    ----------
    str
        Value of this variable if defined *or* the passed default otherwise.
    '''

    return environ.get(name, default)


@type_check
def get_var_or_none(name: str) -> StrOrNoneTypes:
    '''
    String value of the environment variable with the passed name if defined
    *or* ``None`` otherwise.

    Parameters
    ----------
    name : str
        Name of this variable.

    Returns
    ----------
    StrOrNoneTypes
        Value of this variable if defined *or* ``None`` otherwise.
    '''

    return environ.get(name, None)

# ....................{ SETTERS                           }....................
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

    environ[name] = value

# ....................{ UNSETTERS                         }....................
@type_check
def unset_var(name: str) -> None:
    '''
    Unset the environment variable with the passed name (i.e., remove this
    variable from the shell environment for the current process) if defined
    *or* raise an exception otherwise (i.e., if this variable is undefined).

    Parameters
    ----------
    name : str
        Name of this variable.

    Raises
    ----------
    BetseShellEnvException
        If this variable is undefined.
    '''

    # If this variable is undefined, raise an exception.
    die_unless_var(name)

    # Unset this variable.
    unset_var_if_set(name)


@type_check
def unset_var_if_set(name: str) -> None:
    '''
    Unset the environment variable with the passed name (i.e., remove this
    variable from the shell environment for the current process) if defined
    *or* noop otherwise (i.e., if this variable is undefined).

    Parameters
    ----------
    name : str
        Name of this variable.
    '''

    # Reduce this variable to the empty string *BEFORE* unsetting this
    # variable, thus:
    #
    # * Preventing the subsequent attempt to unset this variable from raising
    #   exceptions if undefined.
    # * Handling edge-case platforms whose kernels fail to support the
    #   os.unsetenv() operation internally invoked by deleting keys from the
    #   "os.environ" dictionary (e.g., AIX). Under such platforms, deleting
    #   environment variables in Python fails to delete these variables from
    #   the environment outside of Python. Although reducing this variable to
    #   the empty string does *NOT* delete this variable, no alternatives
    #   exist.
    environ[name] = ""

    # Unset this variable.
    del environ[name]

# ....................{ CONVERTERS                        }....................
def to_str() -> str:
    '''
    Human-readable string of all environment variables.

    This string is intended to be substantially more human-readable than that
    returned by the :meth:`os._Environ.__str__` special method implicitly
    invoked by the :meth:`str.__init__` method, which compresses the entirety
    of this environment onto a single line.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.iterable.mapping import mappings

    # Return the human-readable string produced by joining on newline a
    # generator comprehension yielding the colon-delimited name and value of
    # all environment variables sorted in lexicographic order.
    #
    # Note that the "environ" object is of non-standard type "os._Environ",
    # which the pprint.pformat() fails to recognize and hence format as a
    # "dict"-compatible mapping. Ergo, passing "environ" directly to the
    # mappings.to_str_flat() function would yield a non-human-readable string.
    # While this could be ameliorated by converting "environ" to a "dict" first
    # (e.g., "mappings.to_str_flat(dict(environ))"), doing so still produces
    # less human-readable output than the current approach.
    return mappings.to_str_flat(environ)
