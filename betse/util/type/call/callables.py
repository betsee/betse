#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level **callable** (e.g., function, lambda, method, property) facilities.
'''

# ....................{ IMPORTS                            }....................
from betse.util.io.log import logs
from betse.util.type.types import (
    type_check,
    CallableTypes,
    BuiltinFunctionType,
    BuiltinMethodType,
    FunctionType,
    FunctionTypes,
    GeneratorType,
    MethodType,
    MethodTypes
)
from functools import wraps

if False: wraps  # silence contemptible IDE warnings

# ....................{ CONSTANTS                          }....................
TYPE_TO_NAME = {
    BuiltinFunctionType: 'builtin',
    BuiltinMethodType:   'method',
    FunctionType:        'function',
    GeneratorType:       'generator',
    MethodType:          'method',
}
'''
Dictionary mapping from the class of each possible callable to the lowercase
human-readable type of that class, ignoring low-level implementation details
(e.g., the distinction between built-in and user-defined methods).
'''

# ....................{ DECORATORS                         }....................
@type_check
def deprecated(func: CallableTypes) -> CallableTypes:
    '''
    Decorate the passed **callable** (e.g., function, method) to log a non-fatal
    warning on the first and _only_ first call of this callable.
    '''

    # Avoid circular import dependencies.
    from betse.util.type import strs

    # True only if this callable has had a deprecation warning logged.
    func.__is_deprecation_logged = False

    # Callable returned by this decorator below.
    @wraps(func)
    def _deprecated_inner(*args, **kwargs) -> object:
        # If this callable has *NOT* had a deprecation warning logged, do so.
        if not func.__is_deprecation_logged:
            # Prevent this warning from being logged more than once.
            func.__is_deprecation_logged = True

            # Capitalized human-readable string describing this callable.
            func_name = strs.uppercase_first_char(to_str(func))

            # Log this warning.
            logs.log_warning('%s deprecated.', func_name)

        # Call this callable.
        return func(*args, **kwargs)

    # Return this callable decorated to log this deprecating warning.
    return _deprecated_inner

# ....................{ TESTERS                            }....................
@type_check
def is_function(func: CallableTypes) -> bool:
    '''
    ``True`` only if the passed callable is a function.
    '''

    return isinstance(func, FunctionTypes)


@type_check
def is_method(func: CallableTypes) -> bool:
    '''
    ``True`` only if the passed callable is a method.
    '''

    return isinstance(func, MethodTypes)

# ....................{ CONVERTERS                         }....................
@type_check
def to_str(func: CallableTypes) -> str:
    '''
    Human-readable string describing the type and name of the passed callable
    (e.g., ``method Pala()``).

    For generality, this string's:

    * First character is guaranteed to be lowercase.
    * Last character is guaranteed to be the ``)`` delimiter.

    Parameters
    ----------
    func: CallableTypes
        Callable to describe.

    Returns
    ----------
    str
        Human-readable string describing this callable.
    '''

    # Type of this callable, defaulting to merely "callable" if unrecognized.
    func_type = TYPE_TO_NAME.get(type(func), 'callable')
    # print('func: {}, type: {}'.format(func.__name__, type(func)))

    # Return a human-readable string describing this callable.
    return '{} {}()'.format(func_type, func.__name__)
