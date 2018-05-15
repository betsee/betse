#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **decorator-based profiling** (i.e., debugging reports of the space
and/or time complexity of decorated callables) facilities.
'''

# ....................{ IMPORTS                            }....................
import time
from betse.util.io.log import logs
from betse.util.type.types import (
    type_check, CallableTypes, CallableOrNoneTypes, StrOrNoneTypes)
from functools import wraps

# ....................{ DECORATORS                         }....................
# Decorators accepting optional arguments require additional... finesse.
@type_check
def log_time_seconds(
    # Callable to be decorated passed positionally if at all by Python itself.
    # If:
    #
    # * Passed, this function is being called as a decorator.
    # * Unpassed, this function is being called as a decorator factory expected
    #   to produce the decorator closure specific to the set of all
    #   subsequently passed optional parameters.
    func: CallableOrNoneTypes = None,
    # Require *ALL* subsequent optional parameters to be passed as
    # non-ambiguous keyword arguments rather than ambiguous positional
    # arguments. If *NOT* done, there would exist no means of deciding whether
    # this function is being called as a decorator or decorator factory in the
    # edge case of a decorator passed a callable as a first optional parameter.
    *,
    # Optional parameters passed to this function when called as a decorator
    # factory rather than as a decorator.
    noun: StrOrNoneTypes = None,
    verb: str = 'completed'
) -> CallableTypes:
    '''
    Decorate the passed callable (e.g., function, lambda, method) to log the
    **wall clock time** (i.e., cumulative time spent in both kernel- and
    userspace while in the body of this callable) of each call to this
    callable, denominated in possibly fractional seconds.

    Parameters
    ----------
    func : CallableOrNoneTypes
        Callable to be decorated if this function is being called as a
        decorator *or* ``None`` if this function is being called as a
        decorator factory.
    noun : StrOrNoneTypes
        Human-readable string describing this callable to be logged with the
        time spent in each call of this callable. If the first character of
        this string is *not* already uppercased, this function implicitly
        uppercases this character for readability. Defaults to ``None``, in
        which case a human-readable string describing this callable is
        defaulted to (e.g., ``Method skromt_og_kolabrenning()``).
    verb : str
        Human-readable string describing the action performed by this callable
        to be logged with the time spent in each call of this callable.
        Defaults to a general-purpose verb.

    See Also
    ----------
    https://stackoverflow.com/a/24617244/2809027
        StackOverflow answer inspiring this clever (albeit obtuse, admittedly)
        decorator design.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.call import callables
    from betse.util.type.text import strs

    # Decorator factory creating a decorator specific to the passed parameters.
    def _log_time_seconds_decorator(func: CallableTypes) -> CallableTypes:
        '''
        Create and return a closure decorating the passed callable to log the
        **wall clock time** (i.e., cumulative time spent in both kernel- and
        userspace while in the body of this callable) of each call to this
        callable, denominated in possibly fractional seconds.
        '''

        # Permit these enclosed variables from the outer scope to be modified.
        nonlocal noun

        # If unpassed, default this noun to the name of this callable.
        if noun is None:
            noun = callables.to_str(func)

        # Uppercase the first character of this string for readability.
        noun = strs.uppercase_char_first(noun)

        # Closure decorating this callable to be returned.
        @wraps(func)
        def _log_time_seconds_decorated(*args, **kwargs) -> object:

            # Current time in fractional seconds.
            start_time = time.time()

            # Call this function, passed all passed parameters and preserving the
            # return value as is.
            return_value = func(*args, **kwargs)

            # Cumulative time in fractional seconds spent in this call.
            end_time = time.time() - start_time

            # Log this time, rounded to two decimal places for readability.
            logs.log_info('%s %s in %.2f seconds.', noun, verb, end_time)

            # Return this return value.
            return return_value

        # Return this decorated callable.
        return _log_time_seconds_decorated

    # Return...
    return (
        # If this function is called as a decorator factory, this factory.
        _log_time_seconds_decorator if func is None else
        # Else, this function is called as a decorator. In this case, return a
        # decorator specific to the passed parameters.
        _log_time_seconds_decorator(func)
    )
