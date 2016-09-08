#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
High-level **Python profiling** (i.e., measuring various metrics pertaining to
Python code, including time and space performance) facilities.
'''

# ....................{ IMPORTS                            }....................
from betse.util.type.types import (
    type_check, CallableTypes, MappingType, SequenceTypes,)
from functools import partial
from timeit import Timer

# ....................{ TIMERS                             }....................
@type_check
def time_callable(
    callable: CallableTypes,
    args: SequenceTypes = None,
    kwargs: MappingType = None,
    repetitions: int = 7,
    iterations: int = 1000,
) -> float:
    '''
    Time the passed callable with the passed positional and keyword arguments
    (if any) the passed number of repetitions of the passed number of
    iterations, returning the time in seconds consumed by the fastest call to
    this callable.

    This function _only_ returns the minimum timing of all observed timings,
    discarding the time consumed by all calls other than the fastest call to
    this callable. Why? Because the timing of the fastest call is the most
    statistically signifant timing. All error in timing is **positive,**
    typically due to overhead associated with low-level kernel operations (e.g.,
    I/O) and unrelated running processes. By definition, error in timing cannot
    be negative. As astutely noted by the Stackoverflow answer below:

        There's no way to get negative error because a computer can't ever
        compute faster than it can compute!

    Since all timing error is positive, the minimum timing is the timing that
    exhibits minimum error. For all intents and purposes, all other timings are
    effectively irrelevant.

    Parameters
    ----------
    callable : CallableTypes
        Callable to be timed.
    args : optional[SequenceTypes]
        Sequence of positional arguments to pass to each call to this callable.
        Defaults to the empty tuple.
    kwargs: optional[MappingType]
        Dictionary of keyword arguments to pass to each call to this callable.
        Defaults to the empty dictionary.
    repetitions : optional[int]
        Number of times to repeat each number of times to call this callable.
        Defaults to a reasonably small value.
    iterations: optional[int]
        Number of times to call this callable for each repetition. The total
        number of calls is thus given by `repetitions * iterations`. Defaults to
        a reasonably large value.

    See Also
    ----------
    https://stackoverflow.com/a/24105845/2809027
        Stackoverflow answer strongly inspiring this implementation.
    '''

    # Partial function binding this callable to these arguments. Avoid
    # defaulting unpassed positional and keyword arguments to empty data
    # structures, as doing so appears to substantially skew timings and hence
    # should be avoided if feasible.
    #
    # If both positional and keyword arguments were passed, bind this callable
    # to both.
    callable_bound = None
    if args and kwargs:
        callable_bound = partial(callable, *args, **kwargs)
    # Else if only positional arguments were passed, bind this callable to only
    # positional arguments.
    elif args:
        callable_bound = partial(callable, *args)
    # Else if only keyword arguments were passed, bind this callable to only
    # keyword arguments.
    elif kwargs:
        callable_bound = partial(callable, **kwargs)
    # Else, call this callable as is.
    else:
        callable_bound = callable

    # Return the minimum timing of this callable repeated this number of times.
    return min(Timer(callable_bound).repeat(repetitions, iterations))
