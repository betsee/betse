#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level **Python profiling** (i.e., measuring various metrics pertaining to
Python code, including time and space performance) facilities.
'''

# ....................{ IMPORTS                            }....................
from betse.exceptions import BetseMethodUnimplementedException
from betse.util.io.log import logs
from betse.util.type.types import (
    type_check, CallableTypes, MappingType, SequenceTypes,)
from cProfile import Profile
from enum import Enum
from functools import partial
from io import StringIO
from pstats import Stats
from timeit import Timer

# ....................{ ENUMS                              }....................
ProfileType = Enum('ProfileType', ('NONE', 'CALL', 'SIZE',))
# ProfileType = Enum('ProfileType', ('NONE', 'CALL', 'LINE', 'SIZE',))
'''
Enumeration of all possible types of profiling supported by the
:func:`profile_callable` function.

Attributes
----------
NONE : enum
    Null profiling. Callables profiled under this type of profiling will still
    be called as expected but will *not* be profiled.
CALL : enum
    Call-granularity profiling of callables (e.g., functions, lambdas, methods,
    and `eval` statements). This level of granularity is more coarse-grained
    than that of the ``LINE`` type.
SIZE : enum
    Memory profiling of top-level objects returned by top-level callables,
    typically instances of the
    :class:`betse.science.phase.simphaser.SimPhase` class returned by public
    methods of the :class:`betse.science.simrunner.SimRunner` class. This
    profiling type requires installation of the optional third-party dependency
    :mod:`pympler`.
'''
# LINE : enum
#     Line-granularity profiling of callables (e.g., functions, lambdas, methods,
#     and `eval` statements). This level of granularity is more fine-grained than
#     that of the ``CALL`` type, but requires installation of the optional
#     third-party dependency :mod:`line_profiler`.

# ....................{ PROFILERS                          }....................
@type_check
def profile_callable(
    call: CallableTypes,
    args: SequenceTypes = None,
    kwargs: MappingType = None,
    is_profile_logged: bool = True,
    profile_filename: str = None,
    profile_type: ProfileType = ProfileType.CALL,
) -> object:
    '''
    Profile the passed callable with the passed positional and keyword arguments
    (if any), returning the value returned by this call and optionally logging
    and serializing the resulting profile to the file with the passed filename.

    Parameters
    ----------
    call : CallableTypes
        Callable to be profiled.
    args : optional[SequenceTypes]
        Sequence of positional arguments to pass to the call to this callable.
        Defaults to `None`, in which no such arguments are passed.
    kwargs: optional[MappingType]
        Dictionary of keyword arguments to pass to the call to this callable.
        Defaults to `None`, in which no such arguments are passed.
    is_profile_logged : optional[bool]
        `True` only if logging the profile of this call *after* this call.
        Defaults to `True`.
    profile_filename : optional[str]
        Absolute or relative path of the file to serialize this profile to. If
        this file already exists, this file will be silently overwritten
        _without_ explicit warning or error. (While silently overwriting is
        typically discouraged, doing so is sensible in this case. Profiles are
        frequently and trivially created, rendering their preservation
        insignificant by compare to usability concerns.)
        Defaults to `None`, in which case no such file is serialized.
    profile_type : optional[ProfileType]
        Type of profiling to perform. Defaults to :data:`ProfileType.CALL`, in
        which case call-granularity profiling is performed.

    Returns
    ----------
    object
        Value returned by this call.
    '''

    # Private module function performing this type of profiling. Due to type
    # checking, this type is guaranteed to be a valid key of this dictionary.
    profiler = _PROFILE_TYPE_TO_PROFILER[profile_type]

    # Default unpassed positional and keyword arguments.
    if   args is None:   args = ()
    if kwargs is None: kwargs = {}

    # Profile a call of this callable with this profiler, returning the value
    # that this call returns.
    return profiler(
        call=call,
        args=args,
        kwargs=kwargs,
        is_profile_logged=is_profile_logged,
        profile_filename=profile_filename,
    )

# ....................{ PROFILERS ~ none                   }....................
def _profile_callable_none(
    call, args, kwargs, is_profile_logged, profile_filename) -> object:
    '''
    Call the passed callable with the passed positional and keyword arguments
    *without* profiling this call, returning the value returned by this call.

    See Also
    ----------
    :func:`profile_callable`
        Further details on function signature.
    '''

    return call(*args, **kwargs)

# ....................{ PROFILERS ~ call                   }....................
def _profile_callable_call(
    call, args, kwargs, is_profile_logged, profile_filename) -> object:
    '''
    Profile the passed callable in a call-oriented deterministic manner with the
    passed positional and keyword arguments (if any), returning the value
    returned by this call and optionally logging and serializing the resulting
    profile to the file with the passed filename.

    See Also
    ----------
    :func:`profile_callable`
        Further details on function signature.
    '''

    # Log this fact.
    logs.log_debug('Call-granularity profiling enabled.')

    # Call-granularity profile of the subsequent call to this callable.
    profile = Profile()

    # Value returned by calling this callable with these arguments, storing a
    # profile of this call into this "profile" object.
    return_value = profile.runcall(call, *args, **kwargs)

    # If the caller requested this profile be logged...
    if is_profile_logged:
        # Maximum number of slowest callables to log, arbitrarily defined to be
        # twice the default number of rows in the average Linux terminal.
        CALLABLES_MAX = 48

        # String buffer describing these slowest callables.
        calls_sorted = StringIO()

        # Statistics harvested from this profile into this string buffer.
        calls = Stats(profile, stream=calls_sorted)

        # For readability, strip the dirnames from all pathnames in statistics
        # output, reducing these pathnames to basenames.
        calls.strip_dirs()

        # Sort all profiled callables by cumulative time (i.e., total time spent
        # in a callable including all time spent in calls to callables called by
        # that callable).
        calls.sort_stats('cumtime')

        # Write the slowest sublist of these callables to this string buffer.
        calls.print_stats(CALLABLES_MAX)

        # Sort all profiled callables by "total" time (i.e., total time spent in
        # a callable excluding all time spent in calls to callables called by
        # that callable).
        calls.sort_stats('tottime')

        # Write the slowest sublist of these callables to this string buffer.
        calls.print_stats(CALLABLES_MAX)

        # Log this string buffer.
        logs.log_info(
            'Slowest %d callables profiled by both '
            'cumulative and total time:\n%s',
            CALLABLES_MAX, calls_sorted.getvalue())

    # If the caller requested this profile be serialized to a file...
    if profile_filename is not None:
        # Log this serialization.
        logs.log_info(
            'Writing Python-formatted profile to "%s".', profile_filename)

        # Serialize this profile to this file.
        profile.dump_stats(profile_filename)

    # Return the value returned by this call.
    return return_value

# ....................{ PROFILERS ~ line                   }....................
#FIXME: Implement this function properly. Sadly, Python does *NOT* provide an
#out-of-the-box solution for line-based profiling. To do so, either (...both?)
#of the following third-party packages will need to be dynamically detected,
#imported, and leveraged:
#
#* "line", profiling in a line- rather than call-based manner.
#  * "line_profiler", a third-party C extension profiling each line (rather than
#    function as cProfile does). Infrequently updated, but slightly more
#    frequently than "statprof", which is effectively dead. Sadly, the
#    "line_profiler" API is less than ideal, requiring that each individual
#    functions to be non-transitively line-profiled be explicitly decorated with
#    a "line_profiler"-specific decorator -- which, frankly, is crazy and
#    probably not reasonably supportable. See the live repository at:
#    https://github.com/rkern/line_profiler
#  * "statprof", a third-party C extension operating rather differently than
#    either "lineprof" or cProfile. Rather than deterministically instrumenting
#    each line or function call (respectively), "statprof" non-deterministically
#    wakes up at predefined intervals, records a stack trace, and then goes back
#    to sleep. On application completion, "statprof" then tallies up each stack
#    trace and outputs a command-line table of the most expensive lines. Pretty
#    sweet idea. Unsurprisingly, it also appears to be the fastest profiler.
#    Sadly, "statprof" was last updated in 2015; see the dead repository at:
#    https://github.com/bos/statprof.py
#
#Note that the "pprofile" package *CANNOT* be supported, due to its encumberment
#under a GPL-2 license.
#
#In short, there currently appears to be no BSD-compatible, well-maintained line
#profiling solution for Python. Hence, this function is disabled.

def _profile_callable_line(
    call, args, kwargs, is_profile_logged, profile_filename) -> object:
    '''
    Profile the passed callable in a line-oriented deterministic manner with the
    passed positional and keyword arguments (if any), returning the value
    returned by this call and optionally logging and serializing the resulting
    profile to the file with the passed filename.

    See Also
    ----------
    :func:`profile_callable`
        Further details on function signature.
    '''

    # Log this fact.
    logs.log_debug('Line-granularity profiling enabled.')

    # Raise an exception unless this optional dependency is available.
    # libs.die_unless_runtime_optional('????')

    #FIXME: Excise this after properly implementing this function.
    raise BetseMethodUnimplementedException()

# ....................{ PROFILERS ~ size                   }....................
def _profile_callable_size(
    call, args, kwargs, is_profile_logged, profile_filename) -> object:
    '''
    Profile the space consumed by the value returned by calling the passed
    callable with the passed positional and keyword arguments (if any) and
    return this value.

    See Also
    ----------
    :func:`profile_callable`
        Further details on function signature.
    '''

    # Avoid circular import dependencies.
    from betse.lib import libs
    from betse.util.type.obj import objsizes

    # Log this fact.
    logs.log_debug('Memory profiling enabled.')

    # Raise an exception unless the optional "pympler" dependency is available.
    libs.die_unless_runtime_optional('pympler')

    # Value returned by calling this callable with these arguments.
    return_value = call(*args, **kwargs)

    # Log the profiling to be subsequently performed. Since doing so can recurse
    # through the full object tree in the worst case, inform the end user of
    # this potentially slow operation *BEFORE* doing so.
    logs.log_debug(
        'Profiling memory for %r() return value of type %r...',
        call, type(return_value))

    # Human-readable string synopsizing this object's memory consumption.
    size_profile = objsizes.get_size_profile(
        obj=return_value,

        # Maximum depth of the recursion tree induced by calling this function.
        # To avoid combinatorial explosion, this maximum is necessarily capped
        # to a small positive integer.
        vars_depth=2,

        # Maximum number of the largest instance variables of the value returned
        # by calling this callable to log, arbitrarily defined to be twice the
        # default number of rows in the average Linux terminal.
        vars_max=48,
    )

    # Log this synopsis.
    logs.log_info(
        'Largest top-level object variables '
        'profiled by total space consumption:\n\n'
        'Top-level object %s', size_profile)

    #FIXME: Implement support for serializing this profile to disk.
    # # If the caller requested this profile be serialized to a file...
    # if profile_filename is not None:
    #     # Notify the caller that doing so is currently unsupported..
    #     logs.log_warning(
    #         'Memory profile not saved to "%s", '
    #         'as this feature is currently unimplemented.',
    #         profile_filename)

    # Return the value returned by this call.
    return return_value

# ....................{ GLOBALS ~ private                  }....................
# Technically, the same effect is also achievable via getattr() on the current
# module object. Doing so is complicated by artificial constraints Python
# imposes on doing so (e.g., obtaining the current module object is obscure) and
# the PEP 20 doctrine of "Explicit is better than implicit." We beg to disagree.
# Nonetheless, the explicit approach remains preferable in this edge-case.
_PROFILE_TYPE_TO_PROFILER = {
    ProfileType.CALL: _profile_callable_call,
    # ProfileType.LINE: _profile_callable_line,
    ProfileType.NONE: _profile_callable_none,
    ProfileType.SIZE: _profile_callable_size,
}
'''
Dictionary mapping from each supported type of profiling to the module function
performing this type of profiling.

This private global is intended for use _only_ by the public
:func:`profile_callable` function.
'''

# ....................{ TIMERS                             }....................
@type_check
def time_callable(
    call: CallableTypes,
    args: SequenceTypes = None,
    kwargs: MappingType = None,
    repetitions: int = 7,
    iterations: int = 1000,
) -> float:
    '''
    Time the passed callable with the passed positional and keyword arguments
    (if any) called the passed number of repetitions of the passed number of
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
    call : CallableTypes
        Callable to be timed.
    args : optional[SequenceTypes]
        Sequence of positional arguments to pass to each call to this callable.
        Defaults to `None`, in which no such arguments are passed.
    kwargs: optional[MappingType]
        Dictionary of keyword arguments to pass to each call to this callable.
        Defaults to `None`, in which no such arguments are passed.
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
        callable_bound = partial(call, *args, **kwargs)
    # Else if only positional arguments were passed, bind this callable to only
    # positional arguments.
    elif args:
        callable_bound = partial(call, *args)
    # Else if only keyword arguments were passed, bind this callable to only
    # keyword arguments.
    elif kwargs:
        callable_bound = partial(call, **kwargs)
    # Else, call this callable as is.
    else:
        callable_bound = call

    # Return the minimum timing of this callable repeated this number of times.
    return min(Timer(callable_bound).repeat(repetitions, iterations))
