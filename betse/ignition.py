#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
High-level application initialization common to both the CLI and GUI.
'''

#FIXME: Implement a new CLI option "-p" and "--profile" enabling CLI-driven
#profiling. Python provides a phenomenal C-based profiler out-of-the-box named
#"cProfile", an importable C extension emitting a cProfile dump when the
#profile application halts as well as a command-line table of time spent in each
#application function sorted via a variety of metrics. It's pretty much perfect.
#A variety of upstream third-party utilities then exist to visualize cProfile
#dumps, although the command-line table should suffice for a low-hanging fruit
#implementation. To do so, it appears to suffice to:
#
#    import cProfile
#    cProfile.run('main()', 'cprofile.dmp')
#
#After implementing cProfile support, consider also implementing support for:
#
#* "lineprof", a third-party C extension profiling each line (rather than
#  function as cProfile does). Basically, cProfile on metric steroids.
#* "statprof", a third-party C extension operating rather differently than
#  either "lineprof" or cProfile. Rather than deterministically instrumenting
#  each line or function call (respectively), "statprof" non-deterministically
#  wakes up at predefined intervals, records a stack trace, and then goes back
#  to sleep. On application completion, "statprof" then tallies up each stack
#  trace and outpus a command-line table of the most expensive lines. Pretty
#  sweet idea. Unsurprisingly, it also appears to be the fastest profiler.

#FIXME: Print a non-fatal warning if Numpy is linked against a non-multithreaded
#BLAS implementation. Also, print the name of the BLAS implementation against
#which Numpy is linked with similar "betse info" metadata.
#
#Indeed! It would appear that the metadata we require (...and more!) is
#available via the "numpy.__config__" API. Sure, it's private, but well-
#established at this point. We can't really see it going away. See also:
#
#    https://stackoverflow.com/a/19350234/2809027
#FIXME: Phenomenal Python code for obtaining system metadata, including BLAS,
#LAPACK, Atlas, OpenBLAS, and numpy-specific system metadata:
#
#    https://gist.github.com/sandys/258707dae9b79308594b#file-system_info-py

#FIXME: Consider optimizing frequently used matrix and vector computations with
#Theano, a general-purpose Python mathematical optimization framework. One
#particularly compelling use case for Theano is to portably distribute
#computational work across multiple GPUs. In general, Theano can be used to
#reduce arbitrarily complex symbolic expressions expressed in pure Python to
#dynamically compiled machine code on-the-fly. Fairly amazing, all around. For
#the high-level synopsis, see:
#
#    http://deeplearning.net/software/theano/introduction.html

# ....................{ IMPORTS                            }....................

# ....................{ GLOBALS                            }....................
_IS_IGNITED = False
'''
`True` only if the `init()` function has already been called.

That function uses this private boolean to guard against repeated invocations of
that function from multiple modules  in the same Python process (e.g.,
`betse.science.__init__`, `betse.cli.cliabc`). While that function does
technically support repeated invocations, each additional invocation after the
first inefficiently performs no meaningful work -- and is thus safely ignorable.
'''

# ....................{ INITIALIZERS                       }....................
#FIXME: Update docstring when time and kindness affords.
def init() -> None:
    '''
    Initialize the current application if this application has not already been
    initialized by the current Python process _or_ noop otherwise.

    Specifically, this function:

    * Validates core directories and files required at program startup, creating
      all such directories and files that do _not_ already exist and are
      reasonably creatable.

    To support caller-specific error handling, this function is intended to be
    called immediately _after_ this program begins catching otherwise uncaught
    exceptions.
    '''

    # If this function has already been called, return silently.
    global _IS_IGNITED
    if     _IS_IGNITED:
        return

    # Defer heavyweight imports to their point of use.
    from betse import pathtree
    from betse.lib import libs
    from betse.util.io.log import logconfig
    from betse.util.py import pys

    # Enable the default logging configuration for the current Python process
    # *BEFORE* performing any validation, thus logging any exceptions raised by
    # this validation.
    logconfig.init()

    # Validate core directories and files required at program startup.
    pathtree.init()

    # Validate mandatory dependencies.
    libs.init()

    # Validate the active Python interpreter *AFTER* validating mandatory
    # dependencies. While the former (mostly) comprises unenforced
    # recommendations, the latter comprises enforced requirements and hence is
    # performed first.
    pys.init()

    # Record this function as already called *AFTER* successfully calling all
    # prior initialization.
    _IS_IGNITED = True
