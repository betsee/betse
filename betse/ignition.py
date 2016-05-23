#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
High-level application initialization common to both the CLI and GUI.
'''

#FIXME: Print a non-fatal warning if Numpy is linked against a non-multithreaded
#BLAS implementation. Also, print the name of the BLAS implementation against
#which Numpy is linked with similar "betse info" metadata.
#
#Indeed! It would appear that the metadata we require (...and more!) is
#available via the "numpy.__config__" API. Sure, it's private, but well-
#established at this point. We can't really see it going away. See also:
#
#    https://stackoverflow.com/a/19350234/2809027
#FIXME: Phenomenal Python code for profiling system metadata, including BLAS,
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
# Defer heavyweight imports to their point of use.
from betse import pathtree
from betse.lib import libs
from betse.util.py import pys

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

    # Validate core directories and files required at program startup.
    pathtree.init()

    # Validate mandatory dependencies *AFTER* configuring logging, thereby
    # logging exceptions raised by this validation.
    libs.init()

    # Validate the active Python interpreter *AFTER* validating mandatory
    # dependencies. While the former (mostly) comprises unenforced
    # recommendations, the latter comprises enforced requirements and hence is
    # performed first.
    pys.init()

    # Record this function as already called *AFTER* successfully calling all
    # prior initialization.
    _IS_IGNITED = True
