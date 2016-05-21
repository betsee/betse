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
from betse import pathtree
from betse.lib import libs
from betse.util.py import pys

# ....................{ INITIALIZERS                       }....................
#FIXME: Update docstring when time and kindness affords.
def init() -> None:
    '''
    Initialize the current application.

    Specifically, this method:

    * Validates core directories and files required at program startup, creating
      all such directories and files that do _not_ already exist and are
      reasonably creatable.

    To support caller-specific error handling, this function is intended to be
    called immediately _after_ this program begins catching otherwise uncaught
    exceptions.
    '''

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
