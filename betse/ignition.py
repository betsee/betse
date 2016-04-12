#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
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
from betse.util.io import loggers
from betse.util.py import pythons

# ....................{ INITIALIZERS                       }....................
def init() -> None:
    '''
    Initialize the current application.

    Specifically:

    * Validate core directories and files required at program startup, creating
      all such directories and files that do *not* already exist and are
      feasible to be created.
    * Create the root logger and default handlers for such logger (e.g.,
      printing informational messages to standard output, warnings and errors to
      standard error, and appending debug messages to the default logfile).

    To support caller-specific error handling, this function is intended to be
    called immediately *after* such application begins catching otherwise
    uncaught exceptions.
    '''

    # Validate core directories and files required at program startup.
    pathtree.init()

    #FIXME: After refactoring logging to merely log to standard file handles,
    #call this function *BEFORE* anything else, including pathtree.init().

    # Configure logging *AFTER* creating these directories, as such logging
    # writes to files in such directories.
    loggers.config.init(filename=pathtree.LOG_DEFAULT_FILENAME)
    # self._logger.error('ERROR!')
    # self._logger.warning('WARNING!')
    # self._logger.info('INFO!')
    # self._logger.debug('DEBUG!')

    # Validate mandatory dependencies *AFTER* configuring logging, thereby
    # logging exceptions raised by this validation.
    libs.init()

    # Validate the active Python interpreter *AFTER* validating mandatory
    # dependencies. While the former (mostly) comprises unenforced
    # recommendations, the latter comprises enforced requirements and hence is
    # performed first.
    pythons.init()
