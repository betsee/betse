#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level application initialization common to both the CLI and GUI.
'''

#FIXME: Raise an exception if running with superuser privelages. To avoid
#inevitable security issues, Under *NO* circumstances should either BETSE or
#BETSEE ever be run with elevated permissions. To do so, see the following
#canonical answer:
#    https://stackoverflow.com/a/1026626/2809027

#FIXME: Print memory-related metadata when running "betse info" *AND* log
#non-fatal warnings when BETSE is run under a low-memory environment (e.g.,
#< 4GB available free memory). To do so, note the following canonical API:
#
#    psutil.Process(os.getpid()).get_memory_info()

#FIXME: Consider replacing bottleneck Numpy routines with routines imported from
#the following third-party Numpy-like frameworks:
#
#* "bottleneck", providing optimized routines accepting Numpy arrays --
#  implemented in Cython and hence faster than comparible Numpy routines.
#* "numexpr", a Theano-like framework accepting Numpy arrays -- performing
#  CPU-centric parallelization of expensive array operations. Whereas Theano
#  permits such operations to be conveniently expressed in pure-Python, however,
#  numexpr inconveniently requires such operations be expressed as... wait for
#  it, raw strings. So, that sucks. Nonetheless, worth a possible look.
#* "blaze", a purported second-gen Numpy replacement. We harbour sincere doubts,
#  but everything deserves its millisecond to shine in the light. Ah. We see.
#  Blaze is considerably more heavy-weight than Numpy, and largely serves a
#  completely different marketshare: supercomputing. That's well beyond our
#  means, at the moment. Numpy it is!
#
#In short, "bottleneck" is probably the only framework listed above of interest.

#FIXME: Consider optimizing frequently used matrix and vector computations with
#Theano, a general-purpose Python mathematical optimization framework. One
#particularly compelling use case for Theano is to portably distribute
#computational work across multiple GPUs. In general, Theano can be used to
#reduce arbitrarily complex symbolic expressions expressed in pure Python to
#dynamically compiled machine code on-the-fly. Fairly amazing, all around. For
#the high-level synopsis, see:
#
#    http://deeplearning.net/software/theano/introduction.html
#FIXME: Theano and Torch (a similar framewark) appear to now be subsumed by
#TensorFlow, a Google-backed framework originally implemented in support of
#machine learning workflows at Google (e.g., DeepMind), but sufficiently
#generalized as to support a wide variety of computational needs -- like ours.

#FIXME: The "~/.betse" directory grows fairly large fairly quickly. It'd be
#great to emit non-fatal warnings if its size exceeds some reasonable threshold
#(e.g., 1MB).

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To defer heavyweight and possibly circular imports, the top-level of
# this module may import *ONLY* from standard Python packages. All imports from
# application and third-party packages should be deferred to their point of use.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# ....................{ GLOBALS                            }....................
_IS_IGNITED = False
'''
``True`` only if the :func:`init` function has already been called.

That function uses this private boolean to guard against repeated invocations of
the :func:`init` function from multiple modules in the same Python process
(e.g., :mod:`betse.science.__init__`, :mod:`betse.util.cli.cliabc`). While that
function does technically support repeated calls, each additional call after the
first inefficiently performs no meaningful work and is thus safely ignorable.
'''

# ....................{ IGNITERS                           }....................
def ignite() -> None:
    '''
    Initialize both the current application *and* all mandatory third-party
    dependencies of this application with sane defaults.

    This high-level convenience function intentionally provides *no* means of
    initializing either this application or these dependencies with alternative
    parameters. To do so, callers should call all lower-level initialization
    functions directly (e.g., :func:`init`, :func:`betse.lib.libs.init`).
    '''

    # Defer heavyweight and possibly circular imports.
    from betse.lib import libs

    # Initialize this application.
    init()

    # Initialize all dependencies *AFTER* initializing this application.
    libs.init()

# ....................{ INITIALIZERS                       }....................
def reinit() -> None:
    '''
    (Re-)initialize this application -- but *not* mandatory third-party
    dependencies of this application, which requires external resources (e.g.,
    command-line options, configuration files) to be parsed.

    Specifically:

    * If this application has *not* already been initialized under the active
      Python process, this application will be initilialized.
    * Else, this application has already been initialized under the active
      Python process. In this case, this application will be re-initilialized.

    See Also
    ----------
    :func:`betse.lib.libs.reinit`
        Function (re)-initializing all mandatory third-party dependencies.
    '''

    # Force the init() function to reinitialize this application.
    global _IS_IGNITED
    _IS_IGNITED = False

    # Reinitialize this application.
    init()


#FIXME: Update docstring when time and kindness affords.
def init() -> None:
    '''
    Initialize the current application if this application has not already been
    initialized under the active Python process *or* noop otherwise.

    Specifically, this function (in order):

    #. Validates core directories and files required at program startup,
       creating all such directories and files that do *not* already exist and
       are reasonably creatable.
    #. Validates but does *not* initialize mandatory third-party dependencies of
       this application, which must be initialized independently by the
       :func:`betsee.lib.libs.init` function.

    To support caller-specific error handling, this function is intended to be
    called immediately *after* this application begins catching otherwise
    uncaught exceptions.
    '''

    # If this function has already been called, noop.
    global _IS_IGNITED
    if     _IS_IGNITED:
        return

    # Defer heavyweight and possibly circular imports.
    from betse.lib import libs
    from betse.util.io.log import logconfig
    from betse.util.py import pys

    # Enable the default logging configuration for the current Python process
    # *BEFORE* performing any validation, thus logging any exceptions raised by
    # this validation.
    logconfig.init()

    # Validate mandatory dependencies. Avoid initializing these dependencies
    # here (e.g., by calling libs.init()), which requires the logging
    # configuration to have been finalized (e.g., by parsing CLI options), which
    # has yet to occur this early in the application lifecycle.
    libs.die_unless_runtime_mandatory_all()

    # Validate the active Python interpreter *AFTER* mandatory dependencies.
    # While the former (mostly) comprises unenforced recommendations, the latter
    # comprises enforced requirements and should thus be validated first.
    pys.init()

    # Record this function as having been called *AFTER* successfully doing so.
    _IS_IGNITED = True
