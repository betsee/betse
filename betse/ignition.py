#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level application initialization common to both the CLI and GUI.
'''

#FIXME: Merge both this and the comparable "betsee.guiignition" submodule into
#some combination of the "betse.metaapp", "betse.util.app.meta.metaappabc", and
#"betsee.guimetaapp" submodules as follows:
#
#* Shift the instantiation of application singletons (e.g., calls to the
#  BetseAppMeta() and BetseeAppMeta() constructors) into a single
#  application-specific point: the top-level "__main__" submodule for each
#  application. Specifically:
#  * "betse.__main__" for BETSE.
#  * "betsee.__main__" for BETSEE.
#* Remove the "betse.ignition._IS_INITTED" global entirely. The mere existence
#  of an "app_meta" singleton should now deterministically imply that the
#  application as a whole has been initialized.
#* Refactor all calls to the betse.ignition.init() function to the
#  betse.util.app.metaapp.metaappton.set_app_meta() function instead.
#* Shift the init() function into the "MetaAppABC" superclass as a public
#  method of the same name.
#* Omit all of the following from the new MetaAppABC.init() method:
#  * The passed "app_meta" parameter, which will now be handled elsewhere.
#* Of course, this somewhat suggests that we don't necessarily need a
#  MetaAppABC.init() method at all but that the contents of that method could
#  simply be merged into the MetaAppABC.__init__(). For now, however, let's
#  continue with this possibly overkill approach and see where we end up.
#  * *WAIT.* Actually, we do absolutely still want an init() method to help
#    manage object lifecycle with fine-grained control. What we don't
#    necessarily need is a *PUBLIC* init() method. A private _init() method
#    implicitly called by the __init__() method could be just the thing, as
#    this approach would still enable subclass implementations but would
#    simplify the object creation process for callers -- always a useful thing:
#    e.g.,
#
#    def __init__(self) -> None:
#
#        super().__init__()
#
#        # It is a good idea.
#        self._init()
#
#        # It is a great idea!
#        metappton.set_app_meta(self)
#  * *WAIT.* Everything above is great, except privatizing the init() method
#    part. We probably want to keep that public. The CLI, in particular, is
#    likely to need public access to this method.
#* Follow suite with the reinit() function... or perhaps simply excise this
#  function, which we no longer appear to require or desire.
#* Refactor the ignite() function into a new MetaAppABC.ignite() method.

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

#FIXME: The "~/.betse" directory grows fairly large fairly quickly. It'd be
#great to emit non-fatal warnings if its size exceeds some reasonable threshold
#(e.g., 1MB).

# ....................{ IMPORTS                           }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To defer heavyweight and possibly circular imports, the top-level of
# this module may import *ONLY* from standard Python packages. All imports from
# application and third-party packages should be deferred to their point of
# use.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# ....................{ GLOBALS                           }....................
_IS_INITTED = False
'''
``True`` only if the :func:`init` function has already been called.

That function uses this private boolean to guard against repeated invocations
of the :func:`init` function from multiple modules in the same Python process
(e.g., :mod:`betse.science.__init__`, :mod:`betse.util.cli.cliabc`). While that
function does technically support repeated calls, each additional call after
the first inefficiently performs no meaningful work and is thus safely
ignorable.
'''

# ....................{ IGNITERS                          }....................
def ignite(*args, **kwargs) -> None:
    '''
    Initialize both the current application *and* all mandatory third-party
    dependencies of this application with sane defaults.

    This high-level convenience function intentionally provides *no* means of
    initializing either this application or these dependencies with alternative
    parameters. To do so, callers should call all lower-level initialization
    functions directly (e.g., :func:`init`, :func:`betse.lib.libs.init`).

    Parameters
    ----------
    All passed parameters are passed as is to the :func:`init` function.
    '''

    # Defer heavyweight and possibly circular imports.
    from betse.lib import libs

    # Initialize this application.
    init(*args, **kwargs)

    # Initialize all dependencies *AFTER* initializing this application.
    libs.init()

# ....................{ INITIALIZERS                      }....................
def reinit(*args, **kwargs) -> None:
    '''
    (Re-)initialize this application -- but *not* mandatory third-party
    dependencies of this application, which requires external resources (e.g.,
    command-line options, configuration files) to do so.

    Specifically:

    * If this application has *not* already been initialized under the active
      Python process, this function initializes this application.
    * Else, this application has already been initialized under the active
      Python process. In this case, this function reinitializes this
      application.

    Parameters
    ----------
    All passed parameters are passed as is to the :func:`init` function.

    See Also
    ----------
    :func:`betse.lib.libs.reinit`
        Function (re)-initializing all mandatory third-party dependencies.
    '''

    # Defer heavyweight and possibly circular imports.
    from betse.util.app.meta import metaappton

    # Force the init() function to reinitialize this application.
    global _IS_INITTED
    _IS_INITTED = False

    # Unset the application metadata singleton *BEFORE* the subsequent call to
    # the init() function below sets this singleton. If this is *NOT* done,
    # that call sensibly raises an exception.
    metaappton.unset_app_meta()

    # Reinitialize this application with all passed parameters.
    init(*args, **kwargs)


# To defer heavyweight and possibly circular imports, type-checking of this
# function parameter is deferred to the subsequent call to the
# betse.metaapp.init() function performed internally by this function.
def init(app_meta = None) -> None:
    '''
    Initialize the current application if this application has not already been
    initialized under the active Python process *or* noop otherwise.

    Specifically, this function (in order):

    #. Enables Python's standard handler for segmentation faults.
    #. Globalizes the passed application metadata singleton.
    #. Enables this application's default logging configuration.
    #. Validates (but does *not* initialize) all mandatory third-party
       dependencies of this application, which the :func:`betsee.lib.libs.init`
       function initializes independently.
    #. Validates core directories and files required at program startup,
       creating all such directories and files that do *not* already exist and
       are reasonably creatable.
    #. Validates the active Python interpreter (e.g., to support
       multithreading).
    #. Validates the underlying operating system (e.g., to *not* be a vanilla
       Windows shell environment ala either CMD.exe or PowerShell).

    To support caller-specific error handling, this function is intended to be
    called immediately *after* this application begins catching otherwise
    uncaught exceptions.

    Parameters
    ----------
    app_meta : BetseMetaAppOrNoneTypes
        Caller-specific application metadata singleton (i.e., instance of the
        :class:`betse.metaapp.BetseMetaApp` subclass). Defaults to ``None``, in
        which case this parameter defaults to a vanilla instance of that
        subclass.
    '''

    # Defer heavyweight and possibly circular imports.
    from betse.lib import libs
    from betse.metaapp import BetseMetaApp
    from betse.util.app.meta import metaappton
    from betse.util.io.error import errfault
    from betse.util.io.log import logconfig, logs
    from betse.util.os import oses
    from betse.util.py import pys
    from betse.util.type.obj import objects

    # If this function has already been called...
    global _IS_INITTED
    if     _IS_INITTED:
        # Log this noop. Since this function has already been called, this
        # log message is guaranteed to be logged.
        logs.log_debug('Application reinitialization ignored.')

        # Reduce to a noop.
        return

    # Enable Python's standard handler for segmentation faults *BEFORE*
    # performing any further logic, any of which could conceivably trigger a
    # segmentation fault and hence process termination.
    errfault.handle_faults()

    # If passed no caller-specific singleton, default to a generic singleton.
    if app_meta is None:
        app_meta = BetseMetaApp()

    # Globalizes the passed application metadata singleton *BEFORE* performing
    # any further logic, any of which could conceivably require this singleton.
    # Indeed, the subsequent logconfig.init() call does exactly that.
    metaappton.set_app_meta(app_meta)

    # Enable our default logging configuration for the current Python process
    # *BEFORE* performing any validation, thus logging any exceptions raised by
    # this validation.
    logconfig.init()

    # Log all prior behaviour. Attempting to do so *BEFORE* enabling our
    # default logging configuration above would silently fail, since the
    # standard "logging" API silently squelches debug messages by default.
    logs.log_debug('Default segementation fault handler enabled.')
    logs.log_debug('Application singleton "%s" established.',
        objects.get_class_name_unqualified((app_meta)))

    # Validate mandatory dependencies. Avoid initializing these dependencies
    # here (e.g., by calling libs.init()), which requires the logging
    # configuration to have been finalized (e.g., by parsing CLI options),
    # which has yet to occur this early in the application lifecycle.
    libs.die_unless_runtime_mandatory_all()

    # Validate the active Python interpreter and operating system *AFTER*
    # mandatory dependencies. While the former (mostly) comprises unenforced
    # recommendations, the latter comprises enforced requirements and should
    # thus be validated first.
    oses.init()
    pys.init()

    # Record this function as having been called *AFTER* successfully doing so.
    _IS_INITTED = True
