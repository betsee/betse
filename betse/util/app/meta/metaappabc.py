#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2017-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level **application metadata singleton** (i.e., application-wide object
synopsizing application metadata via read-only properties) hierarchy.
'''

#FIXME: The current approach is inefficient in the case of BETSE being
#installed as a compressed EGG rather than an uncompressed directory. In the
#former case, the current approach (namely, the call to
#resources.get_pathname() performed below) silently extracts the entirety of
#this egg to a temporary setuptools-specific cache directory. That's bad. To
#circumvent this, we'll need to refactor the codebase to directly require only
#"file"-like objects rather than indirectly requiring the absolute paths of
#data resources that are then opened as "file"-like objects.
#
#Specifically, whenever we require a "file"-like object for a codebase
#resource, we'll need to call the setuptools-specific
#pkg_resources.resource_stream() function rather than attempting to open the
#path given by a global below. Ultimately, *ALL* of the codebase-specific
#globals declared below (e.g., "DATA_DIRNAME") should go away.

# ....................{ IMPORTS                           }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To avoid race conditions during setuptools-based installation, this
# module may import *ONLY* from modules guaranteed to exist at the start of
# installation. This includes all standard Python and application modules but
# *NOT* third-party dependencies, which if currently uninstalled will only be
# installed at some later time in the installation. Likewise, to avoid circular
# import dependencies, the top-level of this module should avoid importing
# application modules where feasible.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

from abc import ABCMeta
from betse.exceptions import BetseGitException
from betse.util.type.decorator.decmemo import property_cached
from betse.util.type.types import type_check, ModuleType, StrOrNoneTypes

# ....................{ SUPERCLASSES                      }....................
class MetaAppABC(object, metaclass=ABCMeta):
    '''
    Abstract base class of all **application metadata singleton** (i.e.,
    application-wide object synopsizing application metadata via read-only
    properties) subclasses.

    Pathnames
    ----------
    This superclass details the structure of this application on the local
    filesystem with properties providing the absolute pathnames of
    application-specific files and directories. These pathnames are intended
    for use by both this application and downstream consumers (e.g., BETSEE).
    For portability, these pathnames intentionally support standard and
    non-standard runtime environments -- including :mod:`setuptools`-installed
    script wrappers *and* PyInstaller-frozen executables.

    Caveats
    ----------
    **Neither this superclass nor any subclass of this superclass may be safely
    accessed from the top-level ``setup.py`` script of this or any other
    application.** This application and hence this superclass is *not*
    guaranteed to exist at setuptools-based installation-time for downstream
    consumers (e.g., BETSEE).

    For that same reason, this superclass intentionally avoids refactoring
    metadata constants defined by the top-level application modules (e.g.,
    :mod:`betse.metadata`, :mod:`betse.metadeps`) into object properties (e.g.,
    refactoring the :attr:`betse.metadata.NAME` constant into an object
    :meth:`name` property). Doing so would prevent their reuse from the
    top-level ``setup.py`` scripts defined by downstream consumers, which would
    render these constants effectively useless for their principal use case.

    Attributes
    ----------
    _is_libs_initted : bool
        ``True`` only if the :meth:`init_libs` method has already been called,
        enabling the optimal :meth:`init_libs_if_needed` method to silently
        reduce to a noop when ``True``.
    '''

    # ..................{ INITIALIZERS                      }..................
    def __init__(self, *args, **kwargs) -> None:
        '''
        Initialize both this application metadata singleton and the current
        application (except mandatory third-party dependencies of this
        application) if this is the first creation of such a singleton for this
        application *or* raise an exception otherwise (i.e., if this
        application has already instantiated such a singleton).

        Design
        ----------
        Callers should avoid explicitly calling either the
        :func:`betse.util.app.meta.metaappton.set_app_meta` setter or
        :meth:`init_sans_libs` method, unless absolutely required for
        inconceivable reasons. Why? Because this method already calls those
        callables. Specifically, to enforce the contractual guarantees that:

        * Only a single such singleton be instantiated for the lifetime of this
          application, this method internally calls the
          :func:`betse.util.app.meta.metaappton.set_app_meta` setter preserving
          this guarantee.
        * Instantiating this singleton suffices to initialize this application
          (except third-party dependencies of this application), this method
          internally calls the :meth:`init_sans_libs` method preserving
          this guarantee. For maintainability, callers should avoid doing so.

        Parameters
        ----------
        All parameters are passed as is to the :meth:`init_sans_libs` method.

        Raises
        ----------
        BetseMetaAppException
            If this is *not* the first instantiation of an application
            metadata singleton for the active Python interpreter.

        See Also
        ----------
        :meth:`init_sans_libs`
        :meth:`betse.util.app.meta.metaappton.set_app_meta`
            Further details.
        '''

        # Avoid circular import dependencies.
        from betse.util.app.meta import metaappton

        # Localize all instance variables.
        #
        # Note the init_libs() method to *NOT* have been called yet.
        self._is_libs_initted = False

        # Globalize this singleton *BEFORE* subsequent logic (e.g., the
        # logconfig.init() call performed by the self.init() call), any of
        # which could potentially require this singleton.
        metaappton.set_app_meta(self)

        # Initialize this application.
        self.init_sans_libs(*args, **kwargs)


    @type_check
    def init_sans_libs(self) -> None:
        '''
        Initialize this application *except* mandatory third-party dependencies
        of this application, which requires external resources (e.g.,
        command-line options, configuration files) to have been parsed.

        Specifically, the default implementation of this method (in order):

        #. Enables Python's standard handler for segmentation faults.
        #. Globalizes the passed application metadata singleton.
        #. Enables this application's default logging configuration.
        #. Validates but does *not* initialize all mandatory third-party
           dependencies of this application, which the :meth:`init_libs` method
           does so subsequently.
        #. Validates core directories and files required at program startup,
           creating all such directories and files that do *not* already exist
           and are reasonably creatable.
        #. Validates the active Python interpreter (e.g., to support
           multithreading).
        #. Validates the underlying operating system (e.g., to *not* be a
           vanilla Windows shell environment ala either CMD.exe or PowerShell).

        To support caller-specific error handling, this function is intended to
        be called immediately *after* this application begins catching
        otherwise uncaught exceptions.

        Caveats
        ----------
        **This method is idempotent.** This method has been explicitly designed
        to be safely recallable. Each subsequent invocation of this method
        following the first simply reinitializes this application. While
        typically useless, idempotency is required by low-level automation to
        guarantee consistency across repeated runs (e.g., tests, scripts).

        Design
        ----------
        This method intentionally avoids initializing mandatory dependencies,
        as doing so would require late-time startup logic assumed by such
        initiazilation to have already been performed -- namely, finalization
        of the logging configuration and hence command-line argument parsing.
        By compare, this method is internally called by the :meth:`__init__`
        method called as the first statement of this application. Since no
        startup logic has been performed yet, initialization of dependencies
        is deferred until significantly later in the startup process.
        '''

        # Avoid circular import dependencies.
        from betse.lib import libs
        from betse.util.io.error import errfault
        from betse.util.io.log import logconfig, logs
        from betse.util.os import oses
        from betse.util.py import pys
        from betse.util.test import tests
        from betse.util.type.obj import objects

        # Enable Python's standard handler for segmentation faults *BEFORE*
        # performing any further logic, any of which could conceivably trigger
        # a segmentation fault and hence process termination.
        errfault.handle_faults()

        # Enable our default logging configuration for the current Python
        # process *BEFORE* performing any validation, thus logging any
        # exceptions raised by this validation.
        logconfig.init()

        # Log all prior behaviour. Attempting to do so *BEFORE* enabling our
        # default logging configuration above would silently fail, since the
        # standard "logging" API silently squelches debug messages by default.
        logs.log_debug('Application singleton "%s" established.',
            objects.get_class_name_unqualified((self)))
        logs.log_debug('Default segementation fault handler enabled.')
        logs.log_debug('Testing environment detected: %r', tests.is_testing())

        # Validate mandatory dependencies. Avoid initializing these
        # dependencies now (e.g., by calling init_libs()). Doing so requires
        # finalization of the logging configuration (e.g., by parsing CLI
        # options), which has yet to happen this early in the lifecycle.
        libs.die_unless_runtime_mandatory_all()

        # Validate the active Python interpreter and operating system *AFTER*
        # mandatory dependencies. While the former (mostly) comprises
        # unenforced recommendations, the latter comprises enforced
        # requirements and should thus be validated first.
        oses.init()
        pys.init()

    # ..................{ INITIALIZERS ~ libs               }..................
    @type_check
    def init_libs(
        self, matplotlib_backend_name: StrOrNoneTypes = None) -> None:
        '''
        Initialize all mandatory third-party dependencies of this application
        with sane defaults.

        Specifically, this method (in no particular order):

        * Reconfigures matplotlib with sane defaults specific to the current
          platform and set of all available third-party GUI frameworks.
        * Initializes exactly one available third-party YAML parsing framework
          (e.g., PyYaml, :mod:`ruamel.yaml`).
        * Initializes NumPy.
        * Initializes Pillow.

        Parameters
        ----------
        matplotlib_backend_name : StrOrNoneTypes
            Name of the matplotlib backend to explicitly enable. Defaults to
            ``None``, in which case this method implicitly enables the first
            importable backend known to be both usable and supported by this
            application (in descending order of preference).
        '''

        # Avoid circular import dependencies.
        from betse.lib.matplotlib.matplotlibs import mpl_config
        from betse.lib.numpy import numpys
        from betse.lib.pickle import pickles
        from betse.lib.pil import pils
        from betse.lib.yaml import yamls
        from betse.util.io.log import logs

        # Log this initialization. Since initializing heavyweight third-party
        # dependencies (especially matplotlib) consumes non-trivial time, this
        # message is intentionally exposed to all users by default.
        logs.log_info('Loading third-party dependencies...')

        # Initialize these dependencies in arbitrary order.
        mpl_config.init(backend_name=matplotlib_backend_name)
        numpys.init()
        pickles.init()
        pils.init()
        yamls.init()

        # Note this method to have been called *AFTER* successfully doing so.
        self._is_libs_initted = True


    def init_libs_if_needed(self, *args, **kwargs) -> None:
        '''
        Initialize all mandatory third-party dependencies of the current
        application with sane defaults if these dependencies have yet to be
        initialized (i.e., if the :meth:`init_libs` method has yet to be
        called) *or* silently reduce to a noop otherwise.

        Parameters
        ----------
        All passed parameters are passed as is to the :meth:`init_libs` method.
        '''

        from betse.util.io.log import logs

        # If the init_libs() method has already been called, log this fact
        # *WITHOUT* recalling that method, thus reducing to a noop.
        if self._is_libs_initted:
            logs.log_debug(
                'Ignoring request to reload third-party dependencies...')
        # Else, the init_libs() method has yet to be called. So, do so.
        else:
            self.init_libs(*args, **kwargs)

    # ..................{ PROPERTIES ~ bool                 }..................
    @property_cached
    def is_git_worktree(self) -> bool:
        '''
        ``True`` only if this application has a Git-based **working tree**
        (i.e., top-level directory containing this application's ``.git``
        subdirectory and ``setup.py`` install script), typically due to this
        application having been installed for developer usage.
        '''

        return self.git_worktree_dirname_or_none is not None

    # ..................{ PROPERTIES ~ dir                  }..................
    @property_cached
    def project_dirname(self) -> str:
        '''
        Absolute dirname of this application's **root project directory**
        (i.e., top-level directory containing this application's installable
        ``pyproject.toml`` file or ``setup.py`` script) if found *or* raise an
        exception otherwise (i.e., if this directory is *not* found).

        Equivalently, this is the same as:

        * The root directory archived by release tarballs for this application.
        * The Git-based working tree for this application (i.e., the top-level
          directory containing this application's ``.git`` subdirectory).

        Caveats
        ----------
        **This directory typically does not exist.** This directory is only
        required during installation by non-developers *or* during development
        by developers. Once this application has been installed in a standard
        (i.e., non-editable) fashion by non-developers, this directory is no
        longer required and hence should *not* be assumed to exist.

        Raises
        ----------
        BetseDirException
            If this directory does *not* exist.
        '''

        # Avoid circular import dependencies.
        from betse.util.path import dirs
        from betse.util.py.module import pypackage

        # Absolute dirname of the parent directory of our top-level package.
        package_dirname = pypackage.get_package_project_dirname(self.package)

        # If this directory is not found, fail; else, return this directory.
        return dirs.dir_or_die(package_dirname)


    @property_cached
    def package_dirname(self) -> str:
        '''
        Absolute dirname of this application's root package directory if found
        *or* raise an exception otherwise (i.e., if this directory is *not*
        found).

        This directory typically resides in the ``site-packages`` subdirectory
        of the system-wide standard lib for the active Python interpreter
        (e.g., ``/usr/lib64/python3.6/site-packages/betse``).

        Raises
        ----------
        BetseDirException
            If this directory does *not* exist.
        '''

        # Avoid circular import dependencies.
        from betse.util.path import dirs
        from betse.util.py.module import pymodule

        # Absolute dirname of the directory yielding our top-level package.
        package_dirname = pymodule.get_dirname(self.package)

        # If this directory is not found, fail; else, return this directory.
        return dirs.dir_or_die(package_dirname)


    @property_cached
    def dot_dirname(self) -> str:
        '''
        Absolute dirname of this application's **root dot directory** (i.e.,
        top-level directory containing this application's user-specific files
        and hence residing in the home directory of the current user), silently
        creating this directory if *not* already found.

        This directory contains user-specific files (e.g., logfiles, profile
        files) both read from and written to at application runtime. These are
        typically plaintext files consumable by external users and third-party
        utilities.

        Locations
        ----------
        Denote:

        * ``{package_name}`` the value of the :meth:`package_name` property for
          this application (e.g., ``betse`` for BETSE).
        * ``{username}`` the name of the current user (e.g., ``leycec``).

        Then the dirname returned by this property is:

        * Under Linux, ``~/.{package_name}/``. This property intentionally does
          *not* currently comply with the `XDG Base Directory Specification`_
          (e.g., ``~/.local/share/betse``), which the authors regard as
          `unhelpful if not blatantly harmful <xkcd Standards_>`__.
        * Under OS X, ``~/Library/Application Support/{package_name}``.
        * Under Windows,
          ``C:\\Documents and Settings\\{username}\\Application Data\\{package_name}``.

        .. _XDG Base Directory Specification:
            http://standards.freedesktop.org/basedir-spec/basedir-spec-latest.html
        .. _xkcd Standards:
            https://xkcd.com/927
        '''

        # Avoid circular import dependencies.
        from betse.util.os import oses
        from betse.util.os.shell import shellenv
        from betse.util.path import dirs, pathnames

        # Absolute dirname of this directory.
        dot_dirname = None

        # If the current platform is macOS, set the appropriate directory.
        if oses.is_macos():
            dot_dirname = pathnames.join(
                pathnames.get_home_dirname(),
                'Library',
                'Application Support',
                self.package_name,
            )
        # If the current platform is Windows, set the appropriate directory.
        elif oses.is_windows():
            dot_dirname = pathnames.join(
                shellenv.get_var('APPDATA'), self.package_name)
        # Else, assume the current platform to be POSIX-compatible.
        else:
            #FIXME: Explicitly assert POSIX compatibility here. To do so, we'll
            #want to define and call a new betse.util.os.oses.die_unless_posix()
            #function here.
            dot_dirname = pathnames.join(
                pathnames.get_home_dirname(), '.' + self.package_name)

        # Create this directory if not found.
        dirs.make_unless_dir(dot_dirname)

        # Return this dirname.
        return dot_dirname

    # ..................{ PROPERTIES ~ dir : data           }..................
    @property_cached
    def data_dirname(self) -> str:
        '''
        Absolute dirname of this application's top-level data directory if
        found *or* raise an exception otherwise (i.e., if this directory is
        *not* found).

        This directory typically contains application-internal resources (e.g.,
        media files) required at application runtime.

        Raises
        ----------
        BetseDirException
            If this directory does *not* exist.
        '''

        # Avoid circular import dependencies.
        from betse.util.app import apppath

        # Return the absolute dirname of this application-relative directory if
        # this directory exists *OR* raise an exception otherwise.
        return apppath.get_dirname(package=self.package, dirname='data')

    # ..................{ PROPERTIES ~ dir : git            }..................
    @property_cached
    def git_worktree_dirname(self) -> str:
        '''
        Absolute dirname of this application's Git-based **working tree**
        (i.e., top-level directory containing this application's ``.git``
        subdirectory and ``setup.py`` install script) if this application was
        installed in a developer manner *or* raise an exception otherwise
        (i.e., if this directory is *not* found).

        Raises
        ----------
        BetseDirException
            if this application was installed in a developer manner but this
            directory does *not* exist, a paradox that should never occur.
        BetseGitException
            if this application was *not* installed in a developer manner.
        '''

        # Avoid circular import dependencies.
        from betse.util.path import dirs

        # Absolute dirname of this application's Git-based working tree if
        # this application was installed for development or "None" otherwise.
        git_worktree_dirname = self.git_worktree_dirname_or_none

        # If this application is *NOT* under development, fail.
        if git_worktree_dirname is None:
            raise BetseGitException(
                'Package "{}" not under Git-based development.'.format(
                    self.package_name))

        # If this directory is not found, fail; else, return this directory.
        return dirs.dir_or_die(git_worktree_dirname)


    @property_cached
    def git_worktree_dirname_or_none(self) -> StrOrNoneTypes:
        '''
        Absolute dirname of this application's **Git-based working tree**
        (i.e., top-level directory containing this application's ``.git``
        subdirectory and ``setup.py`` install script) if this application was
        installed in a developer manner *or* ``None`` otherwise.

        Returns
        ----------
        StrOrNoneTypes
            Either:

            * The absolute dirname of this application's Git working tree if
              this application was installed in a developer manner: e.g., by

              * ``python3 setup.py develop``.
              * ``python3 setup.py symlink``.

            * ``None`` if this application was installed in a non-developer
              manner: e.g., by

              * ``pip3 install``.
              * ``python3 setup.py develop``.
        '''

        # Avoid circular import dependencies.
        from betse.util.path import gits

        # Behold! It is a one-liner.
        return gits.get_package_worktree_dirname_or_none(self.package)

    # ..................{ PROPERTIES ~ file                 }..................
    @property_cached
    def log_default_filename(self) -> str:
        '''
        Absolute filename of this application's **default logfile** (i.e.,
        plaintext user-specific file to which all messages, warnings, errors,
        and exceptions are logged by default).
        '''

        # Avoid circular import dependencies.
        from betse.util.path import pathnames

        # Return the absolute path of this file.
        return pathnames.join(self.dot_dirname, self.package_name + '.log')


    @property_cached
    def profile_default_filename(self) -> str:
        '''
        Absolute filename of this application's **default profile dumpfile**
        (i.e., user-specific binary file to which profiled statistics are saved
        by default).
        '''

        # Avoid circular import dependencies.
        from betse.util.path import pathnames

        # Return the absolute path of this file.
        return pathnames.join(self.dot_dirname, self.package_name + '.prof')

    # ..................{ PROPERTIES ~ module : root        }..................
    @property_cached
    def package(self) -> ModuleType:
        '''
        **Root package** (i.e., topmost package for this application, typically
        of the same name as this application and installed into a subdirectory
        of the same name in the ``site-packages`` directory specific to the
        active Python interpreter) for this application.
        '''

        # Avoid circular import dependencies.
        from betse.util.py.module import pypackage

        # Introspection for the glorious victory.
        return pypackage.get_object_type_package_root(obj=self)


    @property_cached
    def package_name(self) -> str:
        '''
        Name of this application's root package (e.g., ``betse`` for BETSE).
        '''

        # Avoid circular import dependencies.
        from betse.util.py.module import pymodule

        # By the power of Grayskull...
        return pymodule.get_name_qualified(module=self.package)

    # ..................{ PROPERTIES ~ module : test        }..................
    @property_cached
    def test_package_name(self) -> str:
        '''
        Name of the root package of this application's ancillary test suite
        (e.g., ``betse_test`` for BETSE).

        Caveats
        ----------
        **This package is typically not installed with this application,** as
        tests are useless (or at least incidental) for most end user purposes.
        Instead, this package is only distributed with tarballs archiving the
        contents of this application's repository at stable releases time. This
        package is *not* guaranteed to exist and, in fact, typically does not.
        '''

        # When our powers combine!
        return self.package_name + '_test'
