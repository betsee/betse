#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2017-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level **application metadata singleton** (i.e., application-wide object
synopsizing application metadata via read-only properties) hierarchy.
'''

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
# from betse.util.io.log import logs
from betse.util.type.decorator.deccls import abstractproperty
from betse.util.type.decorator.decmemo import property_cached
from betse.util.type.types import (
    type_check, ModuleType, StrOrNoneTypes)

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

# ....................{ SUPERCLASSES                      }....................
class MetaAppABC(object, metaclass=ABCMeta):
    '''
    Abstract base class of all **application metadata singleton** (i.e.,
    application-wide object synopsizing application metadata via read-only
    properties) subclasses.

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
    '''

    # ..................{ PROPERTIES ~ subclass             }..................
    # Abstract read-only properties required to be defined by subclasses.

    @abstractproperty
    def package(self) -> ModuleType:
        '''
        **Root Python package** (i.e., top-level package for this application,
        typically of the same name as this application and installed into a
        subdirectory of the same name in the ``site-packages`` directory
        specific to the active Python interpreter).
        '''

        pass

    # ..................{ PROPERTIES                        }..................
    # Concrete read-only properties.

    @property
    def package_name(self) -> str:
        '''
        Unqualified name of this application's root package (e.g., ``betse`` if
        this application is BETSE).
        '''

        return self.package.__name__

    # ..................{ PROPERTIES ~ dir : git            }..................
    @property_cached
    def package_dirname(self) -> str:
        '''
        Absolute pathname of this application's root package directory if found
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
        from betse.util.py import pymodule

        # Absolute pathname of the directory yielding the top-level "betse" package.
        package_dirname = pymodule.get_dirname(self.package)

        # If this directory is not found, fail.
        dirs.die_unless_dir(package_dirname)

        # Return this directory's pathname.
        return package_dirname

    # ..................{ PROPERTIES ~ dir : git            }..................
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
        from betse.util.path import dirs

        # Absolute path of this directory.
        data_dirname = self.get_pathname('data')

        # If this directory is not found, raise an exception.
        dirs.die_unless_dir(data_dirname)

        # Return the absolute path of this directory.
        return data_dirname

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

        # Absolute pathname of this application's Git-based working tree if
        # this application was installed for development or "None" otherwise.
        git_worktree_dirname = self.git_worktree_dirname_or_none

        # If this application is *NOT* under development, fail.
        if git_worktree_dirname is None:
            raise BetseGitException(
                'Package "{}" not under Git-based development.'.format(
                    self.package_name))

        # If this directory is not found, fail.
        dirs.die_unless_dir(git_worktree_dirname)

        # Return this directory's pathname.
        return git_worktree_dirname


    @property_cached
    def git_worktree_dirname_or_none(self) -> StrOrNoneTypes:
        '''
        Absolute dirname of this application's Git-based **working tree**
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

    # ..................{ PROPERTIES ~ private              }..................

    # ..................{ GETTERS ~ path                    }..................
    @type_check
    def get_pathname(self, pathname: str) -> str:
        '''
        Absolute path of the passed pathname relative to the absolute path of
        the root package for this application.

        Specifically, this method returns:

        * If this application is a PyInstaller-frozen executable binary, the
          concatenation of (in order):

          #. The absolute path of the temporary directory containing all
             application data resources extracted from this binary by this
             executable's bootloader, as specified by the PyInstaller-specific
             private attribute ``_MEIPASS`` injected into the canonical
             :mod:`sys` module by the PyInstaller bootloader embedded in this
             binary. "And it's turtles all the way down."
          #. The passed relative pathname.

        * If this application is a :mod:`setuptools`-installed script wrapper,
          the result of querying :mod:`setuptools` for the absolute path of the
          passed relative pathname. In this case, this path will have been
          preserved as is in the :mod:`setuptools`-installed copy of this
          application in the package tree for the active Python interpreter.
        * Else, the concatenation of (in order):

          #. The absolute path of the directory providing this root package.
          #. The passed relative pathname.

          In this case, this application is typically either a
          :mod:`setuptools`-symlinked script wrapper *or* was invoked via the
          secretive ``python3 -m {package.__name__}`` command.

        Parameters
        ----------
        pathname : str
            Relative pathname of the path to be canonicalized.

        Returns
        ----------
        Absolute path of this pathname relative to the absolute path of the
        root package for this application.
        '''

        # Avoid circular import dependencies.
        from betse.lib.setuptools import resources
        from betse.util.path import pathnames
        from betse.util.py import pyfreeze

        # If this application is frozen by PyInstaller, canonicalize this path
        # relative to the directory to which this application is unfrozen.
        if pyfreeze.is_frozen_pyinstaller():
            app_pathname = pathnames.join(
                pyfreeze.get_app_dirname_pyinstaller(), pathname)
        # Else if this application is a setuptools-installed script wrapper,
        # canonicalize this path by deferring to the setuptools resource API.
        elif resources.is_dir(
            module_name=self.package_name, dirname=pathname):
            app_pathname = resources.get_pathname(
                module_name=self.package_name, pathname=pathname)
        # Else, the current application is either a setuptools-symlinked script
        # wrapper *OR* was invoked via the secretive "python3 -m betse"
        # command. In either case, this directory's path is directly obtainable
        # relative to the absolute path of the passed package.
        else:
            app_pathname = pathnames.join(self.package_dirname, pathname)

        # Return this pathname.
        return app_pathname
