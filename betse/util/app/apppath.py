#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **application pathname** (i.e., pathnames relative to the current
application, regardless of the package manager used to install that application
on the local filesystem) hierarchy.
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

# from betse.util.io.log import logs
from betse.util.type.types import type_check, ModuleOrStrTypes

# ....................{ GETTERS                           }....................
@type_check
def get_pathname(package: ModuleOrStrTypes, pathname: str) -> str:
    '''
    Absolute pathname canonicalized from the passed relative pathname in a
    portable manner guaranteed to be relative to the absolute dirname of the
    passed application package if a path with this canonical pathname exists
    *or* raise an exception otherwise (i.e., if no such path exists).

    Specifically, this method returns:

    * If this application is a PyInstaller-frozen executable binary, the
      concatenation of (in order):

      #. The absolute path of the temporary directory containing all
         application data resources extracted from this binary by this
         executable's bootloader as specified by the PyInstaller-specific
         private attribute ``_MEIPASS`` injected into the canonical :mod:`sys`
         module by the PyInstaller bootloader embedded in this binary. "And
         it's turtles all the way down."
      #. The passed relative pathname.

    * If this application is a :mod:`setuptools`-installed script wrapper, the
      result of querying :mod:`setuptools` for the absolute path of the passed
      relative pathname. In this case, this path will have been preserved as is
      in the :mod:`setuptools`-installed copy of this application in the
      package tree for the active Python interpreter.
    * Else, the concatenation of (in order):

      #. The absolute path of the directory providing this root package.
      #. The passed relative pathname.

      In this case, this application is typically either a
      :mod:`setuptools`-symlinked script wrapper *or* was invoked via the
      secretive ``python3 -m {package.__name__}`` command.

    Parameters
    ----------
    package : ModuleOrStrTypes
        Topmost package of the application to canonicalize this pathname for,
        defined as either:

        * The fully-qualified name of this package, in which case this function
          dynamically imports this package.
        * A previously imported package object.
    pathname : str
        Relative pathname of the path to be canonicalized.

    Returns
    ----------
    Absolute pathname of this path relative to the absolute pathname of this
    application package.

    Raises
    ----------
    BetseModuleException
        If this package is a subpackage rather than topmost.
    BetsePathException
        If no path with the absolute pathname to be returned exists.
    BetsePathnameException
        If this pathname is absolute rather than relative.
    '''

    # Avoid circular import dependencies.
    from betse.lib.setuptools import supresource
    from betse.util.path import pathnames, paths
    from betse.util.py import pyfreeze
    from betse.util.py.module import pymodule

    # If this package is *NOT* topmost, raise an exception.
    pymodule.die_unless_topmost(package)

    # If this pathname is absolute rather than relative, raise an exception.
    pathnames.die_if_absolute(pathname)

    # Name of this package.
    package_name = pymodule.get_name_qualified(package)

    # If this application is frozen by PyInstaller, canonicalize this path
    # relative to the directory to which this application is unfrozen.
    if pyfreeze.is_frozen_pyinstaller():
        # Absolute dirname of the directory containing this frozen application.
        app_frozen_dirname = pyfreeze.get_app_dirname_pyinstaller()

        #FIXME: This requires generalization to the passed package, if in fact
        #this can be generalized. If this cannot be generalized, then some
        #other means will be needed to achieve the same or a similar effect.
        app_pathname = pathnames.join(app_frozen_dirname, pathname)
    # Else if this application is a setuptools-installed script wrapper,
    # canonicalize this path by deferring to the setuptools resource API.
    elif supresource.is_dir(
        module_name=package_name, dirname=pathname):
        app_pathname = supresource.get_pathname(
            module_name=package_name, pathname=pathname)
    # Else, the current application is either a setuptools-symlinked script
    # wrapper *OR* was invoked via the secretive "python3 -m betse"
    # command. In either case, this directory's path is directly obtainable
    # relative to the absolute path of the passed package.
    else:
        # Absolute canonical dirname of the directory defining this package.
        package_dirname = pymodule.get_dirname_canonical(package)

        # Canonicalize this path relative to this directory.
        app_pathname = pathnames.join(package_dirname, pathname)

    # If this path is not found, fail; else, return this path.
    return paths.path_or_die(app_pathname)


@type_check
def get_dirname(package: ModuleOrStrTypes, dirname: str) -> str:
    '''
    Absolute dirname canonicalized from the passed relative dirname in a
    portable manner guaranteed to be relative to the absolute dirname of the
    passed application package if a directory with this canonical dirname
    exists *or* raise an exception otherwise (i.e., if no such directory
    exists).

    Raises
    ----------
    BetseDirException
        If no directory with the absolute dirname to be returned exists.

    See Also
    ----------
    :func:`get_pathname`
        Further details.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import dirs

    # Absolute dirname of this directory.
    app_dirname = get_pathname(package=package, pathname=dirname)

    # If this directory is not found, fail; else, return this directory.
    return dirs.dir_or_die(app_dirname)
