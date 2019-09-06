#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level package-specific facilities.

All functions defined by this submodule accept at least a previously imported
package object; most also accept the fully-qualified name of a package.

See Also
----------
:mod:`betse.util.py.module.pyname`
    Related submodule whose functions accept only fully-qualified names.
'''

# ....................{ IMPORTS                           }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To raise human-readable exceptions on missing mandatory
# dependencies, the top-level of this module may import *ONLY* from packages
# guaranteed to exist at installation time -- which typically means *ONLY*
# BETSE packages and stock Python packages.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

from betse.exceptions import BetsePackageException
# from betse.util.io.log import logs
from betse.util.type.types import (
    type_check, ModuleType, StrOrNoneTypes)

# ....................{ GETTERS ~ object                  }....................
@type_check
def get_object_type_package_root(obj: object) -> ModuleType:
    '''
    **Root package** (i.e., topmost package whose package name contains no
    ``.`` delimiters) transitively defining the class of the passed object.

    Parameters
    ----------
    obj : object
        Object to retrieve the root package name of.

    Returns
    ----------
    ModuleType
        Root package transitively defining this object's class.
    '''

    # Avoid circular import dependencies.
    from betse.util.py.module import pymodname

    # Fully-qualified name of the root package of this object's class.
    class_package_root_name = get_object_type_package_root_name(obj)

    # Return the previously imported instance of this root package.
    return pymodname.import_module(module_name=class_package_root_name)


@type_check
def get_object_type_package_root_name(obj: object) -> str:
    '''
    Name of the **root package** (i.e., topmost package whose package name
    contains no ``.`` delimiters) transitively defining the class of the passed
    object if any *or* raise an exception otherwise (i.e., if that class is
    defined by no module, as is the case for standard C-based classes).

    Design
    ----------
    The name of this function is intentionally suffixed by neither
    ``_qualified`` nor ``_unqualified``. Since root package names contain no
    ``.`` delimiters (e.g., :mod:`betse`), the fully-qualified and unqualified
    names for any root package are necessarily identical.

    Parameters
    ----------
    obj : object
        Object to retrieve the root package name of.

    Returns
    ----------
    str
        Name of the root package transitively defining this object's class.

    Raises
    ----------
    BetseTypeException
        If this object's class is *not* defined by a module (e.g., :mod:`str`).
        See the :func:`betse.util.type.cls.classes.get_module_name_qualified`
        function for further details.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.cls import classes
    from betse.util.type.obj import objects
    from betse.util.type.text.string import strs

    # Class of this object.
    cls = objects.get_class(obj)

    # Fully-qualified name of the module defining this class if any *OR* raise
    # the "BetseTypeException" otherwise.
    module_name = classes.get_module_name_qualified(cls)

    # Name of the root package defining this module. Note that:
    #
    # * If this module name is itself this package name, this module name
    #   contains no "." delimiter, in which case this call correctly returns
    #   this module name as is.
    # * Else, this module name contains at least one "." delimiter, in which
    #   case this call correctly returns the prefix of this module name
    #   preceding the first such ".".
    package_root_name = strs.get_prefix_preceding_char_or_text(
        text=module_name, char='.')

    # Return this name.
    return package_root_name

# ....................{ GETTERS ~ path                    }....................
@type_check
def get_package_project_dirname(package: ModuleType) -> str:
    '''
    **Absolute canonical dirname** (i.e., absolute dirname after resolving
    symbolic links) of the **root project directory** (i.e., top-level
    directory containing an installable ``pyproject.toml`` file or ``setup.py``
    script) governing the passed top-level Python package if found *or* raise
    an exception otherwise (i.e., if this directory is *not* found).

    Raises
    ----------
    BetsePackageException
        If no such directory exists for this package.

    See Also
    ----------
    :func:`get_package_project_dirname_or_none`
        Further details.
    '''

    # Avoid circular import dependencies.
    from betse.util.py.module import pymodule

    # Absolute canonical dirname of this package's project directory if found
    # *OR* "None" otherwise.
    package_project_dirname = get_package_project_dirname_or_none(package)

    # If this directory does *NOT* exist, raise an exception.
    if package_project_dirname is None:
        raise BetsePackageException(
            'Package "{}" project directory not found.'.format(
                pymodule.get_name_qualified(package)))
    # Else, this directory exists.

    # Return this dirname.
    return package_project_dirname


@type_check
def get_package_project_dirname_or_none(
    package: ModuleType) -> StrOrNoneTypes:
    '''
    **Absolute canonical dirname** (i.e., absolute dirname after resolving
    symbolic links) of the **root project directory** (i.e., top-level
    directory containing an installable ``pyproject.toml`` file or ``setup.py``
    script) governing the passed top-level Python package if found *or* raise
    an exception otherwise (i.e., if this directory is *not* found).

    Equivalently, this is the same as both:

    * The root directory archived by release tarballs for this application.
    * The Git-based working tree for this application (i.e., the top-level
      directory containing this application's ``.git`` subdirectory).

    Specifically, this function returns non-``None`` only if the parent
    directory of the passed package's directory contains either:

    * A ``pyproject.toml`` file, implying this package to be satisfy `PEP 518`_
      and hence be installable by at least poetry_.
    * A ``setup.py`` script, implying this package to be installable by either
      the standard :mod:`distutils` API or third-party :mod:`setuptools` API.

    :: _PEP 518:
       https://snarky.ca/clarifying-pep-518/
    :: _poetry:
       https://github.com/sdispater/poetry

    Caveats
    ----------
    **This directory typically does not exist.** This directory is only
    required during installation by non-developers *or* during development
    by developers. Once this application has been installed in a standard
    (i.e., non-editable) fashion by non-developers, this directory is no
    longer required and hence should *not* be assumed to exist.

    Parameters
    ----------
    package : ModuleType
        Top-level Python package to be queried.

    Returns
    ----------
    StrOrNoneTypes
        Either:

        * If the parent directory of the passed package's directory contains
          either a ``pyproject.toml`` or ``setup.py`` script, the absolute
          canonical dirname of that parent directory.
        * Else, ``None``.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import files, pathnames
    from betse.util.py.module import pymodule

    # Absolute canonical dirname of the directory providing this package,
    # canonicalized into a directory rather than symbolic link to increase the
    # likelihood of obtaining the actual parent directory of this package.
    package_dirname = pymodule.get_dirname_canonical(package)

    # Absolute dirname of the parent directory of this directory.
    package_project_dirname = pathnames.get_dirname(package_dirname)

    # Absolute filenames of the "pyproject.toml" and "setup.py" files possibly
    # provided by this parent directory.
    pyproject_filename = pathnames.join(
        package_project_dirname, 'pyproject.toml')
    setup_filename = pathnames.join(package_project_dirname, 'setup.py')

    # Return this parent directory's absolute dirname if either of these
    # requisite installation-time files exist *OR* "None" otherwise.
    return (
        package_project_dirname
        if (
            files.is_file(pyproject_filename) or
            files.is_file(setup_filename)
        ) else
        None
    )
