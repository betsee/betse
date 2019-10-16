#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level module and package facilities.

All functions defined by this submodule require that modules be passed as
fully-qualified names (i.e., ``.``-delimited).

All functions defined by this submodule require that modules be at least passed
as previously imported module objects; most also accept fully-qualified names
of modules where applicale.

See Also
----------
:mod:`betse.util.py.module.pymodname`
    Related submodule whose functions accept only fully-qualified names.
'''

# ....................{ IMPORTS                           }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To raise human-readable exceptions on missing mandatory
# dependencies, the top-level of this module may import *ONLY* from packages
# guaranteed to exist at installation time -- which typically means *ONLY*
# BETSE packages and stock Python packages.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

from betse.exceptions import BetseModuleException
from betse.util.io.log import logs
from betse.util.type import types
from betse.util.type.types import (
    type_check,
    ModuleType,
    ModuleOrStrTypes,
    SetType,
    StrOrNoneTypes,
)
from importlib.machinery import ExtensionFileLoader, EXTENSION_SUFFIXES

# ....................{ EXCEPTIONS                        }....................
@type_check
def die_if_builtin(module: ModuleOrStrTypes) -> None:
    '''
    Raise an exception if the passed module is a **builtin** (i.e., compiled
    into the executable file for the active Python interpreter and hence having
    no external file).

    Parameters
    ----------
    module : ModuleOrStrTypes
        Either:

        * The fully-qualified name of this module, in which case this function
          dynamically imports this module.
        * A previously imported module object.

    Raises
    ----------
    BetseModuleException
        If this module is a builtin.

    See Also
    ----------
    :func:`is_builtin`
        Further details.
    '''

    # If this module is a builtin, raise an exception.
    if is_builtin(module):
        raise BetseModuleException(
            'Module "{0}" is a builtin '
            '(i.e., "{0}.__file__" undefined).'.format(
                get_name_qualified(module)))

# ....................{ EXCEPTIONS ~ name                 }....................
@type_check
def die_unless_topmost(module: ModuleOrStrTypes) -> None:
    '''
    Raise an exception unless the passed module is **topmost** (i.e., a
    top-level module whose module name contains no ``.`` delimiters).

    Parameters
    ----------
    module : ModuleOrStrTypes
        Either:

        * The fully-qualified name of this module, in which case this function
          dynamically imports this module.
        * A previously imported module object.

    Raises
    ----------
    BetseModuleException
        If this module is a submodule rather than topmost.
    '''

    # If this module is *NOT* topmost, raise an exception.
    if not is_topmost(module):
        raise BetseModuleException(
            'Module "{}" not topmost '
            '(i.e., contains one or more "." delimiters).'.format(
                get_name_qualified(module)))

# ....................{ TESTERS                           }....................
@type_check
def is_module(obj: object) -> bool:
    '''
    ``True`` only if the passed object is an **in-memory module** (i.e.,
    instance of the :class:`ModuleType` class).
    '''

    # The best things in life are anguish-free.
    return isinstance(obj, ModuleType)


@type_check
def is_builtin(module: ModuleOrStrTypes) -> bool:
    '''
    ``True`` only if the passed module is a **builtin** (i.e., compiled into
    the executable file for the active Python interpreter and hence having no
    external file).

    Parameters
    ----------
    module : ModuleOrStrTypes
        Either:

        * The fully-qualified name of this module, in which case this function
          dynamically imports this module.
        * A previously imported module object.

    Returns
    ----------
    bool
        ``True`` only if this module is a builtin.
    '''

    # Return true only if this module fails to define the special "__file__"
    # attribute, required to be defined by all modules *EXCEPT* builtins.
    return not hasattr(module, '__file__')

# ....................{ TESTERS ~ c                       }....................
@type_check
def is_c(module: ModuleOrStrTypes) -> bool:
    '''
    ``True`` only if the passed module is **implemented in C** (i.e., as either
    an internally builtin module *or* external C extension).

    Parameters
    ----------
    module : ModuleOrStrTypes
        Either:

        * The fully-qualified name of this module, in which case this function
          dynamically imports this module.
        * A previously imported module object.

    Returns
    ----------
    bool
        ``True`` only if this module is implemented in C.

    See Also
    ----------
    :func:`is_builtin`
    :func:`is_c_extension`
        Further details.
    '''

    return is_builtin(module) or is_c_extension(module)


@type_check
def is_c_extension(module: ModuleOrStrTypes) -> bool:
    '''
    ``True`` only if the passed module is a C extension implemented as a
    dynamically linked shared library specific to the current platform.

    Parameters
    ----------
    module : ModuleOrStrTypes
        Either:

        * The fully-qualified name of this module, in which case this function
          dynamically imports this module.
        * A previously imported module object.

    Returns
    ----------
    bool
        ``True`` only if this module is a C extension.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import pathnames

    # Resolve this module's object.
    module = resolve_module(module)

    # If this module was loaded by a PEP 302-compliant C extension loader, this
    # module *MUST* be a C extension.
    if isinstance(getattr(module, '__loader__', None), ExtensionFileLoader):
        return True
    # Else if this module is a builtin, this module has no external file and
    # hence *CANNOT* be a C extension.
    elif is_builtin(module):
        return False

    # Else, fallback to filetype matching heuristics.
    #
    # Absolute path of the file defining this module.
    module_filename = get_filename(module)

    # "."-prefixed filetype of this path if any or "None" otherwise.
    module_filetype = pathnames.get_filetype_dotted_or_none(module_filename)
    # print('module_filetype: {}'.format(module_filetype))

    #FIXME: Mildly inefficient, as "EXTENSION_SUFFIXES" is a list rather than a
    #set. Since this list is small *AND* since this function is called
    #infrequently, this is currently ignorable. The trivial fix is to define a
    #new private "_EXTENSION_SUFFIXES = set(EXTENSION_SUFFIXES)" global above
    #and leverage that here instead.

    # This module is only a C extension if this path's filetype is that of a
    # C extension specific to the current platform.
    return module_filetype in EXTENSION_SUFFIXES

# ....................{ TESTERS : name                    }....................
@type_check
def is_topmost(module: ModuleOrStrTypes) -> bool:
    '''
    ``True`` only if the passed module is **topmost** (i.e., a top-level module
    whose module name contains no ``.`` delimiters).

    Parameters
    ----------
    module : ModuleOrStrTypes
        Either:

        * The fully-qualified name of this module, in which case this function
          dynamically imports this module.
        * A previously imported module object.

    Returns
    ----------
    bool
        ``True`` only if this a topmost module.
    '''

    # Fully-qualified name of this module.
    module_name = get_name_qualified(module)

    # Return true only if this name contains no "." delimiters.
    return '.' not in module_name

# ....................{ GETTERS ~ attr : global           }....................
#FIXME: Rename this getter to iter_global_names() for orthogonality with
#similar iterators in the "betse.util.type.obj.objiter" submodule.
@type_check
def get_global_names(module: ModuleOrStrTypes) -> SetType:
    '''
    Set of the names of all global variables defined by the passed module.

    This function returns the set of the names of all attributes defined by
    this module, excluding:

    * Special attributes reserved for use by Python (e.g., ``__file__``).
    * Callable attributes (e.g., functions, lambdas).

    Parameters
    ----------
    module : ModuleOrStrTypes
        Either:

        * The fully-qualified name of this module, in which case this function
          dynamically imports this module.
        * A previously imported module object.

    Returns
    ----------
    set
        Set of the names of all global variables defined by this module.
    '''

    # Avoid circular import dependencies.
    from betse.util.py import pyident

    # Resolve this module's object.
    module = resolve_module(module)

    # Return this set via a set comprehension.
    return {
        module_attr_name
        for module_attr_name in dir(module)
        if (
            pyident.is_special(module_attr_name) and
            callable(getattr(module, module_attr_name))
        )
    }

# ....................{ GETTERS ~ name                    }....................
@type_check
def get_name_qualified(module: ModuleOrStrTypes) -> str:
    '''
    **Fully-qualified name** (i.e., absolute canonical ``.``-delimited name) of
    the passed module.

    While trivial, this getter is defined for orthogonality with comparable
    non-trivial getters defined by this submodule.

    Parameters
    ----------
    module : ModuleOrStrTypes
        Either:

        * The fully-qualified name of this module, in which case this function
          dynamically imports this module.
        * A previously imported module object.

    Returns
    ----------
    str
        Fully-qualified name of this module.
    '''

    # Resolve this module's object.
    module = resolve_module(module)

    # Who let the one-liners out?
    return module.__name__

# ....................{ GETTERS ~ path                    }....................
@type_check
def get_filename(module: ModuleOrStrTypes) -> str:
    '''
    Absolute filename of the file providing the passed module or package if
    this is not a builtin module *or* raise an exception otherwise (i.e., if
    this is a builtin module).

    If the passed object signifies:

    * A package (e.g., directory), this is the absolute path of the file
      providing this package's ``__init__`` submodule.
    * A non-package (e.g., module, C extension), this is the absolute path of
      the file providing this non-package as is.

    Caveats
    ----------
    Since effectively *all* modules and packages (except in-memory builtin
    modules) are associated with a corresponding file, there intentionally
    exists no corresponding ``get_filename_or_none()`` function.

    Since *only* packages define the ``__path__`` attribute, there likewise
    exists no corresponding ``get_pathname()`` function. Only the ``__file__``
    attribute retrieved by this function is generally applicable to *all*
    modules and packages regardless of type.

    Parameters
    ----------
    module : str or ModuleType
        Either:

        * The fully-qualified name of this module, in which case this function
          dynamically imports this module.
        * A previously imported module object.

    Returns
    ----------
    str
        Absolute filename of this file providing this module or package.

    Raises
    ----------
    BetseModuleException
        If this is a **builtin module** (i.e., compiled into the executable
        file for the active Python interpreter and hence having no external
        file).
    '''

    # Resolve this module's object.
    module = resolve_module(module)

    # If this module is a builtin (and hence does *NOT* define the standard
    # "__file__" attribute), raise an exception.
    die_if_builtin(module)
    # Else, this module is *NOT* a builtin and hence defines this attribute.

    # Return this attribute's value.
    return module.__file__

# ....................{ GETTERS ~ path : dir              }....................
@type_check
def get_dirname(module: ModuleOrStrTypes) -> str:
    '''
    Absolute dirname of the directory providing the passed module or package.

    Specifically, if the passed module is a:

    * Module, this is the directory containing this module.
    * Non-namespace package, this is the directory containing this package's
      top-level ``__init__`` submodule.

    Caveats
    ----------
    **This dirname is not guaranteed to be canonical.** Since this function
    does *not* resolve symbolic links in this dirname, one or more directory
    components of this dirname may be symbolic links. Consider calling the
    :func:`get_dirname_canonical` function instead if this is undesirable.

    Parameters
    ----------
    module : ModuleOrStrTypes
        Either:

        * The fully-qualified name of this module, in which case this function
          dynamically imports this module as a non-optional side effect.
        * A previously imported module object.

    Returns
    ----------
    str
        Absolute dirname of the directory providing this module or package.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import pathnames

    # Return the dirname of the directory containing the Python script
    # implementing this module. Note that get_filename() returns the absolute
    # filename of the "__init__.py" file if this module is a package,
    # guaranteeing the correct behaviour regardless of whether this is
    # specifically a module or package.
    return pathnames.get_dirname(get_filename(module))


@type_check
def get_dirname_canonical(module: ModuleOrStrTypes) -> str:
    '''
    **Absolute canonical dirname** (i.e., absolute dirname after resolving
    symbolic links) of the directory providing the passed module or package.

    See Also
    ----------
    :func:`get_dirname`
        Further details.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import pathnames

    # Return this dirname canonicalized.
    return pathnames.canonicalize(get_dirname(module))

# ....................{ GETTERS ~ version                 }....................
@type_check
def get_version(module: ModuleOrStrTypes) -> str:
    '''
    Version specifier of the passed module if that module defines a version
    specifier *or* raise an exception otherwise.

    See Also
    ----------
    :func:`get_version_or_none`
        Further details on the passed parameter.
    '''

    # Module version if any or "None" otherwise.
    module_version = get_version_or_none(module)

    # If this version does *NOT* exist, raise an exception.
    if module_version is None:
        raise BetseModuleException(
            'Module "{}" version not found.'.format(module))

    # Return this version.
    return module_version


@type_check
def get_version_or_none(module: ModuleOrStrTypes) -> StrOrNoneTypes:
    '''
    Version specifier of the passed module if that module defines a version
    specifier *or* ``None`` otherwise.

    Parameters
    ----------
    module : ModuleOrStrTypes
        Either:
        * The fully-qualified name of this module, in which case this function
          dynamically imports this module.
        * A previously imported module object.

    Returns
    ----------
    StrOrNoneTypes
        This module's version specifier if any *or* ``None`` otherwise.
    '''

    # Avoid circular import dependencies.
    from betse.util.py.module.pymodname import MODULE_NAME_TO_VERSION_ATTR_NAME

    # Resolve this module's object.
    module = resolve_module(module)

    # Name of the version specifier attribute defined by that module. For sane
    # modules, this is "__version__". Insane modules, however, exist.
    module_version_attr_name = MODULE_NAME_TO_VERSION_ATTR_NAME[
        module.__name__]

    # This attribute defined by this module if any *OR* "None" otherwise.
    module_version = getattr(module, module_version_attr_name, None)

    # If this version is undefined, log a non-fatal warning.
    if module_version is None:
        logs.log_warning('Module "%s" version not found.', module.__name__)

    # Return this version.
    return module_version

# ....................{ RESOLVERS                         }....................
@type_check
def resolve_module(module: ModuleOrStrTypes) -> ModuleType:
    '''
    Dynamically import and return the module with the passed name if the passed
    parameter is a string *or* return the passed module as is otherwise.

    Motivation
    ----------
    This utility function is principally intended to simplify the
    implementation of public functions defined by this submodule and the
    :mod:`betse.util.py.module.pypackage` submodule.

    Parameters
    ----------
    module : ModuleOrStrTypes
        Either:

        * The fully-qualified name of this module, in which case this function
          dynamically imports this module.
        * A previously imported module object.

    Returns
    ----------
    ModuleType
        Module object resolved from the passed parameter.
    '''

    # Avoid circular import dependencies.
    from betse.util.py.module import pymodname

    # Return either...
    return (
        # A module object dynamically imported as this module name...
        pymodname.import_module(module)
        # If a module name rather than object was passed; else...
        if types.is_str(module) else
        # This module object.
        module
    )
