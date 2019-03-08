#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level module and package name-based facilities.

All functions defined by this submodule accept only fully-qualified module
names.

See Also
----------
:mod:`betse.util.py.module.pyname`
    Related submodule whose functions accept imported module objects.
'''

# ....................{ IMPORTS                           }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To raise human-readable exceptions on missing mandatory dependencies,
# the top-level of this module may import *ONLY* from packages guaranteed to
# exist at installation time -- which typically means *ONLY* BETSE packages and
# stock Python packages.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

import importlib, sys
from betse.util.io.log import logs
from betse.util.type.iterable.mapping.mapcls import DefaultDict
from betse.util.type.types import type_check, ModuleType, StrOrNoneTypes
from collections import defaultdict
from importlib import util as importlib_util

# ....................{ GLOBALS ~ dict                    }....................
DISTUTILS_PROJECT_NAME_TO_MODULE_NAME = DefaultDict(
    missing_key_value=lambda self, missing_key: missing_key,
    initial_mapping={
        'Numpy': 'numpy',
        'Pillow': 'PIL',
        'PyYAML': 'yaml',
        'SciPy': 'scipy',
    }
)
'''
Dictionary mapping a relevant :mod:`distutils`-specific project name (e.g.,
``PyYAML``) to the fully-qualified name of the corresponding module or package
providing that project (e.g., ``yaml``).

All modules and packages explicitly unmapped by this dictionary default to the
identity mapping -- that is, mapping the fully-qualified names of those modules
and packages to themselves (e.g., from ``pytest`` to ``pytest``).
'''


MODULE_NAME_TO_VERSION_ATTR_NAME = defaultdict(
    # Default attribute name to be returned for all unmapped modules.
    lambda: '__version__',

    # Module-specific attribute names.
    PIL='PILLOW_VERSION',
)
'''
Dictionary mapping the fully-qualified name of a module or package (e.g.,
:mod:``PIL``) to the name of the attribute (e.g., ``PILLOW_VERSION``) declared
by that module or package, providing that module or package's version
specifier.

All modules and packages explicitly unmapped by this dictionary default to the
canonical ``__version__`` attribute name.
'''

# ....................{ EXCEPTIONS                        }....................
@type_check
def die_unless_module(
    module_name: str, exception_message: StrOrNoneTypes = None) -> None:
    '''
    Raise an exception with the passed message (defaulting to a message
    synthesized from the passed module name) if the module with the passed name
    is *not* importable by the active Python interpreter.

    Raises
    ----------
    ImportError
        If this module is unimportable. To permit callers to transparently
        handle importation errors in the standard way, this standard exception
        rather than an application-specific exception (e.g.,
        :class:`betse.exceptions.BetseModuleException`) is raised.

    See Also
    ----------
    :func:`is_module`
        Further details.
    '''

    # If this module is unimportable, raise an exception.
    if not is_module(module_name):
        # If no exception message was passed, synthesize one from this name.
        if not exception_message:
            exception_message = 'Module "{}" not found.'.format(module_name)

        # Raise this exception. (See the docstring for further details.)
        raise ImportError(exception_message)

# ....................{ TESTERS                           }....................
@type_check
def is_module(module_name: str) -> bool:
    '''
    ``True`` only if the module with the passed fully-qualified name is
    importable under the active Python interpreter.

    If this module is a **submodule** (i.e., contains a ``.`` character), all
    parent modules of this module will be imported as a side effect of this
    function call. Likewise, if this module is *not* importable via standard
    mechanisms (e.g., the OS X-specific :mod:`PyObjCTools` package), this
    module itself may also be imported as a side effect.
    '''

    # Depending on context, this function behaves in one of three distinct
    # ways:
    #
    # * If this module's name is a key in the canonical dictionary
    #   "sys.modules" and has thus already been imported at least once under
    #   the active Python process, then...
    #   * If the "sys.modules[module_name].__spec__" attribute is set to a
    #     non-None value, that value is returned.
    #   * Else, the "ValueError" exception is raised.
    # * Else if this module is loadable by iteratively querying all module
    #   loaders in "sys.meta_path" (the canonical list of such loaders), a new
    #   spec is created describing this module and returned.
    # * Else this module is unloadable. In this case:
    #   * If this module is a submodule (i.e., this module's name contains a
    #     "." delimiter) and any parent module of this submodule is unloadable,
    #     the "ImportError" exception is raised.
    #   * Else, None is returned.
    #
    # Since this function only returns a single boolean, these return values
    # and exceptions are converted to simple boolean values.
    try:
        return importlib_util.find_spec(module_name) is not None
    # If this module is a submodule (i.e., this module's name contains a "."
    # delimiter) and any parent module of this submodule is unloadable, this
    # submodule itself is unloadable.
    except ImportError:
        return False
    # If this module appears to have been imported at least once under the
    # active Python process but has no "__spec__" attribute, inspect deeper.
    # This exception does *NOT* necessarily imply this module to not exist.
    # This module may exist even if this exception is thrown -- namely for
    # modules dynamically defined at runtime rather than loaded from external
    # files.
    #
    # Unfortunately, this exception does imply that conventional alternatives
    # to the prior function call (e.g., testing tuples generated by
    # pkgutil.iter_modules()) will also fail to find this module. As a
    # fallback, attempt to manually import this module. Since doing so
    # implicitly imports the "__init__.py" files of all parent packages of this
    # module and hence may have unhelpful side effects, we do so only if the
    # prior call failed.
    except ValueError:
        try:
            importlib.import_module(module_name)
            return True
        except ImportError:
            return False


@type_check
def is_imported(*module_names: str) -> bool:
    '''
    ``True`` only if all modules with the passed fully-qualified names have
    already been imported in the active Python process.
    '''

    # all(). It is awesome.
    return all(module_name in sys.modules for module_name in module_names)

# ....................{ IMPORTERS                         }....................
@type_check
def import_module(
    module_name: str, exception_message: StrOrNoneTypes = None) -> ModuleType:
    '''
    Dynamically import and return the module, package, or C extension with the
    passed fully-qualified name if importable *or* raise an exception with the
    passed message otherwise.
    '''

    # If this module is unimportable, raise an exception.
    die_unless_module(module_name, exception_message)

    # Else, import and return this module.
    return importlib.import_module(module_name)


@type_check
def unimport_module_if_imported(*module_names: str) -> None:
    '''
    Dynamically unimport each of the modules, packages, or C extensions with
    the passed fully-qualified names that have been previously imported under
    the active Python interpreter.

    Specifically, for each passed name:

    * If this name is that of a previously imported module, this module is
      removed from the canonical :attr:`sys.modules` cache of all previously
      imported modules. On the next attempt to import this module, Python will
      re-import and re-cache that module to this cache.
    * Else, this module is silently ignored.
    '''

    # For each of the passed module names...
    for module_name in module_names:
        # If this module has already been imported, unimport this module.
        if is_imported(module_name):
            logs.log_debug('Unimporting module "{}"...'.format(module_name))
            del sys.modules[module_name]
        # Else, this module has *NOT* yet been imported. Ignore this module.
