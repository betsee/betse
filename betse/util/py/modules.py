#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level module facilities.
'''

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To raise human-readable exceptions on missing mandatory dependencies,
# the top-level of this module may import *ONLY* from packages guaranteed to
# exist at installation time -- which typically means *ONLY* BETSE packages and
# stock Python packages.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

import collections, importlib, sys
from betse.exceptions import BetseExceptionModule
from betse.util.type import types

# ....................{ GLOBALS ~ dict                     }....................
SETUPTOOLS_PROJECT_TO_MODULE_NAME = {
    'Matplotlib': 'matplotlib',
    'Numpy': 'numpy',
    'Pillow': 'PIL',
    'PyYAML': 'yaml',
    'SciPy': 'numpy',
    'setuptools': 'setuptools',
    'six': 'six',
    'yamale': 'yamale',
}
'''
Dictionary mapping each relevant `setuptools`-specific project name (e.g.,
`PyYAML`) to the fully-qualified name of the corresponding top-level module or
package providing that project (e.g., `yaml`).

This dictionary is principally used to validate the importability of runtime
dependencies at startup. For consistency, the size of this dictionary should be
greater than or equal to the size of the `betse.metadata.DEPENDENCIES_RUNTIME`
list.
'''


MODULE_TO_VERSION_ATTR_NAME = collections.defaultdict(
    # Default attribute name to be returned for all unmapped modules.
    lambda: '__version__',

    # Module-specific attribute names.
    PIL = 'PILLOW_VERSION',
)
'''
Dictionary mapping the fully-qualified name of a module or package (e.g., `PIL`)
to the name of the attribute (e.g., `PILLOW_VERSION`) declared by that
module or package, providing that module or package's version specifier.

All modules and packages unmapped by this dictionary default to the canonical
`__version__` attribute name.
'''

# ....................{ EXCEPTIONS                         }....................
def die_unless_module(
    module_name: str, exception_message: str = None) -> None:
    '''
    Raise an exception with the passed message (defaulting to a message
    synthesized from the passed module name) if the module with the passed name
    is _not_ importable by the active Python interpreter.

    If this module is a **submodule** (i.e., if this module's name contains one
    or more `.` characters), all transitive parent packages of this module will
    be iteratively imported as an unavoidable side effect of this function call.
    '''

    # If this module is unimportable, raise an exception.
    if not is_module(module_name):
        # If no exception message was passed, synthesize one from this name.
        if not exception_message:
            exception_message = 'Module "{}" not found.'.format(module_name)
        assert types.is_str(exception_message), (
            types.assert_not_str(exception_message))

        # Raise this exception.
        raise BetseExceptionModule(exception_message)

# ....................{ TESTERS                            }....................
def is_module(module_name: str) -> bool:
    '''
    `True` only if the module with the passed fully-qualified name is importable
    under the active Python interpreter.

    If this module is a **submodule** (i.e., contains a `.` character), all
    parent modules of this module will be imported as a side effect of this
    function call. Likewise, if this module is _not_ importable via standard
    mechanisms (e.g., the OS X-specific `PyObjCTools` package), the module
    itself may also be imported as a side effect.
    '''
    assert types.is_str_nonempty(module_name), (
        types.assert_not_str_nonempty(module_name, 'Module name'))

    # Depending on context, this function behaves in one of three distinct ways:
    #
    # * If this module's name is a key in the canonical dictionary "sys.modules"
    #   and has thus already been imported at least once under the current
    #   process, then...
    #   * If the "sys.modules[module_name].__spec__" attribute is set to a non-
    #     None value, such value is returned.
    #   * Else, the "ValueError" exception is raised.
    # * Else if this module is loadable by iteratively querying all module
    #   loaders in "sys.meta_path" (the canonical list of such loaders), a new
    #   spec is created describing this module and returned.
    # * Else, None is returned.
    #
    # Since this function only returns a single boolean, such return values and
    # exceptions are converted to simple boolean values.
    try:
        return importlib.util.find_spec(module_name) is not None
    # Unfortunately, this exception does *NOT* necessarily imply this module to
    # not exist. This module may exist even if this exception is thrown,
    # particularly for modules defined dynamically at runtime rather than loaded
    # from an external file. Further inspection is warranted.
    #
    # Unfortunately, this exception does imply that conventional alternatives to
    # the prior function call (e.g., testing tuples generated by
    # pkgutil.iter_modules()) will also fail to find this module. As a fallback,
    # attempt to manually import this module. Since doing so implicitly imports
    # the "__init__.py" files of all parent packages of this module and hence
    # may have unhelpful side effects, we do so only if the prior call failed.
    except ValueError:
        try:
            importlib.import_module(module_name)
            return True
        except ImportError:
            return False


def is_imported(*module_names) -> bool:
    '''
    `True` only if all modules with the passed fully-qualified names have
    already been imported under the current Python process.
    '''

    for module_name in module_names:
        assert types.is_str_nonempty(module_name), (
            types.assert_not_str_nonempty(module_name, 'Module name'))

        if module_name not in sys.modules:
            return False

    return True

# ....................{ GETTERS                            }....................
def get_dirname(mod) -> str:
    '''
    Get the absolute path of the directory containing the file from which the
    passed module was previously imported.
    '''
    assert hasattr(mod, '__file__'), '"{}" not a module.'.format(mod)

    # Avoid circular import dependencies.
    from betse.util.path import paths

    # Get such dirname.
    return paths.get_dirname(mod.__file__)

# ....................{ GETTERS ~ version                  }....................
def get_version(mod) -> str:
    '''
    Get the version specifier of the passed module.

    If that module provides no version specifier, an exception is raised.

    See Also
    ----------
    `get_version_or_none`
        For further details on the passed parameter.
    '''
    module_version = get_version_or_none(mod)

    # If such version does *NOT* exist, raise an exception.
    if module_version is None:
        raise BetseExceptionModule(
            'Module "%s" version not found.', str(mod))

    return module_version


def get_version_or_none(mod) -> str:
    '''
    Get the version specifier of the passed module if that module provides a
    version specifier _or_ `None` otherwise.

    For convenience, the passed module may be either:

    * The fully-qualified name of such module, in which case such module will be
      dynamically imported.
    * A previously imported module instance.
    '''

    # If a module name was passed, dynamically import this module.
    if types.is_str(mod):
        assert types.is_str_nonempty(mod), (
            types.assert_not_str_nonempty(mod, 'Module name'))
        mod = import_module(mod)

    # Name of the version specifier attribute defined by that module. For sane
    # modules, this is "__version__". Insane modules, however, exist.
    version_attr_name = MODULE_TO_VERSION_ATTR_NAME[mod.__name__]

    # Get such attribute from such module if defined or None otherwise.
    return getattr(mod, version_attr_name, None)

# ....................{ IMPORTERS                          }....................
def import_module(module_name: str, exception_message: str = None) -> type(sys):
    '''
    Dynamically import and return the module, package, or C extension with the
    passed fully-qualified name.

    If this module is unimportable, an exception with the passed message is
    raised.
    '''
    assert types.is_str_nonempty(module_name), (
        types.assert_not_str_nonempty(module_name, 'Module name'))

    # If this module is unimportable, raise an exception.
    die_unless_module(module_name, exception_message)

    # Else, import and return this module.
    return importlib.import_module(module_name)
