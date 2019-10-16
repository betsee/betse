#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level module and package name-based facilities.

All functions defined by this submodule require that modules be passed as
fully-qualified names (i.e., ``.``-delimited).

See Also
----------
:mod:`betse.util.py.module.pymodule`
    Related submodule whose functions accept already imported module objects.
'''

# ....................{ IMPORTS                           }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To raise human-readable exceptions on missing mandatory dependencies,
# the top-level of this module may import *ONLY* from packages guaranteed to
# exist at installation time -- which typically means *ONLY* BETSE packages and
# stock Python packages.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

import importlib, sys
from betse.exceptions import BetseModuleException
from betse.util.io.log import logs
from betse.util.type.iterable.mapping.mapcls import DefaultDict
from betse.util.type.types import (
    type_check, MappingOrNoneTypes, ModuleType, StrOrNoneTypes)
from collections import defaultdict
from importlib import util as importlib_util
from importlib.machinery import ModuleSpec

# ....................{ GLOBALS                           }....................
DISTUTILS_PROJECT_NAME_TO_MODULE_NAME = DefaultDict(
    missing_key_value=lambda self, missing_key: missing_key,
    initial_mapping={
        'Numpy':  'numpy',
        'Pillow': 'PIL',
        'PyYAML': 'yaml',
        'SciPy':  'scipy',
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
def die_if_module(module_name: str) -> None:
    '''
    Raise an exception if the module with the passed name **exists** (i.e., is
    importable by the active Python interpreter).

    Parameters
    ----------
    module_name : str
        Fully-qualified name of the module to be validated.

    Raises
    ----------
    BetseModuleException
        If this module exists.

    See Also
    ----------
    :func:`is_module`
        Further details.
    '''

    # If this module is importable, raise an exception.
    if is_module(module_name):
        raise BetseModuleException(
            'Module "{}" already exists.'.format(module_name))


@type_check
def die_unless_module(
    module_name: str, exception_message: StrOrNoneTypes = None) -> None:
    '''
    Raise an exception with the passed message (defaulting to a message
    synthesized from the passed module name) if the module with the passed name
    does *not* **exist** (i.e., is unimportable by the active Python
    interpreter).

    Parameters
    ----------
    module_name : str
        Fully-qualified name of the module to be validated.
    exception_message : StrOrNoneTypes
        Message of the exception to be raised if this module does *not* exist.
        Defaults to ``None``, in which case a message is synthesized from this
        module name.

    Raises
    ----------
    ImportError
        If this module either does not exist *or* does exist but is
        unimportable (e.g., due to module-scoped side effects at importation
        time). To allow callers to transparently handle importation errors in
        the standard way, this function raises a standard exception rather than
        an application-specific exception (e.g.,
        :class:`betse.exceptions.BetseModuleException`).

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
    ``True`` only if the module with the passed fully-qualified name **exists**
    (i.e., is importable under the active Python interpreter).

    Caveats
    ----------
    In common edge cases, **this function may import all parent modules of this
    module as well as this module itself as a side effect.** Specifically, if
    this module is:

    * A **submodule** (i.e., contains a ``.`` character), this function
      necessarily imports *all* parent modules of this module.
    * *Not* importable via standard mechanisms (e.g., the OS X-specific
      :mod:`PyObjCTools` package), this function may import this module itself.

    Parameters
    ----------
    module_name : str
        Fully-qualified name of the module to be tested.

    Returns
    ----------
    bool
        ``True`` only if this module exists.
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
    # Since this function only returns a boolean value, the above non-boolean
    # values and exceptions are converted into a simple boolean.
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
    # This module may exist even if this exception is raised (e.g., for modules
    # dynamically defined at runtime rather than loaded from external files).
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

# ....................{ GETTERS                           }....................
@type_check
def get_parent_module_name_or_none(module_name: str) -> StrOrNoneTypes:
    '''
    Name of the parent module of the child module with the passed name if this
    module has a parent *or* ``None`` otherwise (i.e., if this is a top-level
    module).

    Parameters
    ----------
    module_name : str
        Fully-qualified name of the module to be munged.

    Returns
    ----------
    StrOrNoneTypes
        Either:

        * If this is a top-level module, ``None``.
        * Else, the fully-qualified name of the parent module of this module.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.text.string import strs

    # Return the prefix of this module name preceding the last "." delimiter.
    return strs.get_prefix_or_none(
        text=module_name, anchor='.', is_first=False)

# ....................{ IMPORTERS                         }....................
@type_check
def import_module(
    module_name: str, exception_message: StrOrNoneTypes = None) -> ModuleType:
    '''
    Dynamically import and return the module, package, or C extension with the
    passed fully-qualified name if importable *or* raise an exception with the
    passed message otherwise.

    Parameters
    ----------
    module_name : str
        Fully-qualified name of the module to be validated.
    exception_message : StrOrNoneTypes
        Message of the exception to be raised if this module does *not* exist.
        Defaults to ``None``, in which case a message is synthesized from this
        module name.

    Raises
    ----------
    ImportError
        If this module either does not exist *or* does exist but is
        unimportable (e.g., due to module-scoped side effects at importation
        time).
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

# ....................{ MAKERS                            }....................
@type_check
def make_module(
    module_name: str,
    module_attr_name_to_value: MappingOrNoneTypes = None,
    module_doc: StrOrNoneTypes = None,
    is_importable: bool = False,
) -> ModuleType:
    '''
    Dynamically create and return a new module with the passed name and
    optional passed docstring and attributes, optionally added to the standard
    :attr:`sys.modules` dictionary for external importation elsewhere.

    Caveats
    ----------
    **This function currently creates only non-package modules.** Although this
    function could be probably generalized to support dynamic creation of
    packages as well, packages impose certain real-world constraints *not*
    imposed by non-package modules (e.g., loading package submodules).

    Parameters
    ----------
    module_name : str
        Fully-qualified name of the module to be created.
    module_attr_name_to_value : MappingOrNoneTypes
        Dictionary mapping from the name to value of each module-scoped
        attribute to be declared in this module. Defaults to ``None``, in which
        case this module is empty (i.e., contains *no* attributes) by default.
    module_doc: StrOrNoneTypes
        **Docstring** (i.e., human-readable documentation in reStructuredText
        (reST) format) to be associated with this module. Defaults to ``None``,
        in which case this module is undocumented by default.
    is_importable : bool
        ``True`` only if external callers are allowed to trivially import this
        module via this module name elsewhere in this codebase. Note that, in
        this case, *all* parent packages of this module must already exist.
        Defaults to ``False`` for safety, in which case the returned object is
        the *only* initial reference to this module.

    Raises
    ----------
    BetseModuleException
        If a module with this name already exists.
    ImportError
        If the ``is_importable`` parameter is ``True`` *and* one or more parent
        packages of this module do *not* already exist.

    See Also
    ----------
    https://stackoverflow.com/questions/2931950/dynamic-module-creation
        StackOverflow *question* strongly inspiring this implementation. Note
        that, against all expectations, this question is substantially more
        informative than its answers.
    '''

    # If this module already exists, raise an exception.
    die_if_module(module_name)

    # Name of the parent module of this module if any *OR* None otherwise.
    parent_module_name = get_parent_module_name_or_none(module_name)

    # If this module is to be importable *AND* has a parent, raise an exception
    # if this parent does *NOT* already exist. By virtue of import mechanics,
    # this suffices to also recursively validate that *ALL* transitive parents
    # of this parent are also importable.
    if is_importable and parent_module_name is not None:
        die_unless_module(parent_module_name)

    # Module specification (i.e., metadata object describing this module).
    #
    # Note that module specifications are typically created and consumed by
    # low-level import machinery. In this case, however, we are creating a
    # module specification for the sole purpose of directly passing that same
    # module specification to the importlib.util.module_from_spec() function.
    module_spec = ModuleSpec(
        # Fully-qualified module name.
        module_name,

        # Module loader. Since this function dynamically creates this module,
        # this module is *NOT* loadable.
        #
        # Note that this parameter should technically only be "None" for
        # namespace packages. Unfortunately, there appears to exist no standard
        # loader class that effectively reduces to a noop. Ergo, "None" it is.
        None,
    )

    # Module defined by this specification to be returned.
    #
    # Note that modules may also be created by directly instantiating the
    # low-level "types.ModuleType" class but that doing so is officially
    # discouraged as "...spec is used to set as many import-controlled
    # attributes on the module as possible."
    module = importlib_util.module_from_spec(module_spec)

    # Set the name of the parent module of this module if any *OR* "None"
    # otherwise. Note that this special attribute could technically also be set
    # by passing the positional "parent" argument to the ModuleSpec.__init__()
    # method called above, but that doing so is complicated by:
    #
    # * That method's refusal to accept a keyword "parent" argument.
    # * The "ModuleSpec" class' refusal to allow callers to externally set the
    #   "parent" instance variable.
    module.__package__ = parent_module_name

    # If this module is to be prepopulated with attributes, do so.
    if module_attr_name_to_value is not None:
        module.__dict__.update(module_attr_name_to_value)

    # If this module is to be importable, register this module with the
    # standard "sys.modules" dictionary.
    if is_importable:
        sys.modules[module_name] = module

    # Return this module.
    return module
